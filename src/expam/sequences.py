import gzip
import io
import os
import random
import re
from expam.ext.kmers import reverse_complement_combine
from expam.logger import timeit

COMPRESSION_EXTNS = ['.tar.gz', '.tar', '.gz']
DEFAULT_MODE = "rb"
DEFAULT_OPENER = open
COMP_PARSE = {
    ".tar.gz": {"mode": "rb", "opener": gzip.open},
    ".gz": {"mode": "rb", "opener": gzip.open}
}


def get_opener(file_name):
    for suffix in COMP_PARSE:

        if check_suffix(file_name, suffix):
            file_name = file_name.replace(suffix, '')
            mode, opener = COMP_PARSE[suffix]["mode"], COMP_PARSE[suffix]["opener"]

            break

    else:
        mode, opener = DEFAULT_MODE, DEFAULT_OPENER

    return file_name, (mode, opener)


def check_suffix(string, sfx):
    string = string.strip()

    if string[-len(sfx):] == sfx:
        return True
    return False


def format_name(name, remove_comp=False):
    name = os.path.basename(name)
    # Must not be a digit. If it is, add 'ref' to the start.
    if name.isdigit():
        name = "ref" + name

    # Extract sequence name.
    if ".fna" in name:
        name = name.replace(".fna", "")
    elif ".faa" in name:
        name = name.replace(".faa", "")

    if remove_comp:
        for suffix in COMP_PARSE:
            if check_suffix(name, suffix):
                name = name.replace(suffix, '')

                break

    return name


def extract_id(metadata):
    # Check for kraken format.
    if "kraken" in metadata:
        start = metadata.rfind("|") + 1
    else:
        # Standard fasta format.
        if ">" in metadata:
            start = metadata.find(">") + 1
        elif "@" in metadata:
            start = metadata.find("@") + 1

    end = metadata.find(" ")
    if not end:  # Case where no metadata is provided.
        end = len(metadata)

    # Contig name, including potential contig version number.
    contig_id = metadata[start:end]
    # Concatenate multiple contig versions.
    version_start = contig_id.rfind(".")
    if version_start:
        contig_id = contig_id[:version_start]
    return contig_id


@timeit
def iterate_sequence(url):
    with open(url, "r") as f:
        seqID, contig = "", ""
        for sequence in f:
            # Strip newlines.
            if '\n' in sequence:
                sequence = sequence.strip('\n')
            if ">" in sequence:
                if contig != "":
                    # Return the previous sequence, if there is one.
                    yield (seqID, contig)
                # Find seqID.
                if "|" in sequence:
                    # >x|seqID|... format.
                    seqID = sequence.split("|")[1]
                else:
                    # >seqID ... format.
                    seqID = sequence.split(" ")[0]
            else:
                # Combine sequence.
                contig += sequence
        # End case: yeild last contig.
        if contig != "":
            yield seqID, contig


@timeit
def get_sequence(url):
    """
    Extract the sequence from FASTA source file @ url.

    :param url: String.
    :return: String.
    """
    sequence = io.BytesIO()
    _, (mode, opener) = get_opener(url)

    with opener(url, mode) as f:
        for line in f:
            line = line.strip()
            # Skip information line at start of FASTA format.
            if b">" in line:
                sequence.write(b"N")

            else:
                sequence.write(line)

    return sequence.getvalue()

class MaskStore:
    def __init__(self):
        self.store = {}

    def __setitem__(self, key, value):
        # Place value in store.
        self.store[key] = value

    def __getitem__(self, key):
        return self.store[key]

    def __delitem__(self, key):
        del self.store[key]


class SequenceStore:
    def __init__(self, k, out_dir):
        self._store = {}
        self.stream = io.BytesIO()
        self.k = k
        self.out_dir = out_dir

        self.accession_ids = {}

    def add_sequence(self, file_path):
        # Check it is a valid file.
        if not os.path.isfile(file_path):
            raise ValueError("Directory %s is not a file!" % file_path)

        print("Extracting sequences from %s..." % file_path)

        # Determine opening mode based on compression of sequence file.
        file_name = format_name(file_path)
        file_name, (mode, opener) = get_opener(file_name)

        with opener(file_path, mode=mode) as f:
            # Extract the accession id.
            metadata = f.readline()
            assert b'>' in metadata

            accession_id, tax_id = self.extract_accession_id(metadata)

            if accession_id is None and tax_id is None:
                raise ValueError("Couldn't identify sequence %s!" % file_name)

            for line in f:
                line = line.strip()
                # Skip information line at start of FASTA format.
                if b'>' in line:
                    self.stream.write(b'N')
                else:
                    self.stream.write(line)

        # Convert to numeric sequence and put in store.
        sequence_view = self.stream.getvalue()
        self._store[file_name] = sequence_view

        # Save accession ID for this file.
        self.accession_ids[file_name] = [accession_id, tax_id]

        # Clear stream.
        self.stream.seek(0)
        self.stream.truncate()

        return file_name, len(self._store[file_name])

    def __getitem__(self, key):
        # If it is a branch, this is a valid request,
        # but no information from the Sequence Store is required.
        if key.isdigit():
            return None
        else:
            return self._store.__getitem__(key)

    def __setitem__(self, key, value):
        self._store.__setitem__(key, value)

    def __delitem__(self, key):
        self._store.__delitem__(key)

    def keys(self):
        return self._store.keys()

    @staticmethod
    def extract_accession_id(string):
        def try_convert(byte_string):
            try:
                return byte_string.decode("utf8")

            except AttributeError:
                return byte_string

        accn_id = re.findall(b'(?:^\s*\>)(.*?)(?:\s+|\|)', string)
        tax_id = re.findall(b'(?:\|kraken:taxid\||\|taxid\|)(\d+)(?:\s+|\|)', string)

        accn_id = None if not accn_id else try_convert(accn_id[0])
        tax_id = None if not tax_id else try_convert(tax_id[0])

        return accn_id, tax_id

    def save_accession_ids(self, suffix=""):
        id_data = "\n".join([
            "%s,%s,%s" % (name, keys[0], keys[1])
            for name, keys in self.accession_ids.items()
        ])

        with open(os.path.join(self.out_dir, "phylogeny", "accession_ids_%s.csv" % suffix), "w") as f:
            f.write(id_data)


def make_reads(in_url, l, n, out_url, e=0.0):
    """
    Make n reads from sequence at url, from a expam maker class.
    """
    maker = Reads(l, out_url)
    maker.make(in_url, n, error=e)


class Reads:
    def __init__(self, l, out_dir):
        # Length of reads.
        self.l = l
        # Directory to save reads.
        self.out_dir = out_dir

    def make(self, f_url, n_reads, error=0.0):
        """
        Produce a number of FASTQ reads from a sequence.

        :param error: Float - Rate that errors will be put into reads (error ~ reads with errors / reads).
        :param f_url: String - Url of the sequence.
        :param n_reads: Int - number of reads to produce.
        """
        if error > 1:
            if error > 100:
                raise ValueError("Error rate must be either < 1 or < 100!")

            error /= 100

        nt_codes = {0: b"A", 1: b"C", 2: b"G", 3: b"T"}
        bad_quality = bytes("N" * 75, encoding="utf8")

        # Check formatting of sequence.
        self.check_fna(f_url)
        # Extract sequence string from f_url.
        sequence = get_sequence(f_url)

        # Extract the file name from the path.
        forward_name, reverse_name = self.get_file_name(f_url)

        # Generate the reads.
        i = 0
        forward_read_stream = io.BytesIO()
        reverse_read_stream = io.BytesIO()

        while i < n_reads:

            # Generate a start index for the read.
            index = random.randint(0, len(sequence) - self.l)
            # Get read.
            read = sequence[index: index + self.l]

            # Check quality.
            if bad_quality in read:
                # This is deemed poor quality. Try again.
                continue

            rid = ("R" + str(abs(hash(str(index))))).encode('utf8')
            meta_line = b">%b Read from %b.\n" % (rid, f_url.encode('utf8'))

            # Generate error.
            has_error = random.randint(0, 100) / 100 < error

            if has_error:
                sub_location = random.randint(0, len(read) - 1)
                nt = nt_codes[random.randint(0, 3)]

                read_string = b"%b%b%b" % (read[:sub_location], nt, read[sub_location + 1:])

            else:
                read_string = b"%b" % read

            forward_reverse_string = reverse_complement_combine(read_string, read_string)
            l = len(read_string)

            forward_string = forward_reverse_string[:l]
            reverse_string = forward_reverse_string[-l:]

            forward_read_stream.write(meta_line)
            forward_read_stream.write(forward_string)
            forward_read_stream.write(b'\n')

            reverse_read_stream.write(meta_line)
            reverse_read_stream.write(reverse_string)
            reverse_read_stream.write(b'\n')

            i += 1

        # Write to file.
        with open(os.path.join(self.out_dir, forward_name), "wb") as f:
            f.write(forward_read_stream.getvalue())

        forward_read_stream.seek(0)
        forward_read_stream.truncate()

        with open(os.path.join(self.out_dir, reverse_name), "wb") as f:
            f.write(reverse_read_stream.getvalue())

        reverse_read_stream.seek(0)
        reverse_read_stream.truncate()

    @staticmethod
    def check_fna(url):
        if not ".fna" in url:
            raise TypeError("Sequence type must be .fna.")

    @staticmethod
    def get_file_name(url):
        path = os.path.normpath(url)
        path = path.split(os.sep)

        forward_name = path[-1].rstrip(".fna") + "_1.fa"
        reverse_name = path[-1].rstrip(".fna") + "_2.fa"

        return forward_name, reverse_name
