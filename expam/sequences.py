import io
import gzip
import os
import random
from os import listdir
from os.path import basename, isfile, join

from expam.c.extract import reverse_complement_combine
from expam.run import timeit

# Parsing of compressed files.
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


def ls(url, ext=""):
    null_ext = not len(ext)

    if isfile(url) and (url[-len(ext):] == ext or null_ext):
        return [url]
    else:
        return [
            f for f in listdir(url)
            if isfile(join(url, f)) and (f[-len(ext):] == ext or null_ext)
        ]


def format_name(name, remove_comp=False):
    name = basename(name)
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
