import io
import re
from expam.database import TIMEOUT
from expam.ext.kmers import bytes_intersect
from expam.process import COMMAND_CODES
from expam.process.jobworker import JobWorker
from expam.sequences import format_name, get_opener


class ReadWorker(JobWorker):
    def __init__(
        self, paired_mode, n_classifiers, logging_dir, status_index, statuses, child_statuses,
        timeout=TIMEOUT, phases=None
    ):
        super(ReadWorker, self).__init__(
            logging_dir, status_index, statuses, child_statuses,
            timeout=timeout, phases=phases
        )
        self.job_name = "reader"

        self.paired_mode = paired_mode
        self.n_classifiers = n_classifiers
        self._c = 0  # Which classifier do we next send data to?

        self.n_streams = 1 if not paired_mode else 2

        self.init_methods_from_dict(
            {
                "classify": {
                    "getter": None,
                    "processor": self.load_reads,
                    "transition": None,
                }
            }
        )

    def load_reads(self, file_urls):
        print("Loading reads from %s..." % ", ".join(file_urls))

        file_names, file_modes, file_openers, file_types = self.init_openers(*file_urls)

        # Open each file.
        _file_handles = [
            opener(file_url, mode)
            for opener, file_url, mode in zip(file_openers, file_urls, file_modes)
        ]

        read_stream = io.BytesIO()

        read_ids = []
        read_seq = None

        while read_ids or read_seq is None:  # read_seq is None is initial case.
            for i in range(self.n_streams):
                file_type = file_types[i]

                if file_type == "fasta":
                    next_read_id = self.scan_fa_to_next_read(_file_handles[i], read_stream)

                else:  # File type is "fastq".
                    next_read_id = self.scan_fq_to_next_read(_file_handles[i], read_stream)

                if next_read_id is not None:
                    if i == 0:
                        read_ids.append(next_read_id)
                    else:
                        read_ids[-1] = bytes_intersect(read_ids[-1], next_read_id)

            read_seq = read_stream.getvalue()

            if len(read_seq) > self.n_streams:
                read_id = read_ids.pop(0)

                self.send_read(min(file_names), read_id, read_seq)

            read_stream.seek(0)
            read_stream.truncate()

        for i in range(self.n_streams):
            # Close file handles.
            _file_handles[i].close()

        # Close data stream.
        del read_stream

    def send_read(self, file_name, read_id, read):
        self.tell_child(
            cid=self._c,
            data=(
                COMMAND_CODES["ACTIVE"],
                (file_name, read_id, read, self.child_id)
            )
        )

        self._c = (self._c + 1) % self.n_classifiers

    def init_openers(self, *file_urls):
        file_names, file_modes, file_openers, file_types = [], [], [], []

        for file_url in file_urls:
            file_name = format_name(file_url)
            file_name, (mode, opener) = get_opener(file_name)

            if self.check_suffix(file_url, ".fasta", ".fa"):
                file_type = "fasta"

            elif self.check_suffix(file_url, ".fastq", ".fq"):
                file_type = "fastq"

            else:
                raise ValueError("Unknown sequence format for %s!" % file_url)

            file_names.append(file_name)
            file_modes.append(mode)
            file_openers.append(opener)
            file_types.append(file_type)

        return file_names, file_modes, file_openers, file_types

    @staticmethod
    def check_suffix(string, *suffixes):
        for suffix in suffixes:
            if re.findall(r"(%s)(?:.tar)*(?:.gz)*$" % suffix, string):
                return True

        return False

    @staticmethod
    def scan_fa_to_next_read(file, read_stream):
        next_read_id = None

        while True:
            line = file.readline()

            if line == b'':
                read_stream.write(b'N')
                break  # EOF.

            if b'>' in line:
                next_read_id = re.match(b'^>(\S+)\s', line).group(1)

                # End my contribution to this read.
                read_stream.write(b'N')
                break

            else:
                read_stream.write(line.strip())

        return next_read_id

    @staticmethod
    def scan_fq_to_next_read(file, read_stream):
        next_read_id = None

        while True:
            line = file.readline()

            if line == b'':
                read_stream.write(b'N')
                break

            if b'@' in line:
                next_read_id = re.match(b'^@(.*)\s', line).group(1)

                read_stream.write(b'N')
                break

            else:
                read_stream.write(line.strip())

                # Skip next two lines.
                for _ in range(2):
                    file.readline()

        return next_read_id

