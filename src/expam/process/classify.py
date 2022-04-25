from multiprocessing import shared_memory
import os
import queue
import re
import numpy as np
from expam.database import CHUNK_SIZE, EXPAM_TEMP_EXT, NULL_VALUE, TIMEOUT, DataTypeConfig, expam_dtypes
from expam.ext.kmers import get_raw_kmers
from expam.ext.map import classify
from expam.process import COMMAND_CODES
from expam.process.jobworker import JobWorker


class ClassifyWorker(JobWorker):
    def __init__(
        self, k, out_dir, temp_dir, kmer_db, lca_matrix, n_readers, logging_dir, alpha, status_index,
        statuses, child_statuses, timeout=TIMEOUT, phases=None, dtype_config: DataTypeConfig = expam_dtypes
    ):
        super(ClassifyWorker, self).__init__(logging_dir, status_index, statuses, child_statuses, timeout=timeout,
                                             phases=phases)
        self.job_name = "classify"
        self.n_readers = n_readers

        # Configuration parameters.
        self.k = k
        self.alpha = alpha
        self.chunks = (k // CHUNK_SIZE) + (k % CHUNK_SIZE != 0)

        self.out_dir = out_dir
        self.temp_dir = temp_dir

        self._read_array_length = 150
        self.read_array = np.ndarray((self._read_array_length, self.chunks + 1), dtype=dtype_config.keys_dtype)

        self.RAW = "_%s.reads_%s" % ("%d", EXPAM_TEMP_EXT)

        # Load the database.
        self.dtypes: DataTypeConfig = dtype_config
        self.keys_shm = shared_memory.SharedMemory(
            name=kmer_db.keys_shm_name,
            create=False
        )
        self.keys = np.ndarray(
            shape=kmer_db.keys_meta["shape"],
            dtype=kmer_db.keys_meta["dtype"],
            buffer=self.keys_shm.buf
        )
        self.values_shm = shared_memory.SharedMemory(
            name=kmer_db.values_shm_name,
            create=False
        )
        self.values = np.ndarray(
            shape=kmer_db.values_meta["shape"],
            dtype=kmer_db.values_meta["dtype"],
            buffer=self.values_shm.buf
        )

        self.current_file_data = [("", []) for _ in range(n_readers)]

        # LCA Processing.
        self.lca_matrix = lca_matrix

        self.init_methods_from_dict(
            {
                "classify": {
                    "getter": None,
                    "processor": self.classify_read,
                    "transition": self.save_final_data,
                }
            }
        )

    def classify_read(self, read_data):
        new_file_name, read_id, read, reader_id = read_data

        # Continue to get reads without checking for errors; only check for errors when there is a break in workflow.
        while True:
            read_id = read_id.decode("utf8")

            max_n_kmers = len(read) - self.k + 1
            if max_n_kmers > self.read_array.shape[0]:
                self.read_array = np.ndarray((max_n_kmers, self.chunks + 1), dtype=self.dtypes.keys_dtype)

            if new_file_name != self.current_file_data[reader_id][0]:
                if len(self.current_file_data[reader_id][1]) > 0:
                    self.save_temp_data(reader_id)

                self.current_file_data[reader_id] = (new_file_name, [])

            result, best_clade, read_length, result_string = self.process_read(read)
            self.current_file_data[reader_id][1].append(
                [result, read_id, "p" + str(best_clade), str(read_length), result_string]
            )

            # Get a head start on any available work.
            try:
                command_code, data = self.work_queue.get(block=True, timeout=self._timeout)

                if command_code != COMMAND_CODES["ACTIVE"]:
                    self.work_queue.put((command_code, data))
                    break

                new_file_name, read_id, read, reader_id = data
            except queue.Empty:
                break

    def save_temp_data(self, reader_id):
        def format_result(result_list):
            return "\n".join(["\t".join(item) for item in result_list])

        def base_name(file_url):
            return re.match(r"(\S+)\.(?:fa|fq|fasta|fastq)(?:\.tar)*(?:\.gz)*$", file_url).group(1)

        file_name, file_data = self.current_file_data[reader_id]

        my_file_name = os.path.join(self.temp_dir, base_name(file_name) + (self.RAW % self.pid))
        with open(my_file_name, "a+") as f:
            f.write(format_result(file_data))

    def process_read(self, read):
        read_length = len(read)
        n_kmers = get_raw_kmers(read, self.k, self.read_array)

        # Control for reads that don't have long enough good regions for kmers to be extracted.
        if n_kmers == 0:
            return "U", 0, read_length, ""

        # Map the kmers.
        cls, common_clade, kmer_string = classify(self.read_array[:n_kmers], self.keys, self.values,
                                                  self.lca_matrix, NULL_VALUE, self.alpha)

        if cls == -1:
            result = "U"

        elif cls == 1:
            result = "S"

        else:
            result = "C"

        return result, common_clade, read_length, kmer_string

    def save_final_data(self, *args):
        for reader_id in range(self.n_readers):
            if len(self.current_file_data[reader_id][1]) > 0:
                self.save_temp_data(reader_id)

        self.keys_shm.close()
        self.values_shm.close()
