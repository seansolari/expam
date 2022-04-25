import multiprocessing
import queue
import numpy as np
from expam.database import CHUNK_SIZE, TIMEOUT, DataTypeConfig, expam_dtypes
from expam.ext.kmers import get_raw_kmers, remove_duplicates, to_32bit, import_mask, map_kmers
from expam.ext.map import binarysearch
from expam.logger import log_timer
from expam.process import COMMAND_CODES
from expam.process.jobworker import JobWorker
from expam.sequences import MaskStore, SequenceStore


class ExtractWorker(JobWorker):
    def __init__(
        self, k, shm_kmer_params, kmer_ranges, node_map, out_dir, logging_dir,
        status_index, statuses, child_statuses, timeout=TIMEOUT,
        phases=None, dtype_config: DataTypeConfig = expam_dtypes
    ):
        super(ExtractWorker, self).__init__(logging_dir, status_index, statuses, child_statuses,
                                            timeout=timeout, phases=phases)

        self.job_name = "extract"
        self.k = k
        self.kmer_ranges = kmer_ranges

        self.chunks = (k // CHUNK_SIZE) + (k % CHUNK_SIZE != 0)
        self.keys_dtype = dtype_config.keys_dtype
        self.d_type = np.dtype(
            [
                ('index', self.keys_dtype),
                *(
                    ('k%d' % chunk, self.keys_dtype)
                    for chunk in range(self.chunks)
                )
            ]
        )

        # For mapping LCAs.
        self.node_map = node_map

        # Hold the sequences and masks that I have extracted.
        self.sequence_store = SequenceStore(k, out_dir)
        self.mask_store = MaskStore()

        shm_kmer_name, shm_kmer_shape = shm_kmer_params
        allocation_shape = (shm_kmer_shape[0], shm_kmer_shape[1] + 1)

        # Long-term arrays used for extraction.
        self.allocation = np.zeros(
            shape=allocation_shape,  # Account for index col.
            dtype=self.keys_dtype,
            order="C",  # Row-major storage (consecutive k-mer access).
        )

        # Long-term array for sending kmers to children.
        self.out_shape = shm_kmer_shape
        self.out_name = shm_kmer_name
        self.out_allocation = multiprocessing.shared_memory.SharedMemory(shm_kmer_name)
        self.out_array = np.ndarray(
            shape=shm_kmer_shape,
            dtype=self.keys_dtype,
            buffer=self.out_allocation.buf
        )

        # Primary queue for getting work during mapping phase.
        self.map_queue = multiprocessing.Queue()

        # Initialise methods for execution.
        self.init_methods_from_dict(
            {
                "import": {
                    "getter": None,
                    "processor": self.import_sequence,
                    "transition": self.prepare_map,
                },
                "map": {
                    "getter": self.get_my_sequence,
                    "processor": self.send_for_mapping,
                    "transition": self.empty_queue,
                },
                "shutdown": self.close_shm
            }
        )

    @log_timer
    def import_sequence(self, sequence_url):
        def resize_base(size):
            nonlocal base

            base = base[:size]
            update_view()

        def new_view():
            nonlocal base
            return base.ravel().view(dtype=self.d_type)

        def update_view():
            nonlocal base, view
            view = base.ravel().view(dtype=self.d_type)

        # Collect raw sequence.
        self.log("Importing sequence at %s." % sequence_url)
        seq_id, seq_length = self.sequence_store.add_sequence(sequence_url)

        # Get view on required segment of allocation.
        base = self.allocation[:seq_length]
        view = new_view()

        # Create draft mask and get corresponding kmers.
        new_size = get_raw_kmers(
            sequence=self.sequence_store[seq_id],
            k=self.k,
            arr=base
        )
        resize_base(new_size)

        # Sort by kmers, which also sorts the indices
        view.sort(
            order=['k%d' % chunk for chunk in range(self.chunks)],
            kind="stable"
        )
        self.log("Removing duplicates.")
        new_size = remove_duplicates(base)
        resize_base(new_size)

        # Save mask as 32-bit array.
        self.mask_store[seq_id] = to_32bit(view['index'])
        self.send_kmers(
            ranges=self.kmer_ranges,
            kmers=base[:, 1:]
        )

    @staticmethod
    def clear_temp_arrays(arr):
        arr[:] = 0

    @log_timer
    def send_kmers(self, ranges, kmers, node_id=None):
        target = np.ndarray(shape=(kmers.shape[1]), dtype=kmers.dtype)
        target[:] = ranges[0][0]

        total = 0
        n = len(ranges) - 1
        lower_ind = upper_ind = self.get_index(target, kmers, greater=False)

        # Wait for previous data to be sent.
        while self.loop_condition(not np.all(self.out_array == 0)):
            self.timeout()

        # Place kmers in out_array.
        self.out_array[:kmers.shape[0], :] = kmers[:]

        # Separate data between children.
        for i, (lower_bound, upper_bound) in ranges.items():
            target[:] = upper_bound
            upper_ind = self.get_index(target, kmers)

            # Upper bound can only attain equality for maximum hash value.
            if i == n:
                while (upper_ind < kmers.shape[0]) and self.kmers_equal(target, kmers[upper_ind]):
                    upper_ind += 1

            # Don't send empty work.
            if lower_ind != upper_ind:

                message_data = (
                    (lower_ind, upper_ind),
                    (self.out_name, self.out_shape),
                )
                if node_id is not None:
                    message_data += (node_id,)

                self.tell_child(
                    cid=i,
                    data=(
                        COMMAND_CODES["ACTIVE"],
                        message_data
                    )
                )

                total += upper_ind - lower_ind

            # Beginning next search.
            lower_ind = upper_ind

        self.log("%d / %d kmers sent!" % (total, kmers.shape[0]))

    def get_index(self, target, kmers, greater=True):
        def _check_ind():
            nonlocal ind, kmers

            if (greater and ind == kmers.shape[0]) or ((not greater) and ind == 0):
                return True

            return False

        # Find closest element in list. Note that this element
        # may be less than the bound.
        _, ind = binarysearch(target, kmers)

        if _check_ind():
            return int(ind)

        while (not _check_ind()) and self.kmers_compare(kmers[ind], target, greater):
            if greater:
                ind += 1

            else:
                ind -= 1

        return int(ind)

    @staticmethod
    def kmers_equal(one, two):
        return np.all(one == two)

    @staticmethod
    def kmers_compare(one, two, greater):
        for i in range(one.shape[0]):
            if one[i] < two[i]:
                return greater

            elif one[i] > two[i]:
                return not greater

            elif i + 1 == one.shape[0]:
                return not greater

    def prepare_map(self, *args):
        # Save accession ids from my sequences.
        self.sequence_store.save_accession_ids(suffix=str(self._status_id))

        # Prepare mapping jobs.
        self.load_queue()

    def load_queue(self):
        for seq_id in self.sequence_store.keys():
            self.map_queue.put(
                (
                    COMMAND_CODES["ACTIVE"],
                    seq_id
                )
            )

    def get_my_sequence(self, *args):

        # Try to get local work first.
        try:
            work = self.map_queue.get(block=True, timeout=self._timeout)

            return work

        except queue.Empty:
            pass

        # Now look for commands from parent.
        try:
            work = self.work_queue.get(block=True, timeout=self._timeout)

            return work

        except queue.Empty:
            return None

    def send_for_mapping(self, seq_id):
        node_id = self.node_map[seq_id]

        sequence = self.sequence_store[seq_id]
        mask = self.mask_store[seq_id]
        n_kmers = mask.shape[0]

        # Get view on required segment of allocation.
        base = self.allocation[:n_kmers]
        view = base.ravel().view(dtype=np.dtype(
            [
                ('index', self.keys_dtype),
                *(
                    ('k%d' % chunk, self.keys_dtype)
                    for chunk in range(self.chunks)
                )
            ]
        ))

        # Place 32-bit mask in 64-bit array.
        import_mask(
            data=mask,
            dest=view['index']
        )
        map_kmers(
            sequence=sequence,
            k=self.k,
            arr=base
        )

        # Send these kmers to be mapped.
        self.send_kmers(
            ranges=self.kmer_ranges,
            kmers=base[:, 1:],
            node_id=node_id
        )

    def empty_queue(self, *args):
        try:
            _ = self.map_queue.get(block=True, timeout=self._timeout)

            raise ValueError("Mapping queue was not empty!")

        except queue.Empty:

            return

    def close_shm(self):
        try:
            self.out_allocation.unlink()
            self.out_allocation.close()

        except FileNotFoundError:
            pass

