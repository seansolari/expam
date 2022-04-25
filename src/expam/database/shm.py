from multiprocessing.shared_memory import SharedMemory
from os import ftruncate


class ResizableSharedMemory(SharedMemory):
    def __init__(self, name=None, create=False, size=0):
        super(ResizableSharedMemory, self).__init__(name, create, size)

    @property
    def size(self):
        """Size in bytes."""
        return self._size

    @size.setter
    def size(self, val):
        self._size = val

    def resize(self, newsize):
        fd = self._fd

        try:
            ftruncate(fd, newsize)
            return 1, newsize
        except OSError as e:
            return 0, e

    def remap(self, newsize):
        if newsize != self.size:
            del self._buf  # Clear buffer exposure.
            self._mmap.resize(newsize)  # Remap memory.
            self._buf = memoryview(self._mmap)  # Re-expose buffer.
            # Update size.
            self.size = newsize

