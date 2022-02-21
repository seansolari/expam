import random

import expam.c.sets
import numpy as np

def split_kmers(arr):
    mask = np.random.randint(2, size=arr.shape[0], dtype=np.uint64)

    new = np.compress(1-mask, arr)
    arr = np.compress(mask, arr)

    return arr, new

def main():
    n = random.randint(1000, 10000)
    print("Generating array of length %d!" % n)

    one = np.arange(n, dtype=np.uint64)

    one, two = split_kmers(one)
    two, three = split_kmers(two)
    three, four = split_kmers(three)
    one, five = split_kmers(one)
    five, six = split_kmers(five)

    base = np.zeros(n, dtype=np.uint64)

    base = base.reshape((base.shape[0], 1))
    one = one.reshape((one.shape[0], 1))
    two = two.reshape((two.shape[0], 1))
    three = three.reshape((three.shape[0], 1))
    four = four.reshape((four.shape[0], 1))
    five = five.reshape((five.shape[0], 1))
    six = six.reshape((six.shape[0], 1))

    print("Filling array...")
    expam.c.sets.filler(256, 0, [base, one, two, three, four, five, six])

    print("Filled properly?", np.array_equal(base, np.arange(n).reshape((n, 1))))

if __name__ == '__main__':
    main()
