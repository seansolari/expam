#!/usr/bin/env python3
from functools import partial
from multiprocessing import Pool


def f(a, b, c):
    return (a, b, c), a + b + c

if __name__ == "__main__":
    data = [(1, 2), (3, 4), (5, 6)]

    def combine(t, s):
        print(", ".join(str(v) for v in t), "-->", s)
    
    def com(v):
        print(v)

    with Pool(2) as mp:
        g = partial(f, c = 0.5)
        result = mp.starmap_async(g, data, callback=com)
