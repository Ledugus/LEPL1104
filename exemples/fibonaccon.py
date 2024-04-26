import numpy as np
from time import perf_counter


def fibonnacci(n):
    f = np.zeros(n, dtype=np.int64)
    f[0] = 1
    f[1] = 2
    for x in range(2, n):
        f[x] = f[x - 1] + f[x - 2]
    return f[-1]


def fibonnaccon_with_pre_alloc(n):

    f = [0] * n
    f[0] = 1
    f[1] = 2
    for x in range(2, n):
        f[x] = f[x - 1] + f[x - 2]
    return f[-1]


def fibonnaccon(n):
    f = [1, 2]
    for x in range(2, n):
        new = f[x - 1] + f[x - 2]
        f.append(new)
    return f[-1]


print(fibonnacci(5))
print(fibonnaccon(5))
n = 200000
a = perf_counter()
fibonnacci(n)
b = perf_counter()
print("fibonnacci : ", b - a)

a = perf_counter()
fibonnaccon(n)
b = perf_counter()
print("fibonnaccon : ", b - a)

a = perf_counter()
fibonnaccon_with_pre_alloc(n)
b = perf_counter()
print("fibonnaccon_with_pre_alloc : ", b - a)
