from numpy import *
import numpy as np


def f(x, u):
    return 5 * (x - u) + 1


def runge_kutta(f, x0, u0, h, n):
    x = x0
    u = u0
    for i in range(int(np.floor(n))):
        k1 = h * f(x, u)
        k2 = h * f(x + h / 2, u + k1 / 2)
        k3 = h * f(x + h / 2, u + k2 / 2)
        k4 = h * f(x + h, u + k3)
        u = u + (k1 + 2 * k2 + 2 * k3 + k4) / 6
        x = x + h
    return u


for pas in range(0, 5):
    h = 1 / (2**pas)
    print("h = ", h)
    print("u(4) = ", runge_kutta(f, 0, 1, h, 4 / h))
