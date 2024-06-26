# -------------------------------------------------------------------------
#
# PYTHON for DUMMIES
# Problème 3
#
# Script de test
#  Vincent Legat
#
# -------------------------------------------------------------------------
#
from time import perf_counter
from numpy import *
import matplotlib.pyplot as plt
import matplotlib

# ============================================================
# FONCTIONS A MODIFIER [begin]
#


def b(t, T, i, p):
    if p == 0:
        return (T[i] <= t) * (t < T[i + 1])
    u = 0.0 if T[i + p] == T[i] else (t - T[i]) / (T[i + p] - T[i]) * b(t, T, i, p - 1)
    u += (
        0.0
        if T[i + p + 1] == T[i + 1]
        else (T[i + p + 1] - t) / (T[i + p + 1] - T[i + 1]) * b(t, T, i + 1, p - 1)
    )
    return u


def bspline(X, Y, t):
    p = 3
    m = len(X)
    T = arange(-p, m + p + 1)

    n = len(T) - 1

    B = zeros((n - p, len(t)))
    for i in range(0, n - p):
        B[i, :] = b(t, T, i, p)
    X = [*X, *X[0:p]]
    Y = [*Y, *Y[0:p]]
    x = X @ B
    y = Y @ B
    return x, y


def deBoor_spline(X, Y, t):
    p = 3
    m = len(X)
    T = arange(-p, m + p + 1)
    i = searchsorted(T[1:], t)
    X = array([*X, *X[0:p]])
    Y = array([*Y, *Y[0:p]])

    d_x = zeros((p + 1, len(t)))
    d_y = zeros((p + 1, len(t)))
    for j in range(0, p + 1):
        d_x[j] = X[j + i - p]
        d_y[j] = Y[j + i - p]
    for r in range(1, p + 1):
        for j in range(p, r - 1, -1):
            alpha = (t - T[j + i - p]) / (T[j + 1 + i - r] - T[j + i - p])
            d_x[j] = (1.0 - alpha) * d_x[j - 1] + alpha * d_x[j]
            d_y[j] = (1.0 - alpha) * d_y[j - 1] + alpha * d_y[j]

    return d_x[p], d_y[p]


#
# FONCTIONS A MODIFIER [end]
# ============================================================
#
def main():
    # -1- Approximation d'un rectangle :-)
    #
    X = array([0, 3, 3, 0])
    Y = array([0, 0, 2, 2])
    t = linspace(0, len(X), len(X) * 1000 + 1)

    a = perf_counter()
    for i in range(1000):
        x, y = bspline(X, Y, t)
    print(perf_counter() - a)

    a = perf_counter()
    for i in range(1000):
        x, y = deBoor_spline(X, Y, t)
    print(perf_counter() - a)
    #
    # -2- Un joli dessin :-)
    #

    matplotlib.rcParams["toolbar"] = "None"
    plt.rcParams["figure.facecolor"] = "white"

    fig = plt.figure("Approximation avec des B-splines")
    plt.plot(X, Y, ".r", markersize=10)
    plt.plot([*X, X[0]], [*Y, Y[0]], "--r")
    plt.plot(x, y, "-b")
    plt.axis("equal")
    plt.axis("off")
    plt.show()


def b_spline_show():

    p = 3
    X = array([1, 1, 1, 1])
    Y = array([0, 0, 2, 2])
    t = linspace(-p, len(X) + p, len(X) * 1000 + 1)
    m = len(X)
    T = arange(-p, m + p + 1)
    n = len(T) - 1
    for p in range(4):
        B = zeros((n - p, len(t)))
        for i in range(0, n - p):

            B[i, :] = b(t, T, i, p)
            plt.plot(t, B[i, :])

        X = ones(n - p)
        plt.plot(t, X @ B)
        plt.show()


if __name__ == "__main__":
    b_spline_show()
