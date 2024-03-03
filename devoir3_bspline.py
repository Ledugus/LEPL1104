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


def deBoor(
    k,
    x,
    t,
    c,
):
    """Evaluates S(x).

    Arguments
    ---------
    k: Index of knot interval that contains x.
    x: Position.
    t: Array of knot positions, needs to be padded as described above.
    c: Array of control points.
    """
    p = 3
    d = [c[j + k - p] for j in range(0, p + 1)]

    for r in range(1, p + 1):
        for j in range(p, r - 1, -1):
            alpha = (x - t[j + k - p]) / (t[j + 1 + k - r] - t[j + k - p])
            d[j] = (1.0 - alpha) * d[j - 1] + alpha * d[j]

    return d[p]


def deBoor_spline(X, Y, t):
    p = 3
    m = len(X)
    T = arange(-p, m + p + 1)
    i = searchsorted(T[1:], t)
    X = [*X, *X[0:p]]

    Y = [*Y, *Y[0:p]]

    d = [c[j + k - p] for j in range(0, p + 1)]

    for r in range(1, p + 1):
        for j in range(p, r - 1, -1):
            alpha = (x - t[j + k - p]) / (t[j + 1 + k - r] - t[j + k - p])
            d[j] = (1.0 - alpha) * d[j - 1] + alpha * d[j]

    x = deBoor(i, t, T, X)

    y = deBoor(i, t, T, Y)

    return x, y


#
# FONCTIONS A MODIFIER [end]
# ============================================================


def main():
    #
    # -1- Approximation d'un rectangle :-)
    #

    X = [0, 3, 3, 0]
    Y = [0, 0, 2, 2]
    t = linspace(0, len(X), len(X) * 1000 + 1)
    a = perf_counter()
    x, y = bspline(X, Y, t)
    #
    # -2- Un joli dessin :-)
    #

    # matplotlib.rcParams["toolbar"] = "None"
    plt.rcParams["figure.facecolor"] = "white"

    fig = plt.figure("Approximation avec des B-splines")
    plt.plot(X, Y, ".r", markersize=10)
    plt.plot([*X, X[0]], [*Y, Y[0]], "--r")
    plt.plot(x, y, "-b")
    plt.axis("equal")
    plt.axis("off")
    plt.show()


if __name__ == "__main__":
    main()
