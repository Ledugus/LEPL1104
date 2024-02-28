#
# PYTHON for DUMMIES 23-24
# Problème 2
#
# Splines cubiques périodiques
#
#  Vincent Legat
#
# -------------------------------------------------------------------------
#


from numpy import *
import numpy as np
from numpy.linalg import solve


# ============================================================
# FONCTIONS A MODIFIER [begin]
#
#


def spline(x, h, U):
    n = np.size(U)
    X = np.arange(0, n + 1) * h
    i = np.searchsorted(X[1:], x)

    d = np.full(n, 4)
    A = np.diag(d) + np.diag(np.ones(n - 1), -1) + np.diag(np.ones(n - 1), 1)
    A[-1, 0] = 1
    A[0, -1] = 1

    b = U[range(-1, n - 1)] - 2 * U[0:n] + U[np.append(np.arange(1, n), [0])]
    dd_u = solve(A, b) / ((h**2) / 6)
    U = np.append(U, U[0])
    dd_u = np.append(dd_u, dd_u[0])

    n += 1
    a = dd_u[: n - 1] / (6 * h)
    b = dd_u[1:n] / (6 * h)
    c = U[: n - 1] / h - dd_u[: n - 1] * h / 6
    d = U[1:n] / h - dd_u[1:n] * h / 6

    r = X[i + 1] - x
    s = x - X[i]

    return r * (a[i] * r * r + c[i]) + s * (b[i] * s * s + d[i])


#
# FONCTIONS A MODIFIER [end]
# ============================================================
#
# -1- Interpolation d'un cercle :-)
#


def main():
    from matplotlib import pyplot as plt

    plt.rcParams["toolbar"] = "None"
    plt.rcParams["figure.facecolor"] = "lavender"

    n = 4
    h = 3 * np.pi / (2 * (n + 1))
    T = np.arange(0, 3 * np.pi / 2, h)
    X = np.cos(T)
    Y = np.sin(T)

    fig = plt.figure("Splines cubiques et cercle :-)")
    plt.plot(X, Y, ".r", markersize=10)
    t = np.linspace(0, 2 * np.pi, 100)
    plt.plot(np.cos(t), np.sin(t), "--r")

    t = np.linspace(0, 3 * np.pi / 2, 100)
    x_spline = spline(t, h, X)
    y_spline = spline(t, h, Y)
    plt.plot(x_spline, y_spline, "-b")

    plt.axis("equal")
    plt.axis("off")
    plt.show()


if __name__ == "__main__":
    main()
