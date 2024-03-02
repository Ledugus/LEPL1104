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

from numpy import *

# ============================================================
# FONCTIONS A MODIFIER [begin]
#


def b(t, T, i, p):
    if p == 0:
        return (T[i] <= t) * (t < T[i + 1])
    else:
        u = (
            0.0
            if T[i + p] == T[i]
            else (t - T[i]) / (T[i + p] - T[i]) * b(t, T, i, p - 1)
        )
        u += (
            0.0
            if T[i + p + 1] == T[i + 1]
            else (T[i + p + 1] - t) / (T[i + p + 1] - T[i + 1]) * b(t, T, i + 1, p - 1)
        )
        return u


def bspline(X, Y, t):
    x = t.copy()
    y = t.copy()
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
    t = linspace(0, len(X), len(X) * 100 + 1)

    x, y = bspline(X, Y, t)

    #
    # -2- Un joli dessin :-)
    #

    import matplotlib.pyplot as plt
    import matplotlib

    matplotlib.rcParams["toolbar"] = "None"
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
