# -------------------------------------------------------------------------
#
# PYTHON for DUMMIES
# Problème de Lorenz
#
# Script de test
#  Vincent Legat
#
# -------------------------------------------------------------------------
#

from numpy import *

# ============================================================
# FONCTION A MODIFIER [begin]
#
#


def lorenz(Xstart, Xend, Ustart, n):

    X = linspace(Xstart, Xend, n + 1)
    U = zeros((n + 1, 3))
    U[0] = Ustart

    # implement runge kutta for lorentz equations
    def f(U):
        x, y, z = U
        return array([10 * (y - x), x * (28 - z) - y, x * y - 8 / 3 * z])

    for i in range(n):
        h = (Xend - Xstart) / n
        k1 = h * f(U[i])
        k2 = h * f(U[i] + k1 / 2)
        k3 = h * f(U[i] + k2 / 2)
        k4 = h * f(U[i] + k3)
        U[i + 1] = U[i] + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return X, U


#
# FONCTION A MODIFIER [end]
# ============================================================


def main():

    # ------------------------------------------------------------------------------------
    #
    # Script de test
    #
    #
    # ------------------------------------------------------------------------------------

    from matplotlib import pyplot as plt

    plt.rcParams["figure.facecolor"] = "lavender"
    plt.rcParams["axes.facecolor"] = "lavender"

    plt.figure("Lorenz Equations")
    Xstart = 0
    Xend = 100
    Ustart = [0, 1, 0]
    n = 10000

    X, U = lorenz(Xstart, Xend, Ustart, n)
    plt.plot(U[:, 0], U[:, 2], "-r", linewidth=0.5)
    plt.show()


main()
