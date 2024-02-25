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
def get_ddU(U, h):
    n = len(U)
    A = np.zeros((n, n))
    b = np.zeros(n)
    k = (h**2) / 6
    for i in range(n):
        b[i] = U[i - 1] - 2 * U[i] + U[(i + 1) % n]
        A[i][i - 1] = 1
        A[i][i] = 4
        A[i][(i + 1) % n] = 1
    return solve(A, b) / k


def spline(x, h, U):
    print(
        """


    New spline :"""
    )
    U = np.append(U, U[0])
    n = np.size(U)
    X = np.arange(0, n + 1) * h
    i = np.searchsorted(X[1:], x)
    X = X[:-1]

    d = np.full(n, 4)  # vecteur rempli de 4
    A = np.diag(d) + np.diag(np.ones(n - 1), -1) + np.diag(np.ones(n - 1), 1)
    A[-1, 1] = 1
    A[0, -1] = 1
    A[0, 0] = 4
    A[0, 1] = 1
    b = U[range(-1, n - 1)] - 2 * U[0:n] + U[np.append(np.arange(1, n), [0])]
    print("original b", b)
    b[0] = 0
    print(A)
    print(b)
    ddU = solve(A, b) / ((h**2) / 6)
    print("ddU", ddU)
    ddU[0] = ddU[-1]
    n += 1
    A = ddU[: n - 1] / (6 * h)
    B = ddU[1:n] / (6 * h)
    C = U[: n - 1] / h - ddU[: n - 1] * h / 6
    D = U[1:n] / h - ddU[1:n] * h / 6

    print(len(ddU), len(A), len(B), len(C), len(D))
    r = X[i + 1] - x
    s = x - X[i]

    return r * (A[i] * r * r + C[i]) + s * (B[i] * s * s + D[i])


def my_spline(x, h, U):
    #
    # -0- Construction des abscisses (y-compris celle qui correspond au retour au point de depart)
    #     Exemple U = [U_0,U_1,U_2] => n=3 => X = [0,h,2h,3h]
    #
    #

    U = np.append(U, U[0])
    U = np.append(U, U[1])
    n = np.size(U)
    X = np.arange(0, n) * h
    ddU = get_ddU(U, h)

    i = np.zeros(len(x), dtype=int)
    for j in range(1, n):
        i[X[j] <= x] = j
    #
    # A MODIFIER ..... [begin]
    # Ici, on renvoie une interpolation linéaire par morceaux :-)
    # Un peu trop simple non !
    #

    return (
        (ddU[i - 1] / (6 * h)) * (X[i] - x) ** 3
        + (ddU[i] / (6 * h)) * (x - X[i - 1]) ** 3
        + ((U[i - 1] / h) - (ddU[i - 1] * h / 6)) * (X[i] - x)
        + ((U[i] / h) - ((ddU[i] * h) / 6)) * (x - X[i - 1])
    )


#
# A MODIFIER ..... [end]
#

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
