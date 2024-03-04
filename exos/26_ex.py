from numpy import *
from numpy.linalg import solve
from scipy.interpolate import lagrange as lagrange_interp
from numpy.polynomial import Polynomial
import matplotlib.pyplot as plt
import matplotlib


def b(t, T, i, p):
    """Retourne une b-spline de base
    selon les noeuds T,
    d'indice i,
    d'ordre p
    """
    if p == 0:
        return (T[i] <= t) * (t < T[i + 1])
    u = 0.0 if T[i + p] == T[i] else (t - T[i]) / (T[i + p] - T[i]) * b(t, T, i, p - 1)
    u += (
        0.0
        if T[i + p + 1] == T[i + 1]
        else (T[i + p + 1] - t) / (T[i + p + 1] - T[i + 1]) * b(t, T, i + 1, p - 1)
    )
    return u


def bspline(X, Y, t, p):
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


def spline(t, h, U):
    """Retourne la spline cubique périodique,
    selon des points d'interpolation réguliers espacés de h
    """
    n = size(U)

    T = arange(0, n + 1) * h
    i = searchsorted(T[1:], t)

    # Création de la matrice A du système des ddu
    d = full(n, 4)
    A = diag(d) + diag(ones(n - 1), -1) + diag(ones(n - 1), 1)
    A[-1, 0] = 1
    A[0, -1] = 1

    # Création du vecteur b
    b = U[range(-1, n - 1)] - 2 * U[0:n] + U[append(arange(1, n), [0])]
    dd_u = solve(A, b) / ((h**2) / 6)

    # On ajoute les valeurs extrêmes de l'intervalle pour la périodicité
    U = append(U, U[0])
    dd_u = append(dd_u, dd_u[0])

    # a, b, c et d sont les coefficients du polynome
    n += 1
    a = dd_u[: n - 1] / (6 * h)
    b = dd_u[1:n] / (6 * h)
    c = U[: n - 1] / h - dd_u[: n - 1] * h / 6
    d = U[1:n] / h - dd_u[1:n] * h / 6

    r = T[i + 1] - t
    s = t - T[i]

    return r * (a[i] * r * r + c[i]) + s * (b[i] * s * s + d[i])


def lagrange(t, T, U):
    """Retourne l'interpolation polynomiale
    Absisses : T
    Ordonnées : U
    intervalle d'interpolation  : t
    """
    return Polynomial(lagrange_interp(T, U).coef[::-1])(t)


def main():
    # initialisation des variables et des points d'observation
    T = arange(0, 4, 1)
    X = array([0, 1, 0])
    Y = array([0, 0, 1])
    t = linspace(0, len(X), len(X) * 100 + 1)

    # Génération des 3 courbes
    x_b, y_b = bspline(X, Y, t, 2)
    x, y = spline(t, 1, X), spline(t, 1, Y)
    x_pol, y_pol = lagrange(t, T, [*X, X[0]]), lagrange(t, T, [*Y, Y[0]])

    # Paramètres de la fenềtre
    matplotlib.rcParams["toolbar"] = "None"
    plt.rcParams["figure.facecolor"] = "white"
    fig = plt.figure("Approximation avec des B-splines")

    # plot des points d'obervation
    plt.plot(X, Y, ".k", markersize=20)

    # plot des points de contact de la courbe b-spline
    plt.plot(x_b[T * 100], y_b[T * 100], ".b")

    # plot des 3 courbes
    plt.plot(x_b, y_b, "-b")
    plt.plot(x_pol, y_pol, "-r")
    plt.plot(x, y, "-g")

    print("Valeurs en t = 1/2")
    print("b", x_b[50], y_b[50])
    print("pol", x_pol[50], y_pol[50])
    print("cubic", x[50], y[50])

    plt.axis("equal")
    plt.axis("off")
    plt.show()


if __name__ == "__main__":
    main()
