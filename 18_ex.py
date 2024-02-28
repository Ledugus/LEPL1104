from numpy import *
from scipy.sparse import csr_matrix
from scipy.sparse import spdiags
from scipy.sparse.linalg import spsolve
from scipy.integrate import quad
import matplotlib.pyplot as plt


def u(x):
    return 1 / ((x - 0.3) ** 2 + 0.01) + 1 / ((x - 0.9) ** 2 + 0.04) - 6


def my_alphonse(u, X):
    n = len(X) - 1
    A = zeros((n, n))
    b = zeros(n)

    return solve(A, b)


def alphonse(X):
    """Calcul de la meilleure approximation au sens de la norme L2
    de la fonction humps au moyen d’une fonction lin´eaire par morceaux
    sur l’intervalle ]X(0),X(-1)["""
    n = len(X)  # nombre de points
    e = n - 1  # nombre d’´el´ements
    h = diff(X)  # taille des ´el´ements (c’est un vecteur !)
    #
    # -1- Assemblage de la matrice et du membre de droite
    # La matrice est calcul´ee analytiquement et est tridiagonale
    # Le membre de droite est int´egr´e num´eriquement via quad
    # -- Important : utilisation des matrices creuses !
    #
    hlow = array([*h, 0]) / 6
    hup = array([0, *h]) / 6
    print(h)
    print(hlow)
    print(hup)
    A = spdiags([hlow, 2 * (hup + hlow), hup], [-1, 0, 1], n, n)
    B = zeros((n, 1))
    for elem in range(e):
        Xleft, Xright = X[elem], X[elem + 1]
        a = quad(lambda x: u(x) * (Xright - x) / (Xright - Xleft), Xleft, Xright)[0]
        B[elem] += a
        b = quad(lambda x: u(x) * (x - Xleft) / (Xright - Xleft), Xleft, Xright)[0]
        B[elem + 1] += b
    #
    # -2- R´esolution du syst`eme
    #
    print(B)
    print(csr_matrix(A))
    return spsolve(csr_matrix(A), B)


def plot():
    x = linspace(0, 3, 100)
    X = linspace(0, 3, 5)
    plt.plot(x, u(x), "-b")
    plt.plot(X, alphonse(X), "-g")
    plt.show()


if __name__ == "__main__":
    plot()
