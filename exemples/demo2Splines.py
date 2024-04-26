#
# Splines cubiques
# Vincent Legat - 2018
# Ecole Polytechnique de Louvain
#

import matplotlib
from matplotlib import pyplot as plt
from numpy import *

n = 21
L = 1
X = linspace(-L, L, n)  # ; print(X)  # Pour debugger step by step
U = sin(2 * pi * X)  # ; print(U)  # Pour debugger step by step

x = linspace(-L, L, 10 * n)

#
# -1- Accéder aux données d'un tableau
#     Assez proche de la syntaxe de MATLAB
#     Attention : on numérote à partir de zero en Python
#       X[0:21:2] Python = MATLAB X(1:2:20)
#

print(X)
print(X[:])
print(X[0:10])
print(X[0:21:2])

#
# -2- Une copie et une vue d'un tableau
#     Attention : ceci est différent de MATLAB !
#     =====> beaucoup de bugs vicieux possibles....
#     Y = une vue du même espace mémoire
#     Z = une vraie copie du mémoire
#

Y = X[0:21:2]
Y[0] = 69
print("X[0] = %3d - Y[0] = %3d" % (X[0], Y[0]))
X[0] = 456
print("X[0] = %3d - Y[0] = %3d" % (X[0], Y[0]))
X[0] = -1
print("X[0] = %3d - Y[0] = %3d" % (X[0], Y[0]))
Z = copy(X)
Z[0] = 69
print("X[0] = %3d - Y[0] = %3d - Z[0] = %3d" % (X[0], Y[0], Z[0]))
X[0] = 456
print("X[0] = %3d - Y[0] = %3d - Z[0] = %3d" % (X[0], Y[0], Z[0]))
X[0] = -1
print("X[0] = %3d - Y[0] = %3d - Z[0] = %3d" % (X[0], Y[0], Z[0]))

#
# -3- Calcul des splines cubiques
#

from scipy.interpolate import CubicSpline as spline

uSpline1 = spline(X[0:10], U[0:10])(x)
uSpline2 = spline(X[0:21:2], U[0:21:2])(x)

#
# -4- Et le joli dessin :-)
#


def plot_spline():
    plt.plot(x, uSpline1, "-b", label="spline sur les 10 premiers point")
    plt.plot(x, uSpline2, "-r", label="spline sur un point sur deux")
    plt.plot(X[0:10], U[0:10], ".r", markersize=20, label="10 premiers points")
    plt.plot(X[0:21:2], U[0:21:2], ".b", markersize=10, label="1 point sur 2")
    plt.legend(loc="upper right")
    plt.show()


def quart_de_cercle():
    k = arange(-3, 4)
    X = sin(k * pi / 6)
    Y = cos(k * pi / 6)
    T = k * pi / 6

    x = linspace(-3, 4, 40 * n)
    sinus = sin(x)
    cosinus = cos(x)
    u_spline = spline(X, Y)(x)
    sin_spline = spline(T, X)(x)
    cos_spline = spline(T, Y)(x)

    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.plot(sin_spline, cos_spline, label="sin")
    plt.plot(x, u_spline, label="sin")
    plt.legend()
    plt.show()


quart_de_cercle()
