# -------------------------------------------------------------------------
#
# PYTHON for DUMMIES 23-24
# Problème 4
#
# Script de test
#  Vincent Legat
#  Nathan Coppin
#  Nicolas Roisin
#
# Largement inspiré d'un code préliminaire de Nicolas Roisin :-)
# Ou les méthodes numériques pour obtenir la solution du projet P2 !
#
# -------------------------------------------------------------------------
#

from numpy import *
import numpy as np
from scipy.special import roots_legendre
import matplotlib.pyplot as plt


# ============================================================
# FONCTIONS A MODIFIER [begin]
#
# ------------------------------------------------------------------------------------
#
# Calcul des composantes radiales et verticales
# du champ magnétique généré
# par une série de boucles circulaires de courants !
#
#    X,Z : tableaux numpy des coordonnées x,z des points en [m]
#    Rsource : liste des rayons des boucles en [m]
#    Zsource : liste de hauteurs des boucles en [m]
#    Isource : liste des courants des boucles en [A]
#    data : structure contenant les paramètres matériels
#
#  La fonction renvoie une liste [Bx,Bz] où Bx,Bz sont des tableaux/variables
#  du même type que X,Z.   Le résultat est exprimé en [Tesla] ou [kg /s2 A]
#


def inductanceMegneticField(X, Z, Rsource, Zsource, Isource, data):
    pi = np.pi
    mu0 = data.mu0

    dtheta = 2 * pi / data.nTheta
    x_i = data.Xcircle
    y_i = data.Ycircle

    def biot_savart(x_i, y_i, rboucle, zboucle, iboucle):
        x_i = x_i * rboucle
        y_i = y_i * rboucle
        distance_cube = np.sqrt((X - x_i) ** 2 + y_i**2 + (Z - zboucle) ** 2) ** 3
        dx_i = -dtheta * y_i
        dy_i = dtheta * x_i
        Bx = (Z - zboucle) * dy_i / distance_cube
        Bz = (-y_i * dx_i - (X - x_i) * dy_i) / distance_cube
        k = mu0 * iboucle / (4 * pi)
        return Bx * k, Bz * k

    Bx = np.zeros_like(X)
    Bz = np.zeros_like(Z)

    for i in range(data.nTheta):
        # pour chaque boucle
        for rboucle, zsource, isource in zip(Rsource, Zsource, Isource):
            # calculer le champ magnétique
            B_x, B_z = biot_savart(x_i[i], y_i[i], rboucle, zsource, isource)
            # additionner les contributions
            Bx += B_x
            Bz += B_z

    return [Bx, Bz]


# ------------------------------------------------------------------------------------
#
# Calcul des points et points d'intégration
# pour une intégration double de Gauss-Legendre
# afin de calculer les flux du champs magnétique
#
#    X0,Xf : intervalle radial d'intégration
#    Z0,Zf = intervalle vertical de moyenne du flux
#    nX : nombre de sous-intervalles en x
#    nZ : nombre de sous-intervalles en z
#    nGL : nombre de points de Gauss-Legendre par intervalle (en x et en z)
#
#  La fonction renvoie une liste [X,Z,W] contenant les abscisses et les poids.
#  Ce seront des tableaux unidimensionnels de taille n*nGL*m*nGL
#  Si n ou m sont nuls, on fera une intégration unidimensionnelle en utilisant
#  X0 ou Z0 et les tableaux auront une taille m*nGL ou n*nGL respectivement
#


def inductanceGaussLegendre(X0, Xf, Z0, Zf, n, m, nGaussLegendre):
    #
    # A COMPLETER / MODIFIER
    #

    xi, we = roots_legendre(nGaussLegendre)

    size = max(n * nGaussLegendre, 1) * max(m * nGaussLegendre, 1)
    X = Xf * ones(size)
    Z = Zf * ones(size)
    W = ones(size)

    #
    # A COMPLETER / MODIFIER
    #

    return [X, Z, W]


# ------------------------------------------------------------------------------------
#
# Calcul des points et poids d'intégration
# pour une intégration double de Simpson
# afin de calculer les flux du champs magnétique
#
#    X0,Xf : intervalle radial d'intégration
#    Z0,Zf = intervalle vertical de moyenne du flux
#    nX : nombre de sous-intervalles en x
#    nZ : nombre de sous-intervalles en z
#
#  La fonction renvoie une liste [X,Z,W] contenant les abscisses et les poids.
#  Ce seront des tableaux unidimensionnels de taille (2nX+1)*(2nZ+1)
#  Si nX ou nZ sont nuls, on fera une intégration unidimensionnelle en utilisant
#  X0 ou Z0 et les tableaux auront une taille (2nZ+1), (2nX/1) ou 1 respectivement
#


def inductanceSimpson(X0, Xf, Z0, Zf, n, m):
    #
    # A COMPLETER / MODIFIER
    #

    X = linspace(X0, Xf, 2 * n + 1)
    Z = linspace(Z0, Zf, 2 * m + 1)
    size = (2 * n + 1) * (2 * m + 1)
    w_x = [2, 4] * n + [1]
    w_x[0] = 1
    w_y = [2, 4] * m + [1]
    w_y[0] = 1
    W = zeros(size)
    m_corr = m if m > 0 else 1
    for i in range(2 * n + 1):
        for j in range(2 * m + 1):
            W[i * (2 * m + 1) + j] = (
                (2 * np.pi * X[i]) * w_x[i] * w_y[j] / (9 * 4 * m_corr)
            )

    X, Z = meshgrid(X, Z)
    print(W)
    return [X.flatten(), Z.flatten(), W]


#
# FONCTIONS A MODIFIER [end]
# ============================================================
