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
        dx_i = -dtheta * y_i
        dy_i = dtheta * x_i
        distance_cube = np.sqrt((X - x_i) ** 2 + y_i**2 + (Z - zboucle) ** 2) ** 3
        Bx = (Z - zboucle) * dy_i / distance_cube
        Bz = (-y_i * dx_i - (X - x_i) * dy_i) / distance_cube
        k = mu0 * iboucle / (4 * pi)
        return Bx * k, Bz * k

    Bx = np.zeros_like(X)
    Bz = np.zeros_like(Z)

    # pour chaque boucle
    for rboucle, zsource, isource in zip(Rsource, Zsource, Isource):
        # pour chaque point du cercle discrétisé dans data
        for i in range(data.nTheta):
            # calculer le champ magnétique en tout points (X, Z)
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
    # on récupère les absisses et poids
    xi, we = roots_legendre(nGaussLegendre)

    # on s'assure que n et m sont non nuls
    n, m = max(n, 1), max(m, 1)

    # on initialise les tableaux de points et poids
    size = max(n * nGaussLegendre, 1) * max(m * nGaussLegendre, 1)
    xi = (xi + 1) / 2
    X = zeros(n * nGaussLegendre)
    Z = zeros(m * nGaussLegendre)
    W = ones(size)

    # taille des intervalles d'intégration
    h_x = (Xf - X0) / n
    h_z = (Zf - Z0) / m

    # on remplit les tableaux de points, même méthode pour X et Z :
    # par intervalle, on remplit les nGL points sur cet intervalle
    for x in range(n):
        X[x * nGaussLegendre : (x + 1) * nGaussLegendre] = X0 + (x + xi) * h_x
    for z in range(m):
        Z[z * nGaussLegendre : (z + 1) * nGaussLegendre] = Z0 + (z + xi) * h_z

    # pour le calcul des poids, on utilise une double boucle qui parcourt les points
    for x in range(n * nGaussLegendre):
        for j in range(m * nGaussLegendre):
            W[j * n * nGaussLegendre + x] = (
                2
                * pi
                * h_x
                * X[x]
                * we[x % nGaussLegendre]
                * we[j % nGaussLegendre]
                / (4 * m)
            )

    # même trick que pour Simpson pour "entremêler" les points 1D -> 2D
    X, Z = meshgrid(X, Z)
    X = X.flatten()
    Z = Z.flatten()
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
    # points d'intégration de Simpson
    X = linspace(X0, Xf, 2 * n + 1)
    Z = linspace(Z0, Zf, 2 * m + 1)

    # taille d'un sous-intervalle
    h_x = (Xf - X0) / (2 * max(n, 1))

    # Construction des poids unidimensionnels pour Simpson composite sur l'intervalle
    w_x = array([2, 4] * n + [1])
    w_x[0] = 1
    w_y = array([2, 4] * m + [1])
    w_y[0] = 1
    w_x = w_x / 3
    w_y = w_y / 3

    # initialisation du tableau des poids
    size = (2 * n + 1) * (2 * m + 1)
    W = zeros(size)

    # on remplit le tableau des points (1D) comme un tableau 2D en ajustant les indices
    for i in range(2 * n + 1):
        for j in range(2 * m + 1):
            W[j * (2 * n + 1) + i] = (2 * np.pi * X[i] * h_x * w_x[i] * w_y[j]) / (
                4 * max(m, 1)
            )
            # cette formule est donnée dans l'énoncé

    # cette fonction retourne des grilles de points 2D, qu'il faut "aplatir" ensuite
    X, Z = np.meshgrid(X, Z)
    return [X.flatten(), Z.flatten(), W]


#
# FONCTIONS A MODIFIER [end]
# ============================================================
