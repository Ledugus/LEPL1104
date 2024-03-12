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

    #
    # A COMPLETER / MODIFIER
    #
    Bx = ones(shape(X))
    Bz = ones(shape(X))

    #
    # A COMPLETER / MODIFIER
    #

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
# Calcul des points et points d'intégration
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

    size = (n + 1) * (m + 1)
    X = Xf * ones(size)
    Z = Zf * ones(size)
    W = ones(size)

    #
    # A COMPLETER / MODIFIER
    #

    return [X, Z, W]


#
# FONCTIONS A MODIFIER [end]
# ============================================================


# ------------------------------------------------------------------------------------
#
# Classe du projet P2
#


class MutualInductanceProject:
    mu0 = 4e-7 * pi  # permeabilité du vide en [H/m]

    Rcoil = 1.0e-2  # rayon des bobines [m]
    Hcoil = 1.0e-2  # épaisseur des bobines [m]
    Zcoil = -0.5e-2  # position verticale de la bobine primaire [m]
    nSpires = 200  # nombre de spires
    I = 0.1  # courant dans la bobine [A]
    Hsecond = 2.0e-2  # position relative de la bobine secondaire [m]

    nZcoil = 20  # discrétisation vertical de la bobine
    nTheta = 8  # discrétisation du cercle
    Theta = linspace(0, 2 * pi, nTheta + 1)[:-1] + pi / nTheta
    Xcircle = cos(Theta)
    Ycircle = sin(Theta)

    def streamPlot(self, X, Z, Bx, Bz, colorStream, colorIplate):
        plt.streamplot(X, Z, Bx, Bz, density=1.0, linewidth=None, color=colorStream)
        x = array([-self.Rcoil, self.Rcoil])
        z = repeat([linspace(0, -self.Hcoil, self.nZcoil)], 2, axis=0)
        plt.plot(x, z, "o-r", linewidth=2)
        plt.plot(x, z + self.Hsecond, "o-b", linewidth=2)
        plt.xlim((-0.03, 0.03))
        plt.ylim((-0.03, 0.03))


def main():

    # ------------------------------------------------------------------------------------
    #
    # Script de test
    #
    #
    # -0- Initialisation du projet
    #
    # ------------------------------------------------------------------------------------

    p = MutualInductanceProject()

    # ------------------------------------------------------------------------------------
    #
    # -1- Calcul de l'inductance, inductance mutuelle, inductance induite
    #     et inductance effective pour Mr. Oestges !
    #
    # ------------------------------------------------------------------------------------

    Zcoil = -linspace(0, p.Hcoil, p.nZcoil)
    Rcoil = p.Rcoil * ones_like(Zcoil)
    Icoil = p.I * ones_like(Zcoil) * p.nSpires / p.nZcoil

    #
    # -1.1- Flux dans la bobine primaire du champs induit par la bobine primaire
    #

    nX = 10
    nZ = 3
    nGL = 2
    X, Z, W = inductanceGaussLegendre(0, p.Rcoil, 0, -p.Hcoil / 2, nX, nZ, nGL)
    Bz = inductanceMegneticField(X, Z, Rcoil, Zcoil, Icoil, p)[1]
    flux = W @ Bz
    L = p.nSpires * (flux / p.I)
    print(f"Using Gauss-Legendre integration rule (nX = {nX} nZ = {nZ}):-)")
    print("Flux across primary        = %14.7e [Tesla m2]" % flux)
    print("Inductance                 = %14.7e [Henry]" % L)

    nX = 10
    nZ = 3
    X, Z, W = inductanceSimpson(0, p.Rcoil, 0, -p.Hcoil / 2, nX, nZ)
    Bz = inductanceMegneticField(X, Z, Rcoil, Zcoil, Icoil, p)[1]
    flux = W @ Bz
    L = p.nSpires * (flux / p.I)
    print("Using Simpson integration rule (nX = %d nZ = %d):-)" % (nX, nZ))
    print("Flux across primary        = %14.7e [Tesla m2]" % flux)
    print("Inductance                 = %14.7e [Henry]" % L)
    nX = 20
    nZ = 5
    X, Z, W = inductanceSimpson(0, p.Rcoil, 0, -p.Hcoil / 2, nX, nZ)
    Bz = inductanceMegneticField(X, Z, Rcoil, Zcoil, Icoil, p)[1]
    flux = W @ Bz
    L = p.nSpires * (flux / p.I)
    print("Using Simpson integration rule (nX = %d nZ = %d):-)" % (nX, nZ))
    print("Flux across primary        = %14.7e [Tesla m2]" % flux)
    print("Inductance                 = %14.7e [Henry]" % L)

    #
    # -1.2- Flux dans la bobine secondaire du champs induit par la bobine primaire
    #

    nX = 20
    nZ = 0
    nGL = 1
    X, Z, W = inductanceGaussLegendre(
        0, p.Rcoil, p.Hsecond, -p.Hcoil + p.Hsecond, nX, nZ, nGL
    )
    Bz = inductanceMegneticField(X, Z, Rcoil, Zcoil, Icoil, p)[1]
    flux = W @ Bz
    M = p.nSpires * (flux / p.I)
    print("Flux across secondary      = %14.7e [Tesla m2]" % flux)
    print("Mutual inductance          = %14.7e [Henry]" % M)

    # ------------------------------------------------------------------------------------
    #
    # -2- Représentation du champ magnétique
    #
    # ------------------------------------------------------------------------------------

    plt.rcParams["toolbar"] = "None"

    n = 20
    X, Z = meshgrid(linspace(-0.03, 0.03, n), linspace(-0.03, 0.03, n))
    Bcoilx, Bcoilz = inductanceMegneticField(X, Z, Rcoil, Zcoil, Icoil, p)

    plt.figure("Mutal induction :-)", figsize=(10, 10))
    p.streamPlot(X, Z, Bcoilx, Bcoilz, "blue", "gray")
    plt.title("Electromagnetic field generated by the primary coil")

    # ------------------------------------------------------------------------------------
    #
    # -3- Bz dans la bobine et la plaque
    #
    # ------------------------------------------------------------------------------------

    plt.figure("Champ magnétique : Bz dans les bobines primaire et secondaire")

    n = 100

    X = linspace(0, p.Rcoil, n)
    Z = p.Zcoil * ones_like(X)
    Bz = inductanceMegneticField(X, Z, Rcoil, Zcoil, Icoil, p)[1]
    plt.plot(X, Bz, "-b")
    Z = (p.Zcoil - p.Hcoil / 2) * ones_like(X)
    Bz = inductanceMegneticField(X, Z, Rcoil, Zcoil, Icoil, p)[1]
    plt.plot(X, Bz, "-b")
    Z = (p.Zcoil + p.Hcoil / 2) * ones_like(X)
    Bz = inductanceMegneticField(X, Z, Rcoil, Zcoil, Icoil, p)[1]
    plt.plot(X, Bz, "-b")

    X = linspace(0, p.Rcoil, n)
    Z = p.Zcoil * ones_like(X) + p.Hsecond
    Bz = inductanceMegneticField(X, Z, Rcoil, Zcoil, Icoil, p)[1]
    plt.plot(X, Bz, "-r")
    Z = (p.Zcoil - p.Hcoil / 2) * ones_like(X) + p.Hsecond
    Bz = inductanceMegneticField(X, Z, Rcoil, Zcoil, Icoil, p)[1]
    plt.plot(X, Bz, "-r")
    Z = (p.Zcoil + p.Hcoil / 2) * ones_like(X) + p.Hsecond
    Bz = inductanceMegneticField(X, Z, Rcoil, Zcoil, Icoil, p)[1]
    plt.plot(X, Bz, "-r")

    plt.grid()

    # ------------------------------------------------------------------------------------
    #
    # -4- Points d'intégration de Gauss-Legendre
    #
    # ------------------------------------------------------------------------------------

    #   plt.rcParams['toolbar'] = 'None'
    #   plt.rcParams['toolbar'] = 'None'
    #   plt.rcParams['figure.facecolor'] = 'lavender'
    #   plt.rcParams['axes.facecolor'] = 'lavender'
    plt.figure("Integration Gauss-Legendre nodes")

    nX = 5
    nZ = 4
    nGL = 3
    X, Z, W = inductanceGaussLegendre(0, p.Rcoil, 0, -p.Hcoil, nX, nZ, nGL)
    plt.plot(X, Z, "ob", markersize=5)

    x = array([0, p.Rcoil])
    z = repeat([linspace(0, -p.Hcoil, nZ + 1)], 2, axis=0)
    plt.plot(x, z, "-k")
    x = array([-0.003, p.Rcoil])
    z = repeat([linspace(0, -p.Hcoil, nZ + 1)], 2, axis=0)
    plt.plot(x, z, "--k")
    z = array([0, -p.Hcoil])
    x = repeat([linspace(0, p.Rcoil, nX + 1)], 2, axis=0)
    plt.plot(x, z, "-k")
    z = array([0.003, -p.Hcoil])
    x = repeat([linspace(0, p.Rcoil, nX + 1)], 2, axis=0)
    plt.plot(x, z, "--k")

    X, Z, W = inductanceGaussLegendre(0, p.Rcoil, 0.002, 0.002, nX, 0, nGL)
    plt.plot(X, Z, "or", markersize=5)
    X, Z, W = inductanceGaussLegendre(-0.002, -0.002, 0, -p.Hcoil, 0, nZ, nGL)
    plt.plot(X, Z, "or", markersize=5)
    plt.axis("equal")
    plt.axis("off")

    # ------------------------------------------------------------------------------------
    #
    # -5- Points d'intégration de Simpson
    #     Oui : il faudrait ecrire une fonction et ne pas tout recopier !
    #     Mais, c'est dimanche.... et j'ai autre chose à faire :-)
    #
    # ------------------------------------------------------------------------------------

    plt.figure("Integration Simpson nodes")
    X, Z, W = inductanceSimpson(0, p.Rcoil, 0, -p.Hcoil, nX, nZ)
    plt.plot(X, Z, "ob", markersize=5)

    x = array([0, p.Rcoil])
    z = repeat([linspace(0, -p.Hcoil, nZ + 1)], 2, axis=0)
    plt.plot(x, z, "-k")
    x = array([-0.003, p.Rcoil])
    z = repeat([linspace(0, -p.Hcoil, nZ + 1)], 2, axis=0)
    plt.plot(x, z, "--k")
    z = array([0, -p.Hcoil])
    x = repeat([linspace(0, p.Rcoil, nX + 1)], 2, axis=0)
    plt.plot(x, z, "-k")
    z = array([0.003, -p.Hcoil])
    x = repeat([linspace(0, p.Rcoil, nX + 1)], 2, axis=0)
    plt.plot(x, z, "--k")

    X, Z, W = inductanceSimpson(0, p.Rcoil, 0.002, 0.002, nX, 0)
    plt.plot(X, Z, "or", markersize=5)
    X, Z, W = inductanceSimpson(-0.002, -0.002, 0, -p.Hcoil, 0, nZ)
    plt.plot(X, Z, "or", markersize=5)
    plt.axis("equal")
    plt.axis("off")

    # ------------------------------------------------------------------------------------
    #
    # -6- La jolie courbe pour Claude Oestges et Nicolas Roisin
    #
    # ------------------------------------------------------------------------------------

    nD = 21
    nX = 5
    nZ = 5
    nGL = 3
    distances = linspace(0, 5e-2, nD)
    M = zeros(nD)
    for i in range(nD):
        h = distances[i] + 0.01
        X, Z, W = inductanceGaussLegendre(0, p.Rcoil, h, h - p.Hcoil, nX, nZ, nGL)
        Bz = inductanceMegneticField(X, Z, Rcoil, Zcoil, Icoil, p)[1]
        flux = W @ Bz
        M[i] = p.nSpires * (flux / p.I)

    plt.figure("Inductuance mutuelle en fonction de la distance")
    plt.plot(distances * 1e2, M * 1e3, "ob-")
    ax = plt.gca()
    ax.set_ylabel("Inductance mutuelle : M [mH]")
    ax.set_xlabel("Distance bobines : d [cm]")
    ax.grid()

    plt.show()


main()
