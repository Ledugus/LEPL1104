# -------------------------------------------------------------------------
#
# PYTHON for DUMMIES 23-24
# Problème 7
#
# Calcul de la zone de stabilité d'une méthode de Gear
#  Vincent Legat
#
# -------------------------------------------------------------------------
#

import numpy as np
from numpy import *

# ============================================================
# FONCTIONS A MODIFIER [begin]
#
# -1- La zone de stabilité de Gear d'ordre un (qui est la méthode
#     d'Euler implicite :-) vous est offerte gracieusement
#     par les descendants de la famille Gear !
#
#     Evidemment, c'est plus compliqué pour les autres ordres !
#


def stabilityGear(x, y, order):

    #
    # A COMPLETER / MODIFIER [begin]
    #

    z = x + 1j * y
    coeffs = [
        np.array([1, 1]),
        np.array([2, 4, -1]) / 3,
        np.array([6, 18, -9, 2]) / 11,
        np.array([12, 48, -36, 16, -3]) / 25,
        np.array([60, 300, -300, 200, -75, 12]) / 137,
        np.array([60, 360, -450, 400, -225, 72, -10]) / 147,
    ]

    if order > 6:
        return zeros(order + 1), zeros_like(x)
    coeff = coeffs[order - 1]

    @vectorize
    def solve(z):
        coeff_to_solve = [coeff[0] * z - 1, *coeff[1:]]
        roots = np.roots(coeff_to_solve)
        return max(abs(roots))

    gain = solve(z)
    #
    # A COMPLETER / MODIFIER [end]
    #

    return gain, coeff


#
# FONCTIONS A MODIFIER [end]
# ============================================================
#
# ------------------------------------------------------------------------------------
#
# Script de test
#
#
# -1- Choix de l'ordre et de la précision de la carte de stabilité
#
# ------------------------------------------------------------------------------------


def main():
    order = 3
    n = 100
    x, y = np.meshgrid(np.linspace(-8, 8, n), np.linspace(-8, 8, n))

    gain, coeff = stabilityGear(x, y, order)

    #
    # -2- Impression des coefficients de la méthode de Gear
    #     Observer l'utilisation du module "fractions" qui permet
    #     d'écrire des nombres en virgule flottante sous la forme de fractions !
    #
    #     Retirer la ligne np.set_printoptions
    #     Noter aussi les astuces dans les paramètres :-)
    #     Merci à : https://stackoverflow.com/questions/42209365/numpy-convert-decimals-to-fractions
    #
    #

    import fractions

    np.set_printoptions(
        formatter={"all": lambda x: str(fractions.Fraction(x).limit_denominator())}
    )
    print("==== Coefficients of Gear's method for order = %d ========== " % order)
    print("     ", end="")
    print(coeff)

    #
    # -3- Faire le joli plot
    #

    import matplotlib.pyplot as plt
    import matplotlib

    matplotlib.rcParams["toolbar"] = "None"
    plt.rcParams["figure.facecolor"] = "lavender"
    plt.rcParams["axes.facecolor"] = "lavender"
    plt.figure("Stability of Gear's method : order = %d" % order)
    plt.contourf(x, y, gain, np.arange(0, 1.1, 0.1), cmap=plt.cm.jet_r)
    plt.contour(x, y, gain, np.arange(0, 1.1, 0.1), colors="black", linewidths=0.5)
    ax = plt.gca()
    ax.axhline(y=0, color="r")
    ax.axvline(x=0, color="r")
    ax.yaxis.grid(color="gray", linestyle="dashed")
    ax.xaxis.grid(color="gray", linestyle="dashed")
    ax.set_aspect("equal", "box")
    plt.show()


main()
