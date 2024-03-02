# -------------------------------------------------------------------------
#
# PYTHON for DUMMIES
# Problème 3
#
# Script de test un peu plus rigolo
# Pour introduire un point : faire un clic sur la figure
# Un double clic permet d'obtenir le calcul de la courbe b-spline
# Idem que pour les splines cubiques périodiques !
#
#  Vincent Legat
#
# -------------------------------------------------------------------------
#

import matplotlib
from matplotlib import pyplot as plt
from numpy import *
from devoir3_bspline import bspline


# ====================== callback pour les événements avec la souris ======
#
#  Observer la gestion distincte du clic simple et double :-)
#  Apres un evenement, on redessine la figure avec draw()
#


global X, Y, n


def mouse(event):
    global X, Y, n
    if event.dblclick:
        if n > 2:
            plt.plot([*X, X[0]], [*Y, Y[0]], "--r")
            t = linspace(0, n, n * 1000 + 1)
            x, y = bspline(X, Y, t)
            plt.plot(x, y, "-b")
            X, Y = [], []
            n = 0
    else:
        x = event.xdata
        y = event.ydata
        if x != None and y != None:
            n = n + 1
            X = append(X, [x])
            Y = append(Y, [y])
            print("New data : " + str(x) + "," + str(y))
            plt.plot([x], [y], ".r", markersize=10)
    fig.canvas.draw()


# ============================= mainProgram ===============================

matplotlib.rcParams["toolbar"] = "None"
matplotlib.rcParams["lines.linewidth"] = 1
plt.rcParams["figure.facecolor"] = "lavender"

X, Y = [], []
n = 0
fig = plt.figure("B-spline approximation : p=3")
fig.canvas.mpl_connect("button_press_event", mouse)
plt.ylim((0, 1))
plt.xlim((0, 1.3))
plt.axis("off")

plt.show()
