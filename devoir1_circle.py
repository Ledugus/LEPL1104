#
# PYTHON for DUMMIES 23-24
# Problème 1
#
# Un cercle roule sans glisser autour d'un autre cercle.
#
# Vincent Legat
# (avec l'aide bienveillante de copilotes anonymes)
#
# -------------------------------------------------------------------------
#
import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

# ============================================================
# FONCTIONS A MODIFIER [begin]
#
# -1- Création d'un cercle...
#     radius : rayon du cercle
#     n : nombre de points répartis entre [0,2*pi] pour tracer le cercle
#     il y aura donc n-1 arcs de cercles :-)
#
#     theta : tableau numpy contenant les angles pour chaque point
#     x,y : tableaux numpy contenant les coordonnées x,y pour chaque point
#


def circlesCreate(radius, n):
    theta = np.linspace(0, 2 * np.pi, n)  # à modifier !
    x = np.cos(theta) * radius
    y = np.sin(theta) * radius

    return theta, x, y


#
# -2- Calcul des angles où le petit cercle aura retrouvé sa position initiale
#     ratio : rapport entre le rayon du cercle intérieur et le rayon du cercle extérieur
#
#     thetas : tableau numpy contenant les angles où le cercle intérieur aura retrouvé sa position initiale


def circlesAngles(ratio):
    thetas = np.arange(0, ratio + 1, 1)
    thetas = thetas * ((2 * np.pi) / (ratio + 1))
    return thetas


#
# -3- Animation des deux cercles...
#     theta : angle du point de contact du petit cercle avec le grand cercle
#     x_rolling,y_rolling : tableaux numpy contenant les coordonnées x,y pour les points du cercle roulant
#     ratio : rapport entre le rayon du cercle intérieur et le rayon du cercle extérieur
#
#     x,y : tableaux numpy contenant les coordonnées x,y pour les points du cercle roulant pour l'angle theta
#


def circlesAnimate(theta, x_rolling, y_rolling, ratio):
    x = (
        (ratio + 1) * np.cos(theta)
        + x_rolling * np.cos(theta * (ratio + 1))
        - y_rolling * np.sin(theta * (ratio + 1))
    )
    y = (
        (ratio + 1) * np.sin(theta)
        + y_rolling * np.cos(theta * (ratio + 1))
        + x_rolling * np.sin(theta * (ratio + 1))
    )

    return x, y


#
# FONCTIONS A MODIFIER [end]
# ============================================================
#
# -1- Test des fonctions et animation :-)


def main():
    #
    # -2- Paramètres du problème
    #     ratio : rapport entre le rayon du cercle intérieur et le rayon du cercle extérieur
    #     n_steps : nombre de pas pour une rotation complète du petit cercle
    #     n : nombre de pas pour que le petit cercle revienne a sa position initiale
    #     i_mid : indice du point initial droit du cercle roulant
    #             (le point initial gauche est à l'indice 0)
    #

    ratio = 4
    n_steps = 63
    n = (ratio + 1) * n_steps + 1
    i_mid = n_steps // 2

    #
    # -3- Création des cercles
    #

    theta_inner, x_inner, y_inner = circlesCreate(ratio, n)
    theta_rolling, x_rolling, y_rolling = circlesCreate(1, n_steps)

    #
    # -4- Et zou, un joli plot et une jolie animation :-)
    #     x_left,y_left : tableaux numpy contenant les coordonnées x,y pour le point gauche du cercle roulant
    #     x_right,y_right : tableaux numpy contenant les coordonnées x,y pour le point droit du cercle roulant
    #

    plt.rcParams["toolbar"] = "None"
    fig = plt.figure("Circle rolling around another circle")

    x_left = np.ones_like(x_inner) * np.nan
    y_left = np.ones_like(y_inner) * np.nan
    x_right = np.ones_like(x_inner) * np.nan
    y_right = np.ones_like(y_inner) * np.nan

    #
    # -4.1- Plot des trajectoires et des animations
    #       line1,line2,Line3,Line4,Line5 : objets pour les animations qui permettront de mettre à jour les plots
    #

    plt.fill(x_inner, y_inner, "b-")
    (line1,) = plt.plot([], [], "r-")
    (line2,) = plt.plot([], [], "b-")
    (line3,) = plt.plot([], [], "k-")
    (line4,) = plt.plot([], [], "ko-")
    (line5,) = plt.plot([], [], "bo")
    #
    # -4.2- Plot du cercle roulant dans son configuration initiale
    #       et du grand cercle fixe en bleu
    #       (ces plots ne seront pas mis à jour par l'animation)

    for theta in circlesAngles(ratio):
        x, y = circlesAnimate(theta, x_rolling, y_rolling, ratio)
        plt.fill(x[0 : i_mid + 1], y[0 : i_mid + 1], color="b", alpha=0.5)
        plt.fill(x[i_mid:], y[i_mid:], color="y", alpha=0.5)

    #
    # -4.3- Définition de la fonction callback pour l'animation
    #       (elle sera appelée à chaque pas de temps)
    #       (elle mettra à jour les plots)
    #

    def animate(i):
        [x, y] = circlesAnimate(theta_inner[i], x_rolling, y_rolling, ratio)

        x_left[i] = x[0]
        y_left[i] = y[0]
        x_right[i] = x[i_mid]
        y_right[i] = y[i_mid]

        line1.set_data(x, y)
        line2.set_data(x_left, y_left)
        line3.set_data(x_right, y_right)
        line4.set_data([x[0], x[i_mid]], [y[0], y[i_mid]])
        line5.set_data([x[0]], [y[0]])

        return line1, line2, line3, line4, line5

    #
    # -4.4- Lancement de l'animation
    #       (elle appelle la fonction callback à chaque pas de temps)
    #

    animation = FuncAnimation(fig, animate, frames=n, interval=50)

    size = 2.2 + ratio
    ax = plt.gca()
    ax.set_xlim(-size, size)
    ax.set_ylim(-size, size)
    ax.set_aspect("equal")
    plt.axis("off")

    #
    # -4.5- Sauvegarde de l'animation dans un fichier
    #       (si ffmpeg est installé sur votre ordinateur)
    #
    # animation.save('circle_animation.mp4', writer='ffmpeg')
    #

    plt.show()


if __name__ == "__main__":
    main()
