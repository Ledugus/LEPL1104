from matplotlib import cm, colormaps
from matplotlib import pyplot as plt
from numpy import *
from numpy.polynomial import Polynomial
from scipy.interpolate import CubicSpline as spline


T = array([1, 1, 1, 1])


@vectorize
def interpolate_2D_linear(x, y):
    X = array([-1, 1])
    bottom = spline(X, [T[0], T[3]])(x)
    top = spline(X, [T[1], T[2]])(x)
    return spline(X, [bottom, top])(y)


@vectorize
def interpolate_correction(x, y):
    phi = [
        (1 - x) * (1 - y) / 4,
        (1 - x) * (1 + y) / 4,
        (1 + x) * (1 + y) / 4,
        (1 + x) * (1 - y) / 4,
    ]
    return matmul(T, phi)


def show_3d():
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    # Make data
    x = linspace(-1, 1, 50)
    y = linspace(-1, 1, 50)
    X, Y = meshgrid(x, y)
    Z = interpolate_correction(X, Y)

    # Plot the surface
    ax.plot_surface(X, Y, Z, cmap=cm.plasma)

    # Set an equal aspect ratio
    ax.set_aspect("equal")

    plt.show()


show_3d()
