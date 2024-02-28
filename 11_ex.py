import matplotlib
from matplotlib import pyplot as plt
from numpy import *
from numpy.polynomial import Polynomial
from scipy.interpolate import CubicSpline as spline
from scipy.interpolate import lagrange as lagrange_interp


def sin_interp():
    x = linspace(-1, 1, 200)
    X = linspace(-1, 1, 21)
    noise = array([(-1) ** (k + 1) * (10 ** (-4)) for k in range(21)])
    Y = sin(2 * pi * X)
    Y_pert = Y + noise
    cubic_spline = spline(X, Y)(x)
    cubic_spline_pert = spline(X, Y_pert)(x)
    lagrange = Polynomial(lagrange_interp(X, Y).coef[::-1])(x)
    lagrange_pert = Polynomial(lagrange_interp(X, Y_pert).coef[::-1])(x)
    plt.plot(x, sin(2 * pi * x), "-r", label="sinus")
    plt.plot(x, lagrange, "-b", label="lagrange")
    plt.plot(x, lagrange_pert, "-g", label="lagrange")
    # plt.plot(x, cubic_spline, "-r", label="spline cubique")
    # plt.plot(x, cubic_spline_pert, "-g", label="spline cubique")
    plt.legend(loc="upper right")
    plt.show()


sin_interp()
