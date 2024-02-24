import matplotlib
from matplotlib import pyplot as plt
from numpy import *
from numpy.polynomial import Polynomial
from scipy.interpolate import CubicSpline as spline
from scipy.interpolate import lagrange as lagrange_interp


def plot_spline():
    x = linspace(0, 25, 35 * 100)
    X = arange(5, 25, 5)
    Y = array([70.2, 70.2, 70.3, 71.2])
    cubic_spline = spline(X, Y)(x)
    lagrange = lagrange_interp(X, Y)
    lagrange = Polynomial(lagrange.coef[::-1])(x)
    print(lagrange[0])
    plt.plot(x, lagrange, "-b", label="lagrange")
    plt.plot(x, cubic_spline, "-r", label="spline cubique")
    plt.legend(loc="upper right")
    plt.show()


plot_spline()
