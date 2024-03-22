from numpy import *
from numpy.polynomial import Polynomial
from scipy.interpolate import lagrange as lagrange_interp
from matplotlib import pyplot as plt


def f(x):
    return x * tan(x)


def deriv(Uh_, Uh):
    return (Uh - Uh_) / (2 * h)


def deriv_seconde(Uh_, U0, Uh, h):
    return (Uh - 2 * U0 + Uh_) / ((h) ** 2)


def richardson(Uh, Uh2):
    return ((4 * Uh2) - Uh) / 3


X = [0.8000, 0.8500, 0.9000, 0.9500, 1.0000]
U = [0.8237, 0.9676, 1.1341, 1.3285, 1.5574]

print(deriv_seconde(U[0], U[2], U[4], 0.1))
print(deriv_seconde(U[1], U[2], U[3], 0.05))
print(
    richardson(
        deriv_seconde(U[0], U[2], U[4], 0.1), deriv_seconde(U[1], U[2], U[3], 0.05)
    )
)


def plot_poly_and_function():
    x = linspace(0, 1, 1000)
    polinterp = Polynomial(lagrange_interp(X[::2], U[::2]).coef[::-1])(x)
    print(Polynomial(lagrange_interp(X[::2], U[::2]).coef[::-1]))
    plt.plot(x, polinterp, "-r", label="lagrange")
    plt.plot(x, x * tan(x), "-b", label="x * tan(x)")
    plt.legend(loc="upper right")
    plt.xlim([0, 1])
    plt.show()
