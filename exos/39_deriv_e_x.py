from numpy import *


def f(x):
    return exp(x)


def deriv(f, x, h):
    """
    Calcule la dérivée d'une fonction f(x) en utilisant la méthode de différences finies
    f'(x) = (f(x+h) - f(x-h)) / (2h)
    """
    return (f(x + h) - f(x - h)) / (2 * h)


def richardson(f, x, h):
    """
    Calcule la dérivée d'une fonction f(x) en utilisant la méthode de Richardson
    """
    return (4 * deriv(f, x, h / 2) - deriv(f, x, h)) / 3


def error(f, x, h):
    """
    Calcule l'erreur de la méthode de Richardson
    """
    return abs(deriv(f, x, h) - f(x))


def troncature(f, x, h):
    return h * h * f(x + h) / 6


def main():
    x = 2
    h = 0.1
    print("exp(2) : ", exp(x))
    print("h=0.1")
    print(deriv(f, x, h))
    print("error : ", deriv(f, x, h) - f(x))
    print("troncature théorique: ", troncature(f, x, h))
    h = 0.01
    print("h=0.01")
    print(deriv(f, x, h))
    print("error : ", deriv(f, x, h) - f(x))
    print("troncature théorique: ", troncature(f, x, h))
    h = 0.001
    print("h=0.001")
    print(deriv(f, x, h))
    print("error : ", deriv(f, x, h) - f(x))
    print("troncature théorique: ", troncature(f, x, h))


if __name__ == "__main__":
    main()
