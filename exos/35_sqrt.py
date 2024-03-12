from numpy import *


def integrate_gauss_legendre(f, a, b):
    """Retourne une approximation de l'intégrale de 1/sqrt(x) entre 0 et 1
    en utilisant la méthode de Gauss-Legendre avec 3 points
    """
    x = array([0.774596669241483377035853079956, 0, -0.774596669241483377035853079956])
    c = array(
        [
            0.555555555555555555555555555556,
            0.888888888888888888888888888889,
            0.555555555555555555555555555556,
        ]
    )

    x = (b - a) / 2 * x + (a + b) / 2
    return (b - a) / 2 * sum(c * f(x))


def integrate_gauss_legendre_5(f, a, b):
    """Retourne une approximation de l'intégrale de 1/sqrt(x) entre 0 et 1
    en utilisant la méthode de Gauss-Legendre avec 5 points
    """
    x = array(
        [
            0.906179845938663992797626878299,
            0.538469310105683091036314420700,
            0,
            -0.538469310105683091036314420700,
            -0.906179845938663992797626878299,
        ]
    )
    c = array(
        [
            0.236926885056189087514264040720,
            0.478628670499366468041291514836,
            0.568888888888888888888888888889,
            0.478628670499366468041291514836,
            0.236926885056189087514264040720,
        ]
    )

    x = (b - a) / 2 * x + (a + b) / 2
    return (b - a) / 2 * sum(c * f(x))


print(integrate_gauss_legendre(lambda x: 1 / sqrt(x), 0, 1))
print(integrate_gauss_legendre_5(lambda x: 1 / sqrt(x), 0, 1))
