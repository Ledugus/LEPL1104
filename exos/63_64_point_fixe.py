from numpy import *


def g(x):
    return 2 * sqrt(x - 1)


def g_2(x):
    return x**2 / 4 + 1


def g_3(x):
    try:
        return 4 * x - x**2 / 2 - 4
    except:
        print("Erreur de calcul")
        return x


def point_fixe(g, x0, tol=1.0e-4, max_iter=20):
    x = x0
    n_iter = 0
    for _ in range(max_iter):
        n_iter += 1
        x_new = g(x)
        print(x_new)
        if abs(x_new - x) < tol:
            return x_new, n_iter
        x = x_new
    return x, n_iter


# 63
"""
print("Point fixe problématique")
print(point_fixe(g, 1.5))
print("Point fixe OK")
print(point_fixe(g, 2.5))
print("Point fixe OK")
print(point_fixe(g_2, 1.5))
print("Point fixe problématique")
print(point_fixe(g_2, 2.5))
"""

# 64
print("Point fixe OK")
print(point_fixe(g_3, 1.9))
print("Point fixe OK")
print(point_fixe(g_3, 3.8))
print("Point fixe OK")
print(point_fixe(g_3, 6))
