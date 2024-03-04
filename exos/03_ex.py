import math
import numpy as np


def get_closest_ints(t: float, n: int, k=4):
    if t >= n - 3:
        return list(range(n - 3, n + 1))
    if t <= -n + 3:
        return list(range(-n, -n + 4, 1))
    if t.is_integer():
        t = int(t)
        return list(range(t - 1, t + 3, 1))
    return list(range(math.floor(t) - 1, math.ceil(t) + 2))


def interpolate(t, u, n):
    absisses = get_closest_ints(t, n)
    ordonnees = [u[x] for x in absisses]
    coefs = np.polyfit(absisses, ordonnees, 3)
    result = 0
    for i in range(4):
        result += t**i * coefs[i]
    return result, coefs


def compute(n, u, t):
    #
    # -1- Calcul du noeud qui precede Xestim
    # Gestion des cas critiques (extrapolation `a droite et `a gauche)
    #
    k = int(math.floor(t))
    k = n - 2 if (k + 2 > n) else k
    k = -n + 1 if (k - 1 < -n) else k
    #
    # -2- Estimation avec des polynomes de Lagrange
    # xi = coordonnee locale avec Xlocal = [-1 0 1 2]
    #
    u_local = u[n + k - 1 : n + k + 3]
    print(k)
    print("échantillon", u_local)
    xi = t - k
    phi = (
        np.array(
            [
                -xi * (xi - 1) * (xi - 2),
                3 * (xi + 1) * (xi - 1) * (xi - 2),
                -3 * (xi + 1) * xi * (xi - 2),
                (xi + 1) * xi * (xi - 1),
            ]
        )
        / 6
    )

    print(u_local.shape())
    print(phi[0].shape())
    return u_local @ phi


print(compute(8, np.array([0, 1, 2, 3, 4, 5, 6, 7, 8]), 1))
