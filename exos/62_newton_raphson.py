from numpy import *


def f(x):
    return x * exp(x)


def df(x):
    return exp(x) + x * exp(x)


def newton_raphson(f, df, x0, tol=1.0e-30, max_iter=100):
    x = x0
    n_iter = 0
    for i in range(max_iter):
        n_iter += 1
        x_new = x - f(x) / df(x)
        if abs(x_new - x) < tol:
            return x_new, n_iter
        x = x_new
    return x, n_iter


print(newton_raphson(f, df, 0.2))
print(newton_raphson(f, df, 20))
