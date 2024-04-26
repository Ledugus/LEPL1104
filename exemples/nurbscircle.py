#
# Dessiner un cercle avec des NURBS
# Vincent Legat - 2022
# Ecole Polytechnique de Louvain
#

from numpy import *
import matplotlib
import matplotlib.pyplot as plt


# =========================================================================


def b(t, T, i, p):
    if p == 0:
        return (T[i] <= t) * (t < T[i + 1])
    else:
        u = (
            0.0
            if T[i + p] == T[i]
            else (t - T[i]) / (T[i + p] - T[i]) * b(t, T, i, p - 1)
        )
        u += (
            0.0
            if T[i + p + 1] == T[i + 1]
            else (T[i + p + 1] - t) / (T[i + p + 1] - T[i + 1]) * b(t, T, i + 1, p - 1)
        )
        return u


# ============================= mainProgram ===============================


T = [0, 0, 0, 1, 1, 2, 2, 3, 3, 4]
a = sqrt(3)
X = [1, 0, 1 / 2, 1, 3 / 2, 2, 1]
Y = [0, 0, a / 2, a, a / 2, 0, 0]
W = [1, 0.5, 1, 0.5, 1, 0.5, 1]

p = 2
n = len(T) - 1
W = array(W)
t = arange(T[p], T[n - p], 0.001)
B = zeros((n - p, len(t)))
for i in range(0, n - p):
    B[i, :] = b(t, T, i, p)
w = W @ B
x = W * X @ B / w
y = W * Y @ B / w
print(W * X)
print(W * X @ B - X * W @ B)

# matplotlib.rcParams["toolbar"] = "None"
matplotlib.rcParams["lines.linewidth"] = 1
plt.figure("Un joli cercle")
plt.axis("off")
plt.axis("equal")
plt.plot(X, Y, "--g")
plt.plot(X, Y, "or", markersize=5)
plt.plot(x, y, "-b")
plt.plot(x[: len(x) // 2], y[: len(y) // 2], "-y")
plt.show()
