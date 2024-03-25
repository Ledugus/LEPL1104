#
# Explicit Euler and Taylor-4 methods
# Vincent Legat - 2018
# Ecole Polytechnique de Louvain
#

from numpy import *
from matplotlib import pyplot as plt

u = lambda x: 3 * exp(-x / 2) - 2 + x

Xstart = 0
Xend = 3
Ustart = 1

# =================== Explicit Euler =========

print(" ====== Explicit Euler method =====")
error = zeros(4)
for j in range(4):
    n = 3 * pow(2, j)
    h = (Xend - Xstart) / n
    X = linspace(Xstart, Xend, n + 1)
    U = zeros(n + 1)
    U[0] = Ustart
    for i in range(n):
        U[i + 1] = U[i] + h * (X[i] - U[i]) / 2
    plt.plot(X, U, ".-b")
    error[j] = abs(U[-1] - u(Xend))
    print(" ==== Euler  (order=1) h=%5.3f : eh(Xend) = %8.2e " % (h, error[j]))
order = mean(log(error[:-1] / error[1:]) / log(2))
print(" ============= Estimated order : %.4f " % order)


# =================== Taylor order 4 =========

print(" ====== Explicit Taylor order 4 method =====")
error = zeros(4)
for j in range(4):
    n = 3 * pow(2, j)
    h = (Xend - Xstart) / n
    X = linspace(Xstart, Xend, n + 1)
    U = zeros(n + 1)
    U[0] = Ustart
    for i in range(n):
        U[i + 1] = (
            U[i]
            + h * (X[i] - U[i]) / 2
            + (h**2) * (2 - X[i] + U[i]) / 8
            - (h**3) * (2 - X[i] + U[i]) / 48
            + (h**4) * (2 - X[i] + U[i]) / 384
        )
    plt.plot(X, U, ".-r")
    error[j] = abs(U[-1] - u(Xend))
    print(" ==== Taylor (order=4) h=%5.3f : eh(Xend) = %8.2e " % (h, error[j]))

order = mean(log(error[:-1] / error[1:]) / log(2))
print(" ============= Estimated order : %.4f " % order)

plt.show()
