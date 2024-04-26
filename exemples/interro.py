from numpy import *
import matplotlib.pyplot as plt


def bezier_curves(t):
    return (1 - t) ** 3, 3 * (1 - t) * (1 - t) * t, 3 * (1 - t) * t * t, t**3


alpha = 5 / 3
t = linspace(0, 1, 100)
t_croisement = searchsorted(t, sqrt(3) / 2)
print(t_croisement)
b0, b1, b2, b3 = bezier_curves(t)
plt.plot(t, b0)
plt.plot(t, b1)
plt.plot(t, b2)
plt.plot(t, b3)
plt.title("Courbes de Bézier")
x_a = y_a = 3 * b0 + 2 * b1 + alpha * b2
x_s = -b1 + b3
y_s = 2 * b2
plt.figure()
plt.title("Trajectoires")
plt.plot(x_a, y_a)
plt.plot(x_s, y_s)
plt.plot(x_s[t_croisement], y_s[t_croisement], ".r", markersize=10)
plt.plot(x_a[t_croisement], y_a[t_croisement], ".r", markersize=10)
plt.show()
