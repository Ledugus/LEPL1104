from numpy import *
import matplotlib.pyplot as plt

x = linspace(-3, 3, 700)
y = linspace(-3, 3, 700)
X, Y = meshgrid(x, y)
print(X, Y)
result = zeros(X.shape, dtype=complex)
result.real = X
result.imag = Y
for value in result:
    for i in range(len(value)):
        if absolute(1 + value[i] + value[i] ** 2 / 2) > 1:
            value[i] = 0
plt.plot(result.real, result.imag, "g.")
plt.axis("equal")
plt.axis(True)
plt.xlabel("Re")
plt.ylabel("Im")
plt.xlim(-3, 3)
plt.ylim(-3, 3)
plt.show()
