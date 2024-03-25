from numpy import *

Xstart = 0
Xend = 1
Ustart = 0
Uend = 0.5 * (exp(Xend) - sin(Xend) - cos(Xend))
Eexpl = zeros(8)
Eimpl = zeros(8)
for j in range(8):
    n = pow(2, j + 1)
    h = (Xend - Xstart) / n
    X = linspace(Xstart, Xend, n + 1)
    Uexpl = zeros(n + 1)
    Uexpl[0] = Ustart
    Uimpl = zeros(n + 1)
    Uimpl[0] = Ustart
    for i in range(n):
        Uexpl[i + 1] = Uexpl[i] + h * (sin(X[i]) + Uexpl[i])
        Uimpl[i + 1] = (Uimpl[i] + h * sin(X[i + 1])) / (1 - h)
    Eexpl[j] = abs(Uexpl[-1] - Uend)
    Eimpl[j] = abs(Uimpl[-1] - Uend)
Oexpl = log(abs(Eexpl[:-1] / Eexpl[1:])) / log(2)
Oimpl = log(abs(Eimpl[:-1] / Eimpl[1:])) / log(2)
print("orderExpl ", *["%.4f " % val for val in Oexpl])
print("orderImpl ", *["%.4f " % val for val in Oimpl])
