from numpy import *
from numpy.linalg import solve


def edp_euler(m, n, beta):

    x = linspace(0, 1, m)
    x = linspace(0, 1, m)
    X, Y = meshgrid(x, x)

    U = X * X + Y * Y
    first = copy(U)
    U_loop = copy(U)
    U_new = copy(U)
    for _ in range(n):
        U[1:-1, 1:-1] += beta * (U[2:, 2:] + U[:-2, :-2] - U[:-2, 2:] - U[2:, :-2])

    for _ in range(n):
        for i in range(1, m - 1):
            for j in range(1, m - 1):
                U_new[i, j] += beta * (
                    U_loop[i + 1, j + 1]
                    + U_loop[i - 1, j - 1]
                    - U_loop[i + 1, j - 1]
                    - U_loop[i - 1, j + 1]
                )
        U_loop = copy(U_new)
    import matplotlib.pyplot as plt
    import matplotlib

    matplotlib.rcParams["toolbar"] = "None"
    plt.rcParams["figure.facecolor"] = "lavender"
    plt.rcParams["axes.facecolor"] = "lavender"
    myColorMap = matplotlib.cm.jet

    plt.figure("Python as Matlab clone...")
    plt.contourf(x, x, first, 10, cmap=myColorMap)
    plt.contour(x, x, first, 10, colors="k", linewidths=1)
    plt.axis("off")
    plt.axis("equal")
    plt.show()

    plt.figure("Python as Matlab clone...")
    plt.contourf(x, x, U, 10, cmap=myColorMap)
    plt.contour(x, x, U, 10, colors="k", linewidths=1)
    plt.axis("off")
    plt.axis("equal")
    plt.show()

    plt.figure("Python as Matlab clone...")
    plt.contourf(x, x, U_loop, 10, cmap=myColorMap)
    plt.contour(x, x, U_loop, 10, colors="k", linewidths=1)
    plt.axis("off")
    plt.axis("equal")
    plt.show()


edp_euler(40, 40, 1)
