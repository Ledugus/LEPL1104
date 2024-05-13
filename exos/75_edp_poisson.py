from numpy import *
from numpy.linalg import solve
from scipy.sparse import dok_matrix
from scipy.sparse.linalg import spsolve


def poissonSolve(nx, ny):
    m = nx * ny
    h = 2 / (nx - 1)
    A = dok_matrix((m, m), dtype=float32)
    B = zeros(m)
    for i in range(m):
        A[i, i] = 1.0
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            index = i + j * nx
            A[index, index] = 4
            A[index, index - 1] = -1
            A[index, index + 1] = -1
            A[index, index + nx] = -1
            A[index, index - nx] = -1
            B[index] = 1

    U = spsolve((A / (h * h)).tocsr(), B).reshape(ny, nx)

    import matplotlib.pyplot as plt
    import matplotlib

    matplotlib.rcParams["toolbar"] = "None"
    plt.rcParams["figure.facecolor"] = "lavender"
    plt.rcParams["axes.facecolor"] = "lavender"
    myColorMap = matplotlib.cm.jet
    x = linspace(-1, 1, nx)

    y = linspace(-1, 1, ny)

    plt.figure("Python as Matlab clone...")
    plt.contourf(x, y, U, 10, cmap=myColorMap)
    plt.contour(x, y, U, 10, colors="k", linewidths=1)
    plt.axis("off")
    plt.axis("equal")
    plt.show()

    return U


I_0 = poissonSolve(4, 4)
I = poissonSolve(10, 10)
I_1 = poissonSolve(20, 20)
I_2 = poissonSolve(40, 40)
