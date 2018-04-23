from math import sqrt
import numpy as np
from scipy import linalg
from matplotlib import pyplot as plt


def f(x):
    return 16 * (np.pi) ** 2 * np.sin(4 * np.pi * x)

def u(x):                #真解函数
    return np.sin(4 * np.pi * x)

#利用中心差分格式求解两点边值问题
def FD1d_bvp(N, f, a, b, u):
    h = (b - a) / (N - 1)
    x = np.linspace(a, b, N, endpoint=True)
    #创建线性差分方程组系数矩阵
    c1 = - 1 / (h ** 2)
    c2 = 2 / (h ** 2)
    g = np.hstack((c1 * np.ones(N - 2), 0))
    c = np.hstack((0, c1 * np.ones(N - 2)))
    d = np.hstack((1, c2 * np.ones(N - 2), 1))
    A = np.diag(g, -1) + np.diag(d) + np.diag(c, 1)
    #创建线性方程组右端项
    rhs = f(x)
    rhs[0] = u(x[0])
    rhs[N - 1] = u(x[N - 1])
    #求解上述代数系统
    U = linalg.solve(A, rhs)
    return x, U


def FD1d_error(x, U, u_exact):
    N = len(x)
    h = (x[-1] - x[0]) / (N - 1)
    ue = u_exact
    ee = ue - U
    e0 = h * np.sum(ee ** 2)
    e1 = sum((ee[1:] - ee[:-1]) ** 2) / h
    e1 = e1 + e0
    e0 = sqrt(e0)
    e1 = sqrt(e1)
    emax = np.max(abs(ue - U))
    return emax, e0, e1

L = 0
R = 1
N = [6, 11, 21, 41, 81]
X = []
UN = []
for i in range(5):
    x, U = FD1d_bvp(N[i], f, L, R, u)
    u_exact = u(x)
    emax, e0, e1 = FD1d_error(x, U, u_exact)
    X.append(x)
    UN.append(U)

ue = u(X[-1])

plt.plot(X[-1], ue, '-k*',
        X[0], UN[0], '-ro',
        X[1], UN[1], '-gs',
        X[2], UN[2], '-bd',
        X[3], UN[3], '-ch',
        X[4], UN[4], '-mx')
plt.title("The solution plot")
plt.xlabel("x")
plt.ylabel("u")
plt.legend(labels=['exact', 'N=6', 'N=11', 'N=21', 'N=41', 'N=81'])
plt.show()

