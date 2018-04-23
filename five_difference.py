import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

class difference():
    def LaxFriedrichs(self,N,M,C,r,U,P,E,x,t):

        if abs(C*r) > 1:
            print('|C*r|>1,LaxFriedrichs差分格式不稳定!')

        #solution by Layer
        for j in range(N):
            for i in range(1,M):
                U[i,j+1] = (U[i+1,j]+U[i-1,j])/2 - C*r*(U[i+1,j]-U[i-1,j])/2
                P[i,j+1] = np.cos(np.pi*(x[i]+t[j+1]))
                c = abs(U[i,j+1]-np.cos(np.pi*(x[i]+t[j+1])))
                E[i,j+1] = c
        return P,U,E
    def LaxWendoroff(self,N,M,C,r,U,P,E,x,t):
        if abs(C*r) > 1:
            print('|C*r|>1,LaxWendoroff差分格式不稳定!')

        for j in range(N):
            for i in range(1,M):
                U[i,j+1] = U[i,j] - C*r*(U[i+1,j]-U[i-1,j])/2+C**2*r**2*(U[i+1,j] - 2*U[i,j]+U[i-1,j])/2
                P[i,j+1] = np.cos(np.pi*(x[i]+t[j+1]))
                c = abs(U[i,j+1]-np.cos(np.pi*(x[i]+t[j+1])))
                E[i,j+1] = c
        return P,U,E

def PDEHyperbolic(uX,uT,M,N,C,dtype = 'LaxFriedrichs'):

    h = uX/M#变量x的步长
    k = uT/N#变量t的步长
    r = k/h #步长比
    x = np.arange(0,1+h/2,h)
    t = np.arange(0,1+k/2,k)
    U = np.zeros((M+1,N+1))#代数解
    P = np.zeros((M+1,N+1))#真解
    E = np.zeros((M+1,N+1))#误差

    #init_value
    U[:,0] = np.cos(np.pi*x)
    P[:,0] = np.cos(np.pi*x)
    E[:,0] = 0

    #boundary value condition

    U[0,:] = np.cos(np.pi*t)
    E[0,:] = 0
    P[0,:] = np.cos(np.pi*t)
    U[M,:] = -np.cos(np.pi*t)
    P[M,:] = -np.cos(np.pi*t)
    E[M,:] = 0
    
    test = difference()
    if dtype is 'LaxFriedrichs':
        P,U,E = test.LaxFriedrichs(N,M,C,r,U,P,E,x,t)
 
    elif dtype is 'LaxWendoroff':
        P,U,E = test.LaxWendoroff(N,M,C,r,U,P,E,x,t)
    else:
        raise 



    return U,P,E,x,t

uX = 1
uT = 1
M = 90
N = 100
C = -1
P,U,E,x,t = PDEHyperbolic(uX,uT,M,N,C,dtype = 'LaxWendoroff')
U = U.T
print(U)
P = P.T
print(P)
fig = plt.figure()
ax = Axes3D(fig)
x,t = np.meshgrid(x,t)
print(np.shape(x),np.shape(t),np.shape(U))
ax.plot_surface(x,t,U,rstride=10,cstride=10,cmap = cm.viridis)
ax.plot_surface(x,t,U,rstride=10,cstride=10,cmap = cm.viridis)
plt.xlabel('x')
plt.ylabel('t')
ax.set_zlabel('U(x,t)')
plt.show()
