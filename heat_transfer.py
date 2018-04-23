import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

def time_grid(b,NT):
    T = np.linspace(0,b,NT+1)
    tau = b/NT
    return T,tau

def space_grid(a,NS):
    X = np.linspace(0,a,NS+1)
    h = a/NS
    return X,h

def u_initial(x):
    u = np.exp(-(x - 0.25)**2/0.01) + 0.1*np.sin(20*np.pi*x)
    return u

def f_right(X,T):
    f = np.zeros(len(X))
    return f

def forward(N,M,r,tau,X,T,U):

    d = 1 - 2*np.ones(N-2)*r
    c = np.ones(N-3)*r
    A = np.diag(c,-1) + np.diag(d,0) + np.diag(c,1)

    for i in range(1,M):
        rhs = tau*f_right(X,T[i])
        rhs[1] = rhs[1] + r*U[0,i-1]
        rhs[-2] = rhs[-2] + r*U[-1,i-1]
        U[1:-1,i] = np.dot(A,U[1: -1,i-1]) + rhs[1:-1]
    
    return U

def backward(N,M,r,tau,X,T,U):
    d = 1+2*np.ones(N-2)*r
    c = -np.ones(N-3)*r
    A = np.diag(c,-1) + np.diag(d) + np.diag(c,1)
 
    for i in range(1,M): 
        rhs = tau*f_right(X,T[i])
        rhs[1] = rhs[1] + r*U[0,i]
        rhs[-2] = rhs[-2] + r*U[-1,i-1]
        A = np.linalg.inv(A)
        U[1:-1,i] = np.dot(A,(U[1:-1,i-1]+rhs[1:-1]))
    return U

def crank_nicholson(N,M,r,tau,X,T,U):
    d1 = 1+2*np.ones(N-2)*r
    d2 = 1 - 2*np.ones(N-2)*r
    c = 0.5*np.ones(N-3)*r

    A1 = np.diag(-c,-1) + np.diag(d1) + np.diag(-c,1)
    A0 = np.diag(c,-1) + np.diag(d2) + np.diag(c,1)
    
    for i in range(1,M):
        rhs = tau*f[:,T[i]]
        rhs[1] = rhs[1] + 0.5*r*(U[0,i]+U[0:i-1])
        rhs[-2] = rhs[-2] + 0.5*r*(U[-1,i]+U[-1,i-1])
        A1 = np.linalg.inv(A1)
        U[1:-1,i] = np.dot(A1,A0*U[1:-1,i-1]+rhs[1:-1])
        
    return U

def heat_equation_fd1d(NS,NT,a,b,k,dtype = 'forward'):
    X,h = space_grid(a,NS)
    T,tau = time_grid(b,NT)
    N = len(X)
    M = len(T) 
    
    r =k*tau/(h**2)
    if r >= 0.5 and (dtype is 'forward'):
        print('The time error does not satisfy the stability!')
    U = np.zeros((N,M))
    U[:,0] = u_initial(X)
    print(U[:,0])
    U[0,:] = np.zeros(len(T))
    U[-1,:] = np.zeros(len(T))
    if dtype is 'forward':
        U = forward(N,M,r,tau,X,T,U)

    elif dtype is 'backward':

        U = backward(N,M,r,tau,X,T,U)
    elif dtype is 'crank_nicholson':

        U = crank_nicholson(N,M,r,tau,X,T,U)
    else:
        exit()

    
    return X,T,U     


a = 1
b = 0.1
k = 1
NS = 100
NT = 10000

X,T,U = heat_equation_fd1d(NS,NT,a,b,k,dtype = 'forward')

U = U.T

fig = plt.figure()
ax = Axes3D(fig)
x,t = np.meshgrid(X,T)
ax.plot_surface(x,t,U,rstride=10,cstride=10,cmap = cm.viridis)
plt.xlabel('X')
plt.ylabel('T')
ax.set_zlabel('U(X,T)')
#print('X,T',x,t)
#m = int(0.1/dt +1)
#n = int(1/dr+1)
#for i in range(m):
#    for j in range(n):
#        z = U[i,j]
#ax.plot_surface(X,T,z,rstride = 1,cstride = 1,cmap='rainbow')
plt.show()
