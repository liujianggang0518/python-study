import numpy as np

box = [0,1,0,1]
nx = 10
ny = 10
N = (nx+1)*(ny+1)
NC = nx*ny
node = np.zeros((N,2))
X, Y = np.mgrid[box[0]:box[1]:complex(0,nx+1),box[2]:box[3]:complex(0,ny+1)]#return mesh-grid ndarrays all of the same dimensions
node[:,0] = X.flatten()
node[:,1] = Y.flatten()
cell = np.zeros((2*NC,3),dtype = np.int)
sign = np.arange(N).reshape(nx+1,ny+1)#网格每个点的标号，有nx+1列
cell[:NC,0] = sign[1:,1:].flatten()
cell[:NC,1] = sign[1:,1:].flatten()
cell[:NC,2] = sign[0:-1,0:-1].flatten()#每个小正方形下面的小三角形
cell[NC:,0] = sign[0:-1,1:].flatten()
cell[NC:,1] = sign[0:-1,0:-1].flatten()
cell[NC:,2] = sign[1:,1:].flatten()
