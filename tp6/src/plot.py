from matplotlib.pyplot import *
from math import *
import numpy as np

with open("plotpy.dat", "r") as f:
    contenu = f.read().split("\n\n")


x = contenu[0].split()
nx = len(x) - 1
x = np.array([float(x[i]) for i in range(nx+1)])

y = contenu[1].split()
ny = len(y) -1
y = np.array([float(y[i]) for i in range(ny+1)])

x , y = np.meshgrid(x,y)
z = contenu[2].split()
nz = len(z)
z = np.array([float(z[i]) for i in range(nz)]).reshape((ny+1,nx+1))

fig, ax = subplots()
cs = ax.contourf(x,y,z,100)
cbar = fig.colorbar(cs)

print("press \'q\' to quit...");
show()
print("The end")