from cmath import pi
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

pot = np.genfromtxt("data/pot_out.dat")
Nx = 100
Ny = 100
dl = 1
position_space = 0
if position_space==1:
    X = np.arange(-Nx/2, Nx/2, dl)
    Y = np.arange(-Ny/2, Ny/2, dl)
else :
    X = np.arange(-pi/dl, pi/dl, 2*pi/(Nx*dl))
    Y = np.arange(-pi/dl, pi/dl, 2*pi/(Ny*dl))
X, Y = np.meshgrid(X, Y)
Z =  pot #pot[99].reshape(100,100)

fig = plt.figure()
ax = plt.axes()

ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X,Y,Z, cmap ="plasma")
#ax.set_zlim([-2,2])
plt.show()
#ax.scatter(X,Y,Z)

#plt.savefig('img/pot_plot.png')

#FIxer affichage de 0 - > Nx vers -Nx/2 -> Nx/2