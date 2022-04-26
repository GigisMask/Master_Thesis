import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

pot = np.genfromtxt("data/100x100_V0_0.3_10/Mean/mean.dat")

X = np.arange(0, 100, 1)
Y = X
X, Y = np.meshgrid(X, Y)
Z =  pot[99].reshape(100,100)

fig = plt.figure()
ax = plt.axes()

ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X,Y,Z, cmap ="plasma")
#ax.set_zlim([-2,2])
plt.show()
#ax.scatter(X,Y,Z)

#plt.savefig('img/pot_plot.png')

#FIxer affichage de 0 - > Nx vers -Nx/2 -> Nx/2