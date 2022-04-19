import sys
from csv import writer
from ctypes import sizeof
import matplotlib.pyplot as plt
import matplotlib.animation as animate
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
def main(argv):
    x_size = int(argv[0])
    y_size = int(argv[1])
    V0 = float(argv[2])
    corr_len = float(argv[3])
    strV0 = ''
    strCorr_len = ''
    if V0 == int(V0):
        strV0 = str(int(V0))
    else:
        strV0 = str(V0)
    
    if corr_len == int(corr_len):
        strCorr_len = str(int(corr_len))
    else:
        strCorr_len = str(corr_len)
    dataSpec = str(x_size) + "x" + str(y_size) + "_V0_" + strV0 + "_" + strCorr_len

    Writer = animate.writers['ffmpeg']
    writer = Writer(fps=25, metadata=dict(artist='Me'), bitrate=2**16)
    fps = 24

    data = np.genfromtxt("data/" + dataSpec + "/Mean/mean.dat")
    frn = len(data)

    X = np.arange(0, x_size, 1)
    Y = np.arange(0, y_size, 1)
    X, Y = np.meshgrid(X, Y)
    Z = np.zeros((x_size, y_size, frn))

    for i in range(len(data)):
        Z[:, :, i] = data[i].reshape(x_size, y_size)


    def change_plot(frame_number, zarray, plot):
        plot[0].remove()
        plot[0] = ax.plot_surface(X, Y, zarray[:, :, frame_number], cmap="plasma")


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    plot = [ax.plot_surface(X, Y, Z[:, :, 0], color='0.75', rstride=1, cstride=1)]

    ax.set_zlim(0, 1.1)
    ani = animate.FuncAnimation(fig, change_plot, frn,
                                fargs=(Z, plot), interval=1000 / fps)

    ax.axis('off')
    ani.save('img/im.mp4', writer=writer)

    plt.show()

if __name__ == "__main__":
   main(sys.argv[1:])
