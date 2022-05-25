from ast import Str
from json.tool import main
import math
import os
import sys

import matplotlib.pyplot as plt
import numpy as np


def main(argv):
    x_size = int(argv[0])
    y_size = int(argv[1])
    V0 = float (argv[2])
    corr_len = float(argv[3])

    # name and path of the folder
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
    dataSpec = str(x_size) + "x" + str(y_size) + \
        "_V0_" + strV0 + "_" + strCorr_len

    dataPath = "data/" + dataSpec + "/"
    # Check how many files there are in the folder (excludes the debug file)
    files = [name for name in os.listdir(dataPath) if (
        os.path.isfile(dataPath + name) and name != 'debug.dat')]

    # number of frames per file/simulation
    frn = len(np.genfromtxt(dataPath + files[0], dtype=np.complex128))
    fileProc = np.zeros((x_size * y_size, frn, len(files)),
                        dtype=np.complex128)  # processed files
    i = 0
    for file in files:
        fileData = np.genfromtxt(dataPath + file, dtype=np.complex128)
        for j in range(frn):
            # i is the file number, j the frame and : is the data
            fileProc[:, j, i] = fileData[j]
        i += 1

    # average wavefunction of the $frn simulations
    dataMean = np.mean(fileProc, axis=2)
    s = 0
    for i in range(x_size):
        for j in range(y_size):
            s += abs(dataMean[i*y_size + j,0])**2 * (2*np.pi/x_size)*(2*np.pi/y_size)

    if s==1:
        print('Mean wavefunction is normalized')
    else:
        print('ERROR: The mean wavefunction is not normalized = ' + str(s))

    try:
        os.mkdir(dataPath + "/Mean")
    except OSError as error:
        print(error)

    f = open(dataPath + "Mean/mean.dat", "w")
    np.savetxt(f, dataMean)
    f.close()


if __name__ == "__main__":
    main(sys.argv[1:])
