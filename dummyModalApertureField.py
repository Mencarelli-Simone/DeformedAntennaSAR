# Simone Mencarelli
# October 23
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.interpolate import griddata


def aperture_distribution(filename, wavelength):
    """

    """
    rows = 3312
    columns = 11

    # %% reader
    with open(filename, 'r') as file:
        content = file.read()
    content = content.split('\n')

    dataMatrix = np.zeros((rows, columns))

    for i in range(len(content)):
        if "Freq" in content[i]:
            for j in range(rows):
                dataline = content[i + j + 1]
                dataline = dataline.split(" ")
                if len(dataline) < columns:  # eof
                    print('eof')
                    break
                else:
                    while "" in dataline:
                        dataline.remove("")
                    dataline = np.array(dataline[0:columns], dtype=float).reshape((1, columns))
                    dataMatrix[j, :] = dataline.astype('float')
            print('eof')
    # %%
    # raw line vectors
    Xundef = dataMatrix[:, 2]
    Yundef = dataMatrix[:, 3]
    Zundef = dataMatrix[:, 4]
    Xdef = dataMatrix[:, 5]
    Ydef = dataMatrix[:, 6]
    Zdef = dataMatrix[:, 7]
    # %% reconstructing the meshgrid
    xpoints = 23
    ypoints = 144

    Xundef = Xundef.reshape((ypoints, xpoints))  # actually y
    Yundef = Yundef.reshape((ypoints, xpoints))  # actually z
    Zundef = Zundef.reshape((ypoints, xpoints))  # actually x
    Xdef = Xdef.reshape((ypoints, xpoints))
    Ydef = Ydef.reshape((ypoints, xpoints))
    Zdef = Zdef.reshape((ypoints, xpoints))
    # %% plotter
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot_surface(Xundef, Yundef, Zundef)
    plt.show()
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot_surface(Xdef, Ydef, Zdef)
    plt.show()

    # %% resample the deformed xy grid to a uniform grid
    # data to resample and y is z, x is y, and z is x

    deformation = Ydef - Yundef
    xlim = [np.min(Zdef), np.max(Zdef)]
    ylim = [np.min(Xdef), np.max(Xdef)]

    # resample
    # number of points 1/3 wavelength
    no = np.round((xlim[1] - xlim[0]) / (wavelength / 3)).astype('int')
    if no % 2 == 0:
        no += 1
    x = np.linspace(xlim[0], xlim[1], no)  # before .2 is constrained
    no = np.round((ylim[1] - ylim[0]) / (wavelength / 3)).astype('int')
    if no % 2 == 0:
        no += 1
    y = np.linspace(0.2, ylim[1], no)
    X, Y = np.meshgrid(x, y)
    Z = griddata((Zdef.reshape(-1), Xdef.reshape(-1)), deformation.reshape(-1), (X, Y))
    Z[np.isnan(Z)] = 0

    # amplitude mask
    A = griddata((Zdef.reshape(-1), Xdef.reshape(-1)), np.ones_like(Xdef.reshape(-1)),
                 (X, Y), 'nearest')
    A[np.isnan(Z)] = 0

    E = A * np.exp(1j * Z / wavelength * 2 * np.pi)

    return x, y, E
