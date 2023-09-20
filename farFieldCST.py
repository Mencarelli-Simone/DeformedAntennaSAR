# Simone Mencarelli
# September 2023
# This file contains a class with the same interface of the Aperture class in radartools.farField
# The pattern however is loaded from a CST ffs file and provided for any theta phi coordinate
# by means of an interpolator, namely the sphere_interp in interpolator_v2.py

# %% includes
import numpy as np
from interpolator_v2 import sphere_interp


# %% functions

# far field source loader (ffsFileReader script)
def ffsLoader(filename):
    # %% variables to read from file
    # header
    frequencies = 0
    position = [0, 0, 0]
    radiatedPower = 0
    acceptedPower = 0
    stimulatedPower = 0
    frequency = 0
    phiSamples = 0
    thetaSamples = 0

    # data matrix
    dataMatrix = None  # size is going to be initialized after reading header

    # %% Read file in ram
    with open(filename, 'r') as file:
        content = file.read()
    content = content.split('\n')

    # %% parse
    # parse just header
    for i in range(len(content)):
        if "Frequencies" in content[i]:
            frequencies = int(content[i + 1])
            i += 1
        if "Position" in content[i]:
            pos = content[i + 1]
            pos = pos.split(' ')
            position = np.array(pos[0:3], dtype=float).reshape((3, 1))
            i += 1
        if "Radiated/Accepted/Stimulated Power , Frequency" in content[i]:
            radiatedPower = float(content[i + 1])
            acceptedPower = float(content[i + 2])
            stimulatedPower = float(content[i + 3])
            frequency = float(content[i + 4])
        if "Total #phi samples, total #theta samples" in content[i]:
            sam = content[i + 1]
            sam = sam.split(' ')
            phiSamples = int(sam[0])
            thetaSamples = int(sam[1])
            break

    # %% data parsing
    # Phi, Theta, Re(E_Theta), Im(E_Theta), Re(E_Phi), Im(E_Phi)
    dataMatrix = np.zeros((phiSamples * thetaSamples, 6))

    for i in range(len(content)):
        if "Phi, Theta, Re(E_Theta), Im(E_Theta), Re(E_Phi), Im(E_Phi)" in content[i]:
            for j in range(phiSamples * thetaSamples):
                dataline = content[i + j + 1]
                dataline = dataline.split(" ")
                while "" in dataline:
                    dataline.remove("")
                dataline = np.array(dataline[0:6], dtype=float).reshape((1, 6))
                dataMatrix[j, :] = dataline

    # %% turn the data into meshgrids
    # theta component
    E_Theta = dataMatrix[:, 2] + 1j * dataMatrix[:, 3]
    E_Theta = E_Theta.reshape((phiSamples, thetaSamples))
    # phi component
    E_Phi = dataMatrix[:, 4] + 1j * dataMatrix[:, 5]
    E_Phi = E_Phi.reshape((phiSamples, thetaSamples))
    # coordinates
    Phi = dataMatrix[:, 0]
    Phi = Phi.reshape((phiSamples, thetaSamples))
    Theta = dataMatrix[:, 1]
    Theta = Theta.reshape((phiSamples, thetaSamples))

    return Phi, Theta, E_Phi, E_Theta, phiSamples, thetaSamples, radiatedPower, stimulatedPower, acceptedPower


# Aperture class for interfacing cst pattern
class Aperture:
    def init(self, filename):
        """
        initialization method, it requires a far field file
        :param filename: CST ffs file
        :return:
        """
        # load far field
        (Phi, Theta, E_Phi, E_Theta, phiSamples, thetaSamples,
         radiatedPower, stimulatedPower, acceptedPower) = ffsLoader(filename)

        # compute directive gain
        self.G = 2 * np.pi * (np.abs(E_Theta) ** 2 + np.abs(E_Phi) ** 2) / (120 * np.pi * radiatedPower)

        # store relevant parameters
        self.Theta = Theta * np.pi / 180  # I use radians
        self.Phi = Phi * np.pi / 180
        self.phiSamples = phiSamples
        self.thetaSamples = thetaSamples

    def mesh_gain_pattern(self, theta_mesh: np.ndarray, phi_mesh: np.ndarray, cubic=True):
        """
        retruns the gain pattern at the specified meshgrid points in spherical coordinates.
        :param theta_mesh: Theta coordinates
        :param phi_mesh: Phi coordinates
        :param cubic: default True: bicubic interpolation utilised, False: linear interpolation.
        :return:
        """
        # need to create an outpattern for some reason
        outpattern = np.zeros_like(theta_mesh).reshape(-1)
        # The vectors need to be flattened hence reshape(-1)
        outpattern = sphere_interp(theta_mesh.reshape(-1), phi_mesh.reshape(-1),
                                   self.Theta.reshape(-1), self.Phi.reshape(-1),
                                   self.G.reshape(-1), outpattern, cubic)
        # reshape to desired format
        return outpattern.reshape(np.shape(theta_mesh))

    # %% todo testing
