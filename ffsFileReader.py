# Simone Mencarelli
# September 2023
# This script loads in memory an antenna pattern from a CST ffs (far field source) file.
# the content of this file is to be integrated in an antenna class

# %% Includes section
import matplotlib.pyplot as plt
import numpy as np

# %% User input
filename = 'farfield.ffs'

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

# plot
fig, ax = plt.subplots(1)
ax.contourf(Theta, Phi, 10 * np.log10(np.abs(E_Theta ** 2 + E_Phi ** 2)))
plt.show()

# %% we assume the cross polarization has no effect and define the directive gain as

G = 2 * np.pi * (np.abs(E_Theta) ** 2 + np.abs(E_Phi) ** 2) / (120 * np.pi * radiatedPower)
fig, ax = plt.subplots(1)
ax.contourf(Theta, Phi, 10 * np.log10(G))
plt.show()

# make this available for the snr and impulse response
