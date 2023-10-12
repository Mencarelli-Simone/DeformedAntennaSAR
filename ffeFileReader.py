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
filename = 'lyceanem/DeformedAntennav2.ffe'
# header
frequencies = 0
position = [0, 0, 0]
radiatedPower = 1
acceptedPower = 1
stimulatedPower = 1
frequency = 0
phiSamples = 1
thetaSamples = 1

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
    if "#Frequency: " in content[i]:
        id = content[i].find(": ")
        frequency = float(content[i][id + 2:])
    if "#No. of Theta Samples: " in content[i]:
        id = content[i].find(": ")
        thetaSamples = int(content[i][id + 2:])
    if "#No. of Phi Samples: " in content[i]:
        id = content[i].find(": ")
        phiSamples = int(content[i][id + 2:])
        break

# %% data parsing
# Theta, Phi, Re(E_Theta), Im(E_Theta), Re(E_Phi), Im(E_Phi)
dataMatrix = np.zeros((phiSamples * thetaSamples, 9))

for i in range(len(content)):
    if "Re(Etheta)" in content[i]:
        for j in range(phiSamples * thetaSamples):
            dataline = content[i + j + 1]
            dataline = dataline.split(" ")
            if len(dataline) < 9:  # eof
                print('eof')
                break
            else:
                while "" in dataline:
                    dataline.remove("")
                dataline = np.array(dataline[0:9], dtype=float).reshape((1, 9))
                dataMatrix[j, :] = dataline.astype('float')

# %% turn the data into meshgrids
# theta component
E_Theta = dataMatrix[:, 2] + 1j * dataMatrix[:, 3]
E_Theta = E_Theta.reshape((phiSamples, thetaSamples))
# phi component
E_Phi = dataMatrix[:, 4] + 1j * dataMatrix[:, 5]
E_Phi = E_Phi.reshape((phiSamples, thetaSamples))
# coordinates
Phi = dataMatrix[:, 1]
Phi = Phi.reshape((phiSamples, thetaSamples))
Theta = dataMatrix[:, 0]
Theta = Theta.reshape((phiSamples, thetaSamples))

D_tot = dataMatrix[:, 8]
D_tot = D_tot.reshape((phiSamples, thetaSamples))
# radiated power
G = 2 * np.pi * (np.abs(E_Theta) ** 2 + np.abs(E_Phi) ** 2) / (120 * np.pi * radiatedPower)
norm = np.max(D_tot) / np.max(G)
radiatedPower = 1. / norm
# plot
fig, ax = plt.subplots(1)
ax.contourf(Theta, Phi, 10 * np.log10(np.abs(E_Theta ** 2 + E_Phi ** 2)))
plt.show()

# %% we assume the cross polarization has no effect and define the directive gain as

G = 2 * np.pi * (np.abs(E_Theta) ** 2 + np.abs(E_Phi) ** 2) / (120 * np.pi * radiatedPower)

fig, ax = plt.subplots(1)
c = ax.pcolormesh(Theta, Phi, 10 * np.log10(G))
fig.colorbar(c, ax=ax, label='[dB]')
ax.set_title('G')
ax.set_xlabel('$\Theta$')
ax.set_ylabel('$\phi$')
plt.show()
# %%
fig, ax = plt.subplots(1)
c = ax.pcolormesh(Theta, Phi, 10 * np.log10((D_tot)))
fig.colorbar(c, ax=ax, label='[dB]')
ax.set_xlabel('$\Theta$')
ax.set_ylabel('$\phi$')
ax.set_title('Dtot')
plt.show()

# %% total power ? if the radius is 1
dtheta = (np.max(Theta) - np.min(Theta)) / thetaSamples * np.pi / 180
dphi = (np.max(Phi) - np.min(Phi)) / phiSamples * np.pi / 180
jacobian = np.sin(Theta * np.pi / 180)
Z0 = 376.73031366857
Power = np.sum(1 / (2 * Z0) * (np.abs(E_Theta) ** 2 + np.abs(E_Phi) ** 2) * jacobian) * dtheta * dphi
print(Power)
