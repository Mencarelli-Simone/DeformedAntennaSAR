# Simone Mencarelli
# September 2023
# This script is intended to produce 2 far field source files to test the rest of the code in absence of
# the actual pattern from the mechanical analysis
# This uses the Aperture object in radartools.farField to simulate first a uniform aperture and then
# an aperture with some "bending" in the y direction i.e. a shorter aperture with a circular phase shift along
# the y direction.


# %% includes
from radartools.farField import UniformAperture
import numpy as np
import matplotlib.pyplot as plt
from ffsFileWriter import ffsWrite
from dummyModalApertureField import aperture_distribution

# %% user input
# input file
filename = 'random_analysis_results/res_1mode.txt'
# output filenames
distort = 'dummyDistortedMode2.ffs'
# antenna length
L = 2  # m
# antenna width
W = 0.3  # m

# %% get field distribution

x, y, E = aperture_distribution(filename, 3e8 / 10e9)
L = np.max(y) - np.min(y)
W = np.max(x) - np.min(x)

# %% new antenna
antenna = UniformAperture(L, W)
# Create a meshgrid for the aperture with lambda/3 spacing
antenna.set_uniform_mesh(len(y), len(x))
# set the field to phase shifted  illumination (phase and amplitude)
antenna.tanEField = E.T
Y, X = np.meshgrid(antenna.l_mesh, antenna.w_mesh)
# create a theta phi meshgrid ( just the positive hemisphere is fine)
theta = np.linspace(0, np.pi / 2, 361)
phi = np.linspace(0, 2 * np.pi, 361)
T, P = np.meshgrid(theta, phi)

# %%
fig, ax = plt.subplots(1)
c = ax.pcolormesh(Y, X, np.angle(antenna.tanEField))
plt.colorbar(c, ax=ax, label='[radians]')
ax.set_xlabel('y [m]')
ax.set_ylabel('x [m]')
plt.show()

# %% compute the gain ( is going to be slow)
print('computing pattern')
E_theta, E_phi = antenna.mesh_E_field(T, P, polarization='y')  # H pol for the radar
print('plotting')
# %%
# equation 18.6.10 Orfanidis - Electromagnetic Waves and Antennas
gain = 2 * np.pi * (np.abs(E_theta) ** 2 + np.abs(E_phi) ** 2) / (antenna.eta * antenna.get_radiated_power())
fig, ax = plt.subplots(1)
ax.pcolormesh(T, P, 10 * np.log10(gain), vmin=-50, vmax=40)
plt.show()
gmax = 4 * np.pi * L * W / (antenna.c / antenna.freq) ** 2
# test
print(gmax - np.max(gain))
# pass
# %% save this to a ffs file
ffsWrite(T.reshape(-1) * 180 / np.pi, P.reshape(-1) * 180 / np.pi,
         E_theta.reshape(-1), E_phi.reshape(-1), len(phi), len(theta),
         distort,
         radiated_power=antenna.get_radiated_power(),
         accepted_power=antenna.get_radiated_power(),
         stimulated_power=antenna.get_radiated_power(),
         frequency=1e10)

# %% test, plot the main cuts

fig, ax = plt.subplots(1)

# H cut
th = np.linspace(-np.pi / 30, np.pi / 30, 501)
ph = np.array(1) * np.pi / 2
t, p = np.meshgrid(th, ph)
E_theta, E_phi = antenna.mesh_E_field(t, p, polarization='y')  # H pol for the radar
# equation 18.6.10 Orfanidis - Electromagnetic Waves and Antennas
gain = 2 * np.pi * (np.abs(E_theta) ** 2 + np.abs(E_phi) ** 2) / (antenna.eta * antenna.get_radiated_power())

ax.plot(t.reshape(-1) * 180 / np.pi, gain.reshape(-1))
ax.set_xlabel('$\Theta$ [deg]')
ax.set_ylabel('Gain')
ax.set_ylim(0, 8500)
# E cut
# th = np.linspace(-np.pi / 2, np.pi / 2, 501)
ph = np.array(1) * 0
t, p = np.meshgrid(th, ph)
E_theta, E_phi = antenna.mesh_E_field(t, p, polarization='y')  # H pol for the radar
# equation 18.6.10 Orfanidis - Electromagnetic Waves and Antennas
gain = 2 * np.pi * (np.abs(E_theta) ** 2 + np.abs(E_phi) ** 2) / (antenna.eta * antenna.get_radiated_power())
ax.plot(t.reshape(-1) * 180 / np.pi, gain.reshape(-1))
plt.show()
