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

# %% user input

# antenna length
L = 2  # m
# antenna width
W = 0.3  # m
# the frequency is 10 GHz by default                          z
# bending max displacement in the z direction   ____..--''     ^
#                                               __________     .  > y
bend = 5e-2  # m ( 5 cm max)

# %% Flat antenna setup
antenna = UniformAperture(L, W)
# Create a meshgrid for the aperture with lambda/3 spacing
antenna.set_uniform_mesh_resolution(1 / 3 * antenna.c / 10e9, 1 / 3 * antenna.c / 10e9)
# set the field to uniform illumination (phase and amplitude)
antenna.tanEField = np.ones_like(antenna.tanEField)
# create a theta phi meshgrid ( just the positive emisphere is fine)
theta = np.linspace(0, np.pi / 2, 481)
phi = np.linspace(0, 2 * np.pi, 803)
T, P = np.meshgrid(theta, phi)
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
         'dummyReference.ffs',
         radiated_power=antenna.get_radiated_power(),
         accepted_power=antenna.get_radiated_power(),
         stimulated_power=antenna.get_radiated_power(),
         frequency=1e10)

# %% Now lets make a distorted antenna, assuming a circular curvature
from scipy.optimize import fsolve
import scipy as sp

# bending radius
fun = lambda r: r * np.arccos(1 - bend / r) - L
r = float(fsolve(fun, 10))
alpha = L / float(r)
# length shrinking
L1 = np.sin(alpha) * r
# wavelength dispacement
phase = lambda y: ((r - np.sqrt(r ** 2 - (L1 / 2 - y) ** 2)) / (3e8 / 1e10)) % (2 * np.pi)

# new antenna
antenna = UniformAperture(L1, W)
# Create a meshgrid for the aperture with lambda/3 spacing
antenna.set_uniform_mesh_resolution(1 / 3 * 3e8 / 10e9, 1 / 3 * 3e8 / 10e9)
# set the field to phase shifted  illumination (phase and amplitude)
Y, X = np.meshgrid(antenna.l_mesh, antenna.w_mesh)
antenna.tanEField = np.ones_like(antenna.tanEField) * np.exp(1j * phase(Y))

fig, ax = plt.subplots(1)
ax.pcolormesh(Y, X, np.angle(antenna.tanEField))
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
         'dummyDistorted.ffs',
         radiated_power=antenna.get_radiated_power(),
         accepted_power=antenna.get_radiated_power(),
         stimulated_power=antenna.get_radiated_power(),
         frequency=1e10)

# %% test, plot the main cuts

fig, ax = plt.subplots(1)

# H cut
th = np.linspace(-np.pi / 30, np.pi / 30, 501)
ph = np.array(1) * np.pi/2
t, p = np.meshgrid(th, ph)
E_theta, E_phi = antenna.mesh_E_field(t, p, polarization='y')  # H pol for the radar
# equation 18.6.10 Orfanidis - Electromagnetic Waves and Antennas
gain = 2 * np.pi * (np.abs(E_theta) ** 2 + np.abs(E_phi) ** 2) / (antenna.eta * antenna.get_radiated_power())

ax.plot(t.reshape(-1)*180/np.pi, gain.reshape(-1))
ax.set_xlabel('$\Theta$ [deg]')
ax.set_ylabel('Gain [dB]')
# E cut
#th = np.linspace(-np.pi / 2, np.pi / 2, 501)
ph = np.array(1) * 0
t, p = np.meshgrid(th, ph)
E_theta, E_phi = antenna.mesh_E_field(t, p, polarization='y')  # H pol for the radar
# equation 18.6.10 Orfanidis - Electromagnetic Waves and Antennas
gain = 2 * np.pi * (np.abs(E_theta) ** 2 + np.abs(E_phi) ** 2) / (antenna.eta * antenna.get_radiated_power())
ax.plot(t.reshape(-1)*180/np.pi, gain.reshape(-1))
plt.show()
