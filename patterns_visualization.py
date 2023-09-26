# Simone Mencarelli
# September 2023
# Load and visualize the reference and distorted antenna pattern

# %%includes
import numpy as np
from farFieldCST import Aperture
import matplotlib.pyplot as plt
import matplotlib

#matplotlib.use('Qt5Agg')

# %% User input
reference_pattern = 'dummyReference.ffs'
distorted_pattern = 'dummyDistorted.ffs'

# %% load patterns
ant_ref = Aperture(reference_pattern)
ant_dist = Aperture(distorted_pattern)

# %% plot gain in main cuts
theta = np.linspace(-5 * np.pi / 180, 5* np.pi / 180, 2501)
phi_E = np.array(0)
phi_H = np.array(np.pi / 2)
# E cut
T, P = np.meshgrid(theta, phi_E)
gain_r = ant_ref.mesh_gain_pattern(T, P)
gain_d = ant_dist.mesh_gain_pattern(T, P, cubic=True) # need more samples in the pattern
fig, ax = plt.subplots(1)
ax.plot(theta * 180 / np.pi, (gain_r.reshape(-1)), 'r', label='E-cut Nom.')
ax.plot(theta * 180 / np.pi, (gain_d.reshape(-1)), '--r', label='E-cut Dist.')
# H cut
T, P = np.meshgrid(theta, phi_H)
gain_r = ant_ref.mesh_gain_pattern(T, P)
gain_d = ant_dist.mesh_gain_pattern(T, P)
ax.plot(theta * 180 / np.pi, (gain_r.reshape(-1)), 'b', label='H-cut Nom.')
ax.plot(theta * 180 / np.pi, (gain_d.reshape(-1)), '--b', label='H-cut Dist.')
ax.set_xlabel('$\Theta$ [deg]')
ax.set_ylabel('Gain ')
ax.legend()
plt.show()

# %% plot gains in 2-d surfaces
# reference pattern
phi = np.linspace(0, 2 * np.pi, 571)
T, P = np.meshgrid(theta, phi)
gain_r = ant_ref.mesh_gain_pattern(T, P)
gain_d = ant_dist.mesh_gain_pattern(T, P)
# %%
fig, ax = plt.subplots(1)
ax.pcolormesh(T * np.cos(P), T * np.sin(P), 10*np.log10(gain_r))
ax.set_xlabel("$\\theta\  cos \phi$")
ax.set_ylabel("$\\theta\  sin \phi$")
plt.show()
fig, ax = plt.subplots(1)
ax.pcolormesh(T * np.cos(P), T * np.sin(P), 10*np.log10(gain_d))
ax.set_xlabel("$\\theta\  cos \phi$")
ax.set_ylabel("$\\theta\  sin \phi$")
plt.show()