# Simone Mencarelli
# September 2023
# Load and visualize the reference and distorted antenna pattern

# %%includes
import numpy as np
from farFieldCST import Aperture
import matplotlib.pyplot as plt
import matplotlib
from radartools.utils import *

matplotlib.use('Qt5Agg')

# %% User input
reference_pattern = 'lyceanem/NormalAntennav2.ffe'
distorted_pattern = 'lyceanem/DeformedAntennav2.ffe'
# reference_pattern = 'dummyReference.ffs'
# distorted_pattern = 'dummyDistortedMode2.ffs' # ode 2 is mode 1 actually
# %% load patterns
ant_ref = Aperture(reference_pattern)
ant_dist = Aperture(distorted_pattern)
#%%
fw = 3.45
fh = fw*.7
fontsize = 8
# %% plot gain in main cuts
theta = np.linspace((-90) * np.pi / 180, (90)* np.pi / 180, 2501)
phi_E = np.array(0)
phi_H = np.array(np.pi / 2)
# E cut
T, P = np.meshgrid(theta, phi_E)
gain_r = ant_ref.mesh_gain_pattern(T, P)
gain_d = ant_dist.mesh_gain_pattern(T, P, cubic=True) # need more samples in the pattern
fig, ax = plt.subplots(1)
ax.plot(theta * 180 / np.pi, 10*np.log10(gain_r.reshape(-1)), 'r', label='E-cut Nom.')
ax.plot(theta * 180 / np.pi, 10*np.log10(gain_d.reshape(-1)), '--r', label='E-cut Dist.')
# H cut
T, P = np.meshgrid(theta, phi_H)
gain_r = ant_ref.mesh_gain_pattern(T, P)
gain_d = ant_dist.mesh_gain_pattern(T, P)
ax.plot(theta * 180 / np.pi, 10*np.log10(gain_r.reshape(-1)), 'b', label='H-cut Nom.')
ax.plot(theta * 180 / np.pi, 10*np.log10(gain_d.reshape(-1)), '--b', label='H-cut Dist.')
ax.set_xlabel('$\Theta$ [deg]')
ax.set_ylabel('Gain ')
ax.legend(loc='lower left', prop={"size": fontsize})
fig.set_figwidth(fw)
fig.set_figheight(fh)
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label, ] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(fontsize)

plt.tight_layout()
plt.rcParams['svg.fonttype'] = 'none'
fig.savefig("niceplotsforpaper/patterncuts.svg", format="svg", dpi=600)
plt.show()
plt.show()

# %% plot gains in 2-d surfaces
# reference pattern
phi = np.linspace(0, 2 * np.pi, 571)
T, P = np.meshgrid(theta, phi)
gain_r = ant_ref.mesh_gain_pattern(T, P)
gain_d = ant_dist.mesh_gain_pattern(T, P)
# %%
fig, ax = plt.subplots(1)
c = ax.pcolormesh(T * np.cos(P), T * np.sin(P), 10*np.log10(gain_r), cmap=plt.cm.plasma, vmin=-10, vmax=20, rasterized=True)
#ax.pcolormesh(T * np.cos(P), T * np.sin(P), (gain_r))
ax.set_xlabel("$\\theta\  cos \phi$")
ax.set_ylabel("$\\theta\  sin \phi$")

cbar = fig.colorbar(c, ax=ax, label='[dB]')
fig.set_figwidth(fw)
fig.set_figheight(fh)
ax.set_aspect('equal', 'box')
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label, ] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(fontsize)
for t in (cbar.ax.get_yticklabels() + [cbar.ax.yaxis.label]):
    t.set_fontsize(fontsize)
plt.tight_layout()
plt.rcParams['svg.fonttype'] = 'none'
fig.savefig("niceplotsforpaper/ref_pattern.svg", format="svg", dpi=600)
plt.show()
plt.show()
plt.show()
#%%
fig, ax = plt.subplots(1)
c= ax.pcolormesh(T * np.cos(P), T * np.sin(P), 10*np.log10(gain_d), cmap=plt.cm.plasma, vmin=-10, vmax=20, rasterized=True)
ax.set_xlabel("$\\theta\  cos \phi$")
ax.set_ylabel("$\\theta\  sin \phi$")
cbar = fig.colorbar(c, ax=ax, label='[dB]')
fig.set_figwidth(fw)
fig.set_figheight(fh)
ax.set_aspect('equal', 'box')
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label, ] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(fontsize)
for t in (cbar.ax.get_yticklabels() + [cbar.ax.yaxis.label]):
    t.set_fontsize(fontsize)
plt.tight_layout()
plt.rcParams['svg.fonttype'] = 'none'
fig.savefig("niceplotsforpaper/dist_pattern.svg", format="svg", dpi=600)
plt.show()
plt.show()
plt.show()
plt.show()

i = np.unravel_index(np.argmax(gain_r, axis=None), gain_r.shape)

#%% 3d visual
R = 10*np.log10(gain_r)
R = np.where(R<0, 0, R)
R = np.where(R>20, 20, R)
X,Y,Z = meshSph2cart(R,T,P)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# Plot the surface.
ax.plot_surface(X, Y, Z, cmap=plt.cm.plasma)
# Tweak the limits and add latex math labels.
#ax.set_zlim(0, 100)


plt.show()