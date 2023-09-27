# Simone Mencarelli
# September 2023
# Azimuth impulse response for deformed antenna
import matplotlib.pyplot as plt

# %% includes
from radartools.spherical_earth_geometry_radar import *
from radartools.design_functions import *
from farFieldCST import Aperture
import numpy as np
from numpy.fft import fft, ifft, fftshift, ifftshift

# %% User input
# broadside incidence angle
incidence_broadside = 25 * np.pi / 180
altitude = 500e3
# frequency
f = 10e9
# antenna length (or equivalent for 3db beamwidth)
La = 2
# speed of light
c = 299792458.0
# range swath
swath = 100e3
# pattern files
reference_pattern = 'dummyReference.ffs'
distorted_pattern = 'dummyDistorted.ffs'

# %% set scene geometry
radarGeo = RadarGeometry()
# set rotation
looking_angle = incidence_angle_to_looking_angle(incidence_broadside, altitude)
radarGeo.set_rotation(looking_angle, 0, 0)
# set position
radarGeo.set_initial_position(0, 0, altitude)
# set orbital speed
v_s = radarGeo.orbital_speed()
radarGeo.set_speed(v_s)
# that's it
# %% load antenna patterns
ant_ref = Aperture(reference_pattern)
ant_dist = Aperture(distorted_pattern)

# %% calculation
## 1 Doppler Bandwidth (nominal 3dB beamwidth)
Bd = nominal_doppler_bandwidth(La, incidence_broadside, c / f, v_s, altitude)
print('Doppler Bandwidth: ', Bd)
## 2 Doppler Axis Incidence Axis
doppler = np.linspace(-Bd / 2, Bd / 2, 10001)
# incidence at swath edges
r0, rg0 = range_from_theta(incidence_broadside * 180 / np.pi, altitude)
rgNF = np.array((rg0 - swath / 2, rg0 + swath / 2))
rNF = range_ground_to_slant(rgNF, altitude)
rgNF, incNF = range_slant_to_ground(rNF, altitude)
# incidence axis
incidence = np.linspace(incNF[0], incNF[1], 101)  # random length
## 3 Incidence Doppler meshgrid
I, D = np.meshgrid(incidence, doppler)
## 4 Incidence Azimuth Stationary time
I, A, Tk = mesh_doppler_to_azimuth(I, D, c / f, v_s, altitude)
## 5 GCS
X, Y, Z = mesh_incidence_azimuth_to_gcs(I, A, c / f, v_s, altitude)
## 6 LCS
Xl, Yl, Zl = mesh_gcs_to_lcs(X, Y, Z, radarGeo.Bc2s, radarGeo.S_0)
## 7 LCS spherical
R, T, P = meshCart2sph(Xl, Yl, Zl)
## 8 Antenna patterns
T[np.isnan(T)] = 0
P[np.isnan(P)] = 0  # to be safe
G_ref = ant_ref.mesh_gain_pattern(T, P)
G_dist = ant_dist.mesh_gain_pattern(T, P)
## 9 Integrand
# the matched filter amplitude
H = 1 / (stationary_phase_amplitude_multiplier(I, Tk, c / f, v_s, altitude) * G_ref)
denom = integrate.simps(H ** 2, D, axis=0)
# antenna weights
Wa = stationary_phase_amplitude_multiplier(I, Tk, c / f, v_s, altitude) * G_dist
numer = integrate.simps(H * Wa, D, axis=0) ** 2
## 10 Core - SNR
# Boltzman constant
k_boltz = 1.380649E-23  # J/K
# the sin of the incidence angle at each ground range point
sin_theta_i = sin(incidence)
# earth radius
re = 6371e3
# The range at each ground range point
r0 = re * (np.sqrt(cos(incidence) ** 2 + 2 * altitude / re + altitude ** 2 / re ** 2) - cos(incidence))
max_gain = ant_ref.max_gain()
# print(max_gain)
R0_mesh = re * (np.sqrt(cos(I) ** 2 + 2 * altitude / re + altitude ** 2 / re ** 2) - cos(I))
# cosine of earth-centric elevation angle of azimuth circle parallel to the orbital plane on the sphere
cos_theta_e = (re + R0_mesh * cos(I)) / (re + altitude)
vg = v_s / (re + altitude) * re * cos_theta_e[0, :]
# the equation for the reference case snr is
SNR_core_ref = (c / f) ** 2 * max_gain ** 2 * c * Bd * vg / (
        128 * np.pi ** 3 * r0 ** 4 * k_boltz * sin_theta_i * denom)
# the equation for the distorted antenna (general) case for the snr is
SNR_core = (c / f) ** 2 * max_gain ** 2 * c * vg * numer / (
        128 * np.pi ** 3 * r0 ** 4 * k_boltz * sin_theta_i * Bd * denom)
# azimuth resolutions
daz = vg / Bd

# %% Plotting

fig, ax = plt.subplots(1)
ax.plot(incidence * 180 / np.pi, 10 * np.log10(SNR_core_ref), '--', label='reference')
ax.plot(incidence * 180 / np.pi, 10 * np.log10(SNR_core), '--', label='distorted antenna')
ax.legend()
ax.set_xlabel('incidence angle [deg]')
ax.set_ylabel('core SNR [dB]')
plt.show()

fig, ax = plt.subplots(1)
ax.plot(incidence * 180 / np.pi, (SNR_core_ref), '--', label='reference')
ax.plot(incidence * 180 / np.pi, (SNR_core), '--', label='distorted antenna')
ax.legend()
ax.set_xlabel('incidence angle [deg]')
ax.set_ylabel('core SNR [neper]')
plt.show()

rs, rg = range_from_theta(incidence * 180 / np.pi, altitude)

fig, ax = plt.subplots(1)
ax.plot(rg, 10 * np.log10(SNR_core_ref), '--', label='reference')
ax.plot(rg, 10 * np.log10(SNR_core), '--', label='distorted antenna')
ax.legend()
ax.set_xlabel('ground range [m]')
ax.set_ylabel('core SNR [dB]')
plt.show()

fig, ax = plt.subplots(1)
ax.plot(rg, (SNR_core_ref), '--', label='reference')
ax.plot(rg, (SNR_core), '--', label='distorted antenna')
ax.legend()
ax.set_xlabel('ground range [m]')
ax.set_ylabel('core SNR [neper]')
plt.show()
