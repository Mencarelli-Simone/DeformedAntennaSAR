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
# %% load patterns
ant_ref = Aperture(reference_pattern)
ant_dist = Aperture(distorted_pattern)

# %% calculation
# 1 Doppler Bandwidth (nominal 3dB beamwidth)
Bd = nominal_doppler_bandwidth(La, incidence_broadside, c / f, v_s, altitude)
print('Doppler Bandwidth: ', Bd)
# 2 Doppler Axis Incidence Axis
doppler = np.linspace(-Bd / 2, Bd / 2, 10001)
# incidence at swath edges
r0, rg0 = range_from_theta(incidence_broadside * 180 / np.pi, altitude)
rgNF = np.array((rg0 - swath / 2, rg0 + swath / 2))
rNF = range_ground_to_slant(rgNF, altitude)
rgNF, incNF = range_slant_to_ground(rNF, altitude)
# incidence axis
incidence = np.linspace(incNF[0], incNF[1], 101)  # random length
# 3 Incidence Doppler meshgrid
I, D = np.meshgrid(incidence, doppler)
# 4 Incidence Azimuth Stationary time
I, A, Tk = mesh_doppler_to_azimuth(I, D, c / f, v_s, altitude)
# 5 GCS
X, Y, Z = mesh_incidence_azimuth_to_gcs(I, A, c / f, v_s, altitude)
# 6 LCS
Xl, Yl, Zl = mesh_gcs_to_lcs(X, Y, Z, radarGeo.Bc2s, radarGeo.S_0)
# 7 LCS spherical
R, T, P = meshCart2sph(Xl, Yl, Zl)
# 8 Antenna patterns
T[np.isnan(T)] = 0
P[np.isnan(P)] = 0  # to be safe
G_ref = ant_ref.mesh_gain_pattern(T, P)
G_dist = ant_dist.mesh_gain_pattern(T, P)
# 9 SNR
