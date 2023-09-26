# Project progress:
## Abstract
compute the Azimuth impulse response (over range) for
a deformed antenna processing the signal with 
the nominal antenna pattern matched filter.
Also consider effects on SNR and resolution loss. 

## steps / scripts

- **ffsFileReader:** a script to parse CST 
far field source files into an antenna object for
SAR performance evaluation. 
- **farFieldCST:** incorporates the ffs reader into an 
"aperture" object and makes the pattern available for 
any \theta--\phi coordinates mesh by means of a 
spherical interpolator.
  - **interpolator_v2**: contains a JIT(just in time compiled)
  function for the spherical coordinates interpolation of
  patterns.
- **dummyPatterns**: generates two example patterns (saving
to ffs files.) utilizes the aperture object from radartools.farField 
to perform the pattern integration. For the distorted pattern
assumes a flat (shorter) aperture with a circular phase-taper
in the y direction, to approximate a cylindrical bend of 
a uniform radiating flat surface.
  - **radartools.farField**: contains the UniformAperture Class
  used to generate the far fields from a Huigens source aperture.
  - **ffsFileWriter**: contains a function to write a pattern 
  into a ffs file.

- **patterns_visualization**: just a script to visualize the
two patterns loaded from ffs files.

- **deformedAntennaAIR:** (Azimuth impulse response)

- **deformedAntennaSNR:** compute the core SNR i.e. the SNR for 
distributed scatterers normalized with respect to RCS, transmitted
power, noise figure, losses and chirp (noise) bandwidth. For the
nominal case and the distorted antenna case (with the signal and noise
still processed using the nominal antenna pattern)

# dependencies
- numpy 
- matplotlib
- scipy
- numba

# Getting started
rename the dummyDistorted0.fss and dummyReference0.ffs removing the "0" 
at the end and run **deformedAntennaAIR:** and **deformedAntennaSNR:**.\
Alternative ffs patterns can be visualized with **patterns_visualization** 
and input in the **user input** section of the above scripts.

# notes
the *radartools* folder is copied from the design-baseline project.
it contains some functions and objects useful to modelling the 
geometry and set the radar parameters.