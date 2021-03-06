Code vacf computes from an XDATCAR file generated by VASP the following quantities:

1) The velocity auto-correlation function via polynomial interpolation of atomic positions of order mint.

2) The vibrational density of states (VDOS) from the Fourier-transform of the velocity auto-correlation function.

3) The moments of the VDOS via direct integration and from time-derivatives of the velocity.

4) The entropy.

5) Parameters controlling the calculation are located in file 'P1' including:

ratio (of total time steps to equilibration time steps)
potim
mint
nintmax (maximum number of time steps to include in the integration over the velocity auto-correlation function)
dampfactor (controls the gaussian damping of the velocity auto-correlation function)
solid (logical.  true assumes fgas = 0.)
