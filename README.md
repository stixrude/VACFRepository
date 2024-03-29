VACF by L. Stixrude 2019-

Citations.  If you use VACF, please cite the following paper:

Wilson, A. and L. Stixrude, Entropy, dynamics, and freezing of CaSiO3 liquid, Geochimica et Cosmochimica Acta, 302, 1-17, 2021.

Code vacf computes from an XDATCAR file generated by VASP the following quantities:

1) The velocity auto-correlation function via polynomial interpolation of atomic positions of order mint.

2) The vibrational density of states (VDOS) from the Fourier-transform of the velocity auto-correlation function.

3) The moments of the VDOS via direct integration and from time-derivatives of the velocity.

4) The entropy.

Parameters controlling the calculation are located in file 'P1' including:

ratio (of total time steps to equilibration time steps)
potim
mint
nintmax (maximum number of time steps to include in the integration over the velocity auto-correlation function)
dampfactor (controls the gaussian damping of the velocity auto-correlation function)
solid (logical.  true assumes fgas = 0.)
