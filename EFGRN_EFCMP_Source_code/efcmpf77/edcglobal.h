c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c	GLOBAL CONSTANTS
c	================
c
c	NRECMAX = max. number of observation positions
c	NZMAX = max. number of the discrete source depths
c	NRMAX = max. number of the discrete radial diatances
c	NSMAX = max. number of the source rectangles
c	NPSMAX = max. number of discrete point sources each depth
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	integer NZMAX,NRMAX,NSMAX,NPSMAX,NRECMAX,NFIELDS
	parameter(NZMAX=501,NRMAX=500)
	parameter(NSMAX=170000,NPSMAX=170000)
	parameter(NRECMAX=100000)
	parameter(NFIELDS=3)
