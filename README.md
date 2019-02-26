
This directory contains basic routines to start developping 
a code for simulating the phi**4 theory with the HMC algorithm.

To compile use "make". The (very rudimentary) program is run with

phi4 infile

List of program files:

phi4.c        : main program
lattice.h     : contains global variables, the field and its dimensions
hopping.c     : initializes the hopping field
ranlxd.c      : Luescher's random number generator
evolution.c   : contains various checks the program (reversibility, hamiltonian conservation, ...)
evolutionar.c : calculates the values for magnetization (m, m^2, m^4)
hmc.c         : contain all the functions necessary to evolve the field, used in evolutionar routine
measure.c     : calcultes the observables (magnetization, susceptivity and Binder cumulant) with relative errors