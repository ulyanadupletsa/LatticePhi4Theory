#ifndef HMC_H
#define HMC_H

#include "stdio.h"
#include "stdlib.h"
#include "lattice.h"
#include "ranlxd.h"
#include "math.h"
#include "phi4.h"

#ifndef HMC_C
#define HMC_C

extern void gauss();
extern double action(double kappa, double lambda);
extern double hamiltonian (double S);
extern void move_phi(double deltaT);
extern void move_mom(double deltaT, double kappa, double lambda);
extern double force(int i, double kappa, double lambda);
extern void old_phi (double oldphi[V]);
extern void reverse_mom (double revmom[V]);
extern double modulus (double oldphi[V]);
extern void reject_phi (double phi[V], double oldphi[V]);
extern void leapfrog(int step, double deltaT, double kappa, double lambda);
extern double magnetization(double phi[V]);

#endif

#endif
