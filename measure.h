#ifndef MEASURE_H
#define MEASURE_H

#include "stdio.h"
#include "stdlib.h"
#include "lattice.h"
#include "ranlxd.h"
#include "math.h"
#include "phi4.h"
#include "hmc.h"
#include "evolutionar.h"



#ifndef MEASURE_C
#define MEASURE_C

extern void measure(int nconf, int nstep, double kappa);

#endif

#endif

