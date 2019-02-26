#ifndef EVOLUTION_H
#define EVOLUTION_H

#include "stdio.h"
#include "stdlib.h"
#include "lattice.h"
#include "ranlxd.h"
#include "math.h"
#include "phi4.h"
#include "hmc.h"



#ifndef EVOLUTION_C
#define EVOLUTION_C

extern void evolution(hmc_params_t *hmc_params, act_params_t *act_params);
#endif

#endif

