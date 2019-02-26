#ifndef EVOLUTIONAR_H
#define EVOLUTIONAR_H

#include "stdio.h"
#include "stdlib.h"
#include "lattice.h"
#include "ranlxd.h"
#include "math.h"
#include "phi4.h"
#include "hmc.h"



#ifndef EVOLUTIONAR_C
#define EVOLUTIONAR_C

extern void evolutionar(hmc_params_t *hmc_params, act_params_t *act_params);
#endif

#endif

