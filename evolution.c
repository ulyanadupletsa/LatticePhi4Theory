/* 
 *   File evolution.c
 * 
 *   void evolution (hmc_params_t *hmc_params, act_params_t *act_params)
 *      routine that takes as input parameters poiters to the infile data
 *      and then performs various checks of the program, after the 
 *      termalization of the Markov chain; in this part lambda and kappa
 *      are fixed, lambda=1.3282 and kappa=0.18169, while nstep varies;
 *      it does NOT implement the accept/reject step
 *      
 *      1) Time-reversibility
 *      2) <exp(-dH)>
 *      3) <m^2>/V^2 without the accept/reject step
 *      4) <|DH|> 
 *
 *      routines in hmc.c are called in order to solve the molecular 
 *      dynamics equations and to calculate the hamiltonian, the
 *      magnetization and to inizialize to mom[] field
 *
 */ 



#define EVOLUTION_C

#include "stdio.h"
#include "stdlib.h"
#include "lattice.h"
#include "ranlxd.h"
#include "math.h"
#include "phi4.h"
#include "hmc.h"



void evolution(hmc_params_t *hmc_params, act_params_t *act_params){

	int i,j;
	double traj_length, *deltaH, deltaT, act, Hi, Hf, *dH;
	double kappa, lambda;
    int nstep, ntherm, ntraj, nbin, nconf;
    double oldphi[V], deltaphi;
    double sum2, m, *mag2, meanM, meanH, sigmaM, sigmaH;
    double *expValue, meanExp, sigmaExp;


    FILE *DH = fopen("DH.txt", "a");
    FILE *magSquared = fopen("magSquared.txt", "a");
    FILE *expDH = fopen("expDH.txt", "a");


    traj_length = (*hmc_params).traj_length;
    nstep =  (*hmc_params).nstep;
    ntherm = (*hmc_params).ntherm;
    ntraj = (*hmc_params).ntraj;
    nbin = (*hmc_params).nbin;
    nconf = ntraj/nbin;
    deltaT = traj_length/nstep;
    kappa = (*act_params).kappa;
    lambda = (*act_params).lambda;

    mag2 = malloc(nconf*sizeof(double));
    deltaH = malloc(nconf*sizeof(double));
    expValue = malloc(nconf*sizeof(double));
    dH = malloc(nconf*sizeof(double));

    /*Termalization*/
    for (i=0; i<ntherm; i++){
        gauss();
        leapfrog(nstep, deltaT, kappa, lambda);
    }

    /* reversibility check */
    old_phi(oldphi);
    leapfrog(nstep, deltaT, kappa, lambda);
	reverse_mom(mom);
	leapfrog(nstep, deltaT, kappa, lambda);
	deltaphi = modulus(oldphi);
	printf ("Reversibility %e\n", deltaphi);


    /*Updates of field configurations and calculations
      of Hamiltonian variation and squared magnetization*/
    meanM = 0;
    meanH = 0;
    meanExp = 0;

    for(i=0; i<nconf; i++){
        deltaH[i] = 0;
        dH[i] = 0;
    }
  
    for (i=0; i<nconf; ++i){
        sum2 = 0;

        for (j=0; j<nbin; ++j){
            gauss();
            act = action(kappa, lambda);
            Hi = hamiltonian(act);
            leapfrog(nstep, deltaT, kappa, lambda);
            act = action(kappa, lambda);
            Hf = hamiltonian(act);
            m = magnetization(phi);
            sum2 += m*m;
            dH[i] += (Hf-Hi);
            deltaH[i] += fabs(Hf - Hi);
        }

        deltaH[i] /= nbin;
        dH[i] /= nbin;
        meanH += deltaH[i];
        mag2[i] = (sum2/nbin)/(V*V);
        meanM += mag2[i];
    }
    meanM = meanM/nconf;
    meanH = meanH/nconf;
    
    /*Calculations of <exp(-dH)> and of the relative standard deviation*/
    for(i=0; i<nconf; i++){
        expValue[i] = exp(-dH[i]);
        meanExp += expValue[i];
    }
    meanExp = meanExp/nconf;

    sigmaExp = 0;
    for(i=0; i<nconf; i++){
        sigmaExp += (expValue[i]-meanExp)*(expValue[i]-meanExp);
    }
    sigmaExp = sigmaExp/(nconf*(nconf-1));
    fprintf(expDH, "%lf %lf \n", meanExp, sqrt(sigmaExp));
    free(expValue);
    fclose(expDH);


    /*Calculations of the standard deviation for the squared magnetization*/
    sigmaM = 0;
    for (i=0; i<nconf; i++){
        sigmaM += (mag2[i]-meanM)*(mag2[i]-meanM);
    }
    sigmaM = sigmaM/(nconf*(nconf -1));
    fprintf(magSquared,"%d %lf %lf \n", nstep, meanM, sqrt(sigmaM));
    free(mag2);
    fclose(magSquared);


    /*Calculations of the standard deviation for the Hamiltonian variation*/
    sigmaH = 0;
    for (i=0; i<nconf; i++){
        sigmaH += (deltaH[i]-meanH)*(deltaH[i]-meanH);
    }
    sigmaH = sigmaH/(nconf*(nconf -1));
    fprintf(DH,"%d %lf %lf \n", nstep, meanH, sqrt(sigmaH));
    free(deltaH);
    fclose(DH);
}
