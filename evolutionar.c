/* 
 *   File evolutionar.c
 * 
 *   void evolutionar (hmc_params_t *hmc_params, act_params_t *act_params)
 *      routine that takes as input parameters poiters to the infile data
 *      and then calculates powers of magnetization needed in the last
 *      part; after the Markov chain is thermalized; in this part lambda
 *      is fixed, lambda=1.145, and kappa varies in [0.15,0.23]; here is
 *      implemented the ACCEPT/REJECT STEP
 *
 *      It produces 3 output files: magnet1.txt, magnet2.txt and magnet4.txt
 *
 *      routines in hmc.c are called in order to solve the molecular 
 *      dynamics equations and to calculate the hamiltonian, the
 *      magnetization and to inizialize to mom[] field
 *
 *      To acceptance rate is held under control, printing on terminal
 %      the percentage of acceptance after each run
 *
 */ 



#define EVOLUTIONAR_C

#include "stdio.h"
#include "stdlib.h"
#include "lattice.h"
#include "ranlxd.h"
#include "math.h"
#include "phi4.h"
#include "hmc.h"
#include "measure.h"


void evolutionar(hmc_params_t *hmc_params, act_params_t *act_params){

    FILE *magnet1 = fopen("magnet1.txt", "w");
    FILE *magnet2 = fopen("magnet2.txt", "w");
    FILE *magnet4 = fopen("magnet4.txt", "w");

	int i, j, k, nreject;
	double traj_length, deltaT, act, Hi, Hf, volume;
	double kappa, lambda;
    int nstep, ntherm, ntraj, nbin, nconf;
    double oldphi[V], r[2], dH;
    double sum1, sum2, sum4, m, acceptanceRate;


    traj_length = (*hmc_params).traj_length;
    nstep =  (*hmc_params).nstep;
    ntherm = (*hmc_params).ntherm;
    ntraj = (*hmc_params).ntraj;
    nbin = (*hmc_params).nbin;
    nconf = ntraj/nbin;
    deltaT = traj_length/nstep;
    kappa = (*act_params).kappa;
    lambda = (*act_params).lambda;

    /*Termalization*/
    for (i=0; i<ntherm; i++){
        gauss();
        old_phi(oldphi);
        act = action(kappa, lambda);
        Hi = hamiltonian(act);
        leapfrog(nstep, deltaT, kappa, lambda);
        act = action(kappa, lambda);
        Hf = hamiltonian(act);
        dH = Hf - Hi;
        if(exp(-(Hf - Hi))<1){
            ranlxd(r,1);
            if(exp(-(Hf - Hi))<r[0]){
                for(j=0; j<V; ++j){
                    phi[j] = oldphi[j];
                }
            }
        }
    }


    
    volume = L*L*L;
    nreject = 0;
    for (i=0; i<nconf; ++i){
        sum1 = 0;
        sum2 = 0;
        sum4 = 0;

       
        for (j=0; j<nbin; ++j){
            gauss();
            old_phi(oldphi);
            act = action(kappa, lambda);
            Hi = hamiltonian(act);
            leapfrog(nstep, deltaT, kappa, lambda);
            act = action(kappa, lambda);
            Hf = hamiltonian(act);
            dH = Hf - Hi;
            if (exp(-dH)>0){
                ranlxd(r,1);
                if(exp(-dH)<r[0]){
                    for(k=0; k<V; ++k){
                        phi[k] = oldphi[k];
                    }
                    nreject += 1;
                }
            }
            m = magnetization(phi);
            sum1 = sum1 + m;
            sum2 = sum2 + m*m;
            sum4 = sum4 + m*m*m*m;
        }

        sum1 = (sum1/nbin);
        sum2 = (sum2/nbin);
        sum4 = (sum4/nbin);
    
        fprintf(magnet1, "%lf\n", fabs(sum1));
        fprintf(magnet2, "%lf\n", sum2);
        fprintf(magnet4, "%lf\n", sum4);
    }

    acceptanceRate = (double) (ntraj-nreject)/ntraj;
    printf("Acceptance rate = %f\n", acceptanceRate) ;

    fclose(magnet1);
    fclose(magnet2);
    fclose(magnet4);

    measure(nconf, nstep, kappa);
}

