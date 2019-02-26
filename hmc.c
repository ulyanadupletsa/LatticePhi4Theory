/*
 *  File hmc.c
 *
 *   This module contains different functions needed to solve the molecular
 *   dynamics equations
 *----------------------------------------------------------------------------
 *   double action(void)
 *      This routine computes the action S[phi] for the global field phi in
 *      lattice.h with the parameters kappa and lambda from act_params.
 *      S = Sum_x [ -2*kappa*sum_mu phi_x phi_{x+mu}+phi_x^2+lambda(phi_x^2-1)^2]
 *---------------------------------------------------------------------------- 
 *   double hamiltonian (double S)
 *      calculates the hamiltonian of the system given the action S calculated
 *      in the preceding step
 *---------------------------------------------------------------------------- 
 *   void gauss()
 *      implements Box-Muller method to generate random momentum field according to
 *      a gaussian distribution
 *---------------------------------------------------------------------------- 
 *   void move_phi(double deltaT)
 *      function of the leapfrog integrator that updates the phi field
 *----------------------------------------------------------------------------
 *   void move_mom(double deltaT, double kappa, double lambda)
 *      function of the leapfrog integrator that updates the mom field; it
 *      calls
 *          double force(int i, double kappa, double lambda) 
 *              function that calculates the force field associated with the
 *              hamiltonian
 *---------------------------------------------------------------------------- 
 *   void leapfrog(int nstep, double deltaT, double kappa, double lambda)
 *      uses the results of move_mom and move_phi to implement the leapfrog
 *      integrator
 *----------------------------------------------------------------------------
 *   void old_phi (double oldphi[V])
 *      function to store the prevoius values of phi field
 *---------------------------------------------------------------------------- 
 *   void reverse_mom (double revmom[V])
 *      function that changes sign to the mom field configurations, used
 *      in performing reversibility check
 *---------------------------------------------------------------------------- 
 *   double modulus (double oldphi[V])
 *      calculates the difference in absolute value between two successive 
 *      configurations of phi field; it is used in time-reversibility check
 *----------------------------------------------------------------------------  
 *   extern double magnetization(double phi[V])
 *      calculates the magnetization at a given field configuration
 *
*/

#define HMC_C

#include "stdio.h"
#include "stdlib.h"
#include "lattice.h"
#include "ranlxd.h"
#include "math.h"
#include "phi4.h"
#include "hmc.h"


void gauss(){
    int i;
    double x[2];
    double tpi;

    tpi = 8*atan(1);

    i = 0;
    while (i<(V-1)){
        ranlxd(x, 2);
        mom[i] = sqrt(-2*log(1-x[0]))*cos(tpi*(1-x[1]));
        mom[i+1] = sqrt(-2*log(1-x[0]))*sin(tpi*(1-x[1]));
        i += 2;
    }

}

double action(double kappa, double lambda)
{
    int i,mu;
    double phin,S,phi2;

    S=0;

    /* loop over all sites */
    for (i=0;i<V;i++)
    {
	/*sum over neighbors in positive direction*/
	phin=0;
	for (mu=0;mu<D;mu++){
        phin+=phi[hop[i][mu]];
    } 

	phi2=phi[i]*phi[i];
	S+=-2*kappa*phin*phi[i]+phi2+lambda*(phi2-1.0)*(phi2-1.0);
    }

    return S;
}


double hamiltonian (double S){

    int i;
    double sumMom2;
    double H;

    sumMom2 = 0;

    for (i=0; i<V; i++){
        sumMom2 += mom[i]*mom[i];
    }
    
    H = 0.5*sumMom2 + S;
    return H;
}

double force(int i, double kappa, double lambda){

	double phin;
    double force;
    int mu;

 

    phin = 0;
    for (mu=0;mu<D;mu++){
    	phin += phi[hop[i][mu]] + phi[hop[i][D+mu]];
    } 

    force = -2*kappa*phin + 2*phi[i] + 4*lambda*(phi[i]*phi[i]-1.0)*phi[i];

    return force;
}

void move_phi(double deltaT){

	int j;

    for(j=0; j<V; j++){
        phi[j] = phi[j] + deltaT*mom[j];
    }
}



void move_mom(double deltaT, double kappa, double lambda){

	int j;
	double f;
    

    for (j=0; j<V; j++){
    	f = force(j, kappa, lambda);
        mom[j] = mom[j] - 0.5*deltaT*f;
    }

}

void old_phi (double oldphi[V]){
	int i;
	for (i=0; i<V; ++i){
		oldphi[i] = phi[i];
	}
}

void reverse_mom (double revmom[V]){
	int i;
	for(i=0; i<V; ++i){
		revmom[i] = (-1)*mom[i];
	}
}

double modulus (double oldphi[V]){
	int i;
	double sum = 0;

	for(i=0; i<V; ++i){
		sum += fabs(phi[i]-oldphi[i]);
	}

	sum /= V;

	return sum;
}


void leapfrog(int nstep, double deltaT, double kappa, double lambda){
    int k;
    for (k=0; k<nstep; ++k){
                move_mom(deltaT, kappa, lambda);
                move_phi(deltaT);
                move_mom(deltaT, kappa, lambda);
            }
}

extern double magnetization(double phi[V]){
    int k;
    double magnetization;

    magnetization = 0;
    for (k = 0; k < V; k++) {
        magnetization += phi[k];
    }   
    magnetization = fabs(magnetization);

    return magnetization;
}


