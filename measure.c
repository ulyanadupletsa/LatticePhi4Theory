/* 
 *   File measure.c
 * 
 *   This module calculates the observables that are used for
 *   spontaneous symmetry breaking study. It reads through the 
 *   data files, magnet1.txt, magnet2.txt and magnet4.txt, 
 *   calculated in evolutionar.c and calculates
 *      1) <m^2>/V^2
 *      2) chi = [<m^2>-<|m|>^2]/V -> susceptivity
 *      3) binder = [<m^4>]/[<m^2>^2] -> Binder cumulant
 *
 *    Jackknife error method is implemented for 2) and 3) 
 *    and usual standard deviation is calculated for 1)
 * 
 *   It produces 3 output files, containing k value,
 *   observables value and associated error
 */ 



#define MEASURE_C

#include "stdio.h"
#include "stdlib.h"
#include "lattice.h"
#include "ranlxd.h"
#include "math.h"
#include "phi4.h"
#include "hmc.h"

void measure(int nconf, int nstep, double kappa){

	FILE *magnet1 = fopen("magnet1.txt", "r");
    FILE *magnet2 = fopen("magnet2.txt", "r");
    FILE *magnet4 = fopen("magnet4.txt", "r");
    FILE *squaredMagnetization = fopen("squaredMagnetization.txt", "a");
    FILE *susceptivity = fopen("susceptivity.txt", "a");
    FILE *binderCumulant = fopen("binderCumulant.txt", "a");

    int i;
    double *mag, *mag2, *mag4;
    double *chi, *binder, *clusterChi, *clusterBinder;
    double sum, sum2, sum4;
    double mean, mean2, mean4;
    double meansqrtmag, meanChi, meanBinder;
    double sigmaMag2, sigmaChi, sigmaBinder;
    double volume;

    mag = malloc(nconf*sizeof(double));
    mag2 = malloc(nconf*sizeof(double));
    mag4 = malloc(nconf*sizeof(double));

 
    chi = malloc(nconf*sizeof(double));
    binder = malloc(nconf*sizeof(double));

    clusterChi = malloc(nconf*sizeof(double));
    clusterBinder = malloc(nconf*sizeof(double));

    volume = (double) L*L*L;

    sum = 0;
    sum2 =0;
    sum4 = 0;

    for(i=0; i<nconf; i++){
    	fscanf(magnet1, "%lf", &mag[i]);
    	sum += mag[i];
    }
    fclose(magnet1);

    for(i=0; i<nconf; i++){
    	fscanf(magnet2, "%lf", &mag2[i]);
    	sum2 += mag2[i];
    }
    fclose(magnet2);

    for(i=0; i<nconf; i++){
    	fscanf(magnet4, "%lf", &mag4[i]);
    	sum4 += mag4[i];
    }
    fclose(magnet4);

    mean = sum/nconf;
    mean2 = sum2/nconf;
    mean4 = sum4/nconf;


    /*Calculation of susceptivity and Binder cumulant*/
    for(i=0; i<nconf; i++){
    	chi[i] = (mag2[i] - mag[i]*mag[i])/volume;
    	binder[i] = (mag4[i])/((mag2[i])*(mag2[i]));
    }
    free(mag);
    free(mag4);

    
    meanChi = (mean2-mean*mean)/volume;
    meanBinder = (mean4)/(mean2*mean2);
    
    for(i=0; i<nconf; i++){
    	clusterChi[i] = meanChi - (chi[i]-meanChi)/(nconf-1);
    	clusterBinder[i] = meanBinder - (binder[i]-meanBinder)/(nconf-1); 
    }
    free(chi);
    free(binder);

  
    sigmaChi = 0;
    sigmaBinder = 0;

    for(i=0; i<nconf; i++){
    	sigmaChi += (clusterChi[i]-meanChi)*(clusterChi[i]-meanChi);
    	sigmaBinder += (clusterBinder[i]-meanBinder)*(clusterBinder[i]-meanBinder);
    }

    free(clusterChi);
    free(clusterBinder);

    sigmaChi = sigmaChi*(nconf-1)/nconf;
    sigmaBinder = sigmaBinder*(nconf-1)/nconf;


    meansqrtmag = mean2/(volume*volume);
    for(i=0; i<nconf; i++){
        mag2[i] = mag2[i]/(volume*volume);
    }
    sigmaMag2 = 0;
    for(i=0; i<nconf; i++){
        sigmaMag2 += (mag2[i]-meansqrtmag)*(mag2[i]-meansqrtmag);
    }
    free(mag2);
    sigmaMag2 = sigmaMag2/(nconf*(nconf-1));
   

    fprintf(squaredMagnetization, "%d %lf %lf\n", nstep, meansqrtmag, sqrt(sigmaMag2));
    /*fprintf(squaredMagnetization, "%lf %lf %lf \n",kappa, meansqrtmag, sqrt(sigmaMag2));*/
    fprintf(susceptivity, "%lf %lf %lf \n", kappa, meanChi, sqrt(sigmaChi));
    fprintf(binderCumulant, "%lf %lf %f \n", kappa, meanBinder, sqrt(sigmaBinder));
    
    fclose(squaredMagnetization);
    fclose(susceptivity);
    fclose(binderCumulant);
}


