/* 
 *   File phi4.c
 *
 *   Contains the main program and a few other routines from which the
 *   code to simulate the phi**4 theory can be built. Routines for reading 
 *   in the main parameters of the action are provided.
 *
 * 
 *   static int get_val(FILE* fp, char *str, char* fmt,  void* val)
 *      Routine which reads one line from the input file.
 *      Format of the lines is <keyword> <value>.
 *      Checks if the keyword in string str matches,
 *      then gets the value according to the format in fmt
 *
 *     
 *   static int read_input(char *input)
 *      Parses the input file (format as specified in get_val)
 *      and prints the parameters onto the screen. Currently
 *      it reads the basic values for the action and also for the 
 *      future HMC and the seed of the random number generator.
 *
 */ 

#define CONTROL
#include "phi4.h"
#include "string.h"
#include "hmc.h"
#include "evolution.h"
#include "evolutionar.h"
#include "measure.h"
#include "lattice.h"


/*  
 *  data structures to store all the parameters of the algorithm,
 *  and action defined in phi4.h
 *  seed for initializing the random numbers
 */

static hmc_params_t hmc_params;
static act_params_t act_params;
static int seed;


static int get_val(FILE* fp, char *str, char* fmt,  void* val)
{
    char c[128];

    if(1!=fscanf(fp,"%s",c))
    {
	fprintf(stderr,"Error reading input file at %s\n",str);
	exit(1);
    }

    if(strcmp(str,c)!=0)
    {
	fprintf(stderr,"Error reading input file expected %s found %s\n",str,c);
	exit(1);
    }

    if(1!=fscanf(fp,fmt,val))
    {
	fprintf(stderr,"Error reading input file at %s\n",str);
	fprintf(stderr,"Cannot read value format %s\n",fmt);
	exit(1);
    }

    return 0;

}


static int read_input(char *input)
{
    FILE* fp;

    fp=fopen(input,"r");
    if (fp==NULL) {
	fprintf(stderr, "Cannot open input file %s \n",input);
	exit(1);
    }

    get_val(fp, "kappa",       "%lf",&act_params.kappa  );
    get_val(fp, "lambda",      "%lf",&act_params.lambda );
    get_val(fp, "ntherm",      "%i" ,&hmc_params.ntherm );
    get_val(fp, "ntraj",       "%i" ,&hmc_params.ntraj  );
    get_val(fp, "traj_length", "%lf",&hmc_params.traj_length);
    get_val(fp, "nstep",       "%i" ,&hmc_params.nstep  );
    get_val(fp, "seed",        "%i" ,&seed   );
    get_val(fp, "nbin",        "%i" ,&hmc_params.nbin  );



   
    printf("PARAMETERS\n");
    printf("L              %i\n", L);
    printf("DIM            %i\n", D);
    printf("kappa          %f\n", act_params.kappa);
    printf("lambda         %f\n", act_params.lambda);
    printf("ntherm         %i\n", hmc_params.ntherm);
    printf("ntraj          %i\n", hmc_params.ntraj);
    printf("traj_length    %f\n", hmc_params.traj_length);
    printf("nstep          %i\n", hmc_params.nstep);
    printf("nbin           %i\n", hmc_params.nbin);
    printf("END PARAMETERS\n");

    return 0;
}


int main(int argc, char* argv[]){

    hmc_params_t *hmc = &hmc_params;
    act_params_t *act = &act_params; 

    if (argc != 2) {
	fprintf(stderr, "Number of arguments not correct\n");
	fprintf(stderr, "Usage: %s <infile> \n",argv[0]);
	exit(1);
    }

    /* get the parameters from the input file */
    read_input(argv[1]);
    /* initialize the random number generator */
    rlxd_init(1,seed);
    /* initialize the nearest neighbour field */
    hopping(hop);
    /* initialize the phi field */
    ranlxd(phi,V);

/*  After initiliazing the field, one can choose between executing
 *  various checks without implementing the accept/reject step, 
 *  or run the evolutionar() routine which calculates all the 
 *  the observables need to study spaontaneous symmetry breaking 
 *  of the system
 */
    /*evolution(hmc, act);*/
    evolutionar(hmc, act);
  

    return 0;
}

