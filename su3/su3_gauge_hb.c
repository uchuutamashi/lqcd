// Ref:
// G.P. Lepage, LQCD for Novices
#include "su3_utils.c"
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_fit.h>

// globals
const gsl_rng_type *T;
gsl_rng *r; // random generator
su3_matrix links[4*N*N*N*N]; // link variables
su3_matrix rands[50]; // random unitary matrices

FILE *output;

// prototypes
void setup();
void cleanup();

int
main (int argc, char** argv)
{
  setup();

  int O[4]={0,0,0,0};

  for(int n=0;n<Ntherm;n++) {
    //fprintf(output, "%f \n",plaquette(links,0,1,O));
    hb_update(links,r,rands); //thermalize
  }

  // Confinement phase should have small Polyakov VEV
  if(confined){
    while(polyakov_mean(links) > pol_thres){
      hb_update(links,r,rands); //thermalize
    }
  }

  printf("Thermalized! \n");

  for(int n=0;n<Ncf;n++){
    // Decorrelate
    for(int k=0;k<Ncor;k++){
      hb_update(links,r,rands);
    }

    if(n%5==0){printf("%d \n",n);}
    // Measurement

    //fprintf(output, "%f %f \n", wloop_mean(1,1,links),wloop_mean(1,2,links));
    //fprintf(output, "%f %f %f %f %f %f %f %f\n",twloop2_mean(1,links),twloop2_mean(2,links),twloop2_mean(3,links),twloop2_mean(4,links),twloop2_mean(5,links),twloop2_mean(6,links),twloop2_mean(7,links),twloop2_mean(8,links));
    //fprintf(output, "%f %f %f %f %f %f %f %f\n",plaq_corr(1,links),plaq_corr(2,links),plaq_corr(3,links),plaq_corr(4,links),plaq_corr(5,links),plaq_corr(6,links),plaq_corr(7,links),plaq_corr(8,links));
    //printf("%f \n",polyakov_mean(links));
    printf( "%f \n", plaq_mean_s(0,links));
    printf( "%f %f \n", wloop_mean(1,1,links),wloop_mean(1,2,links) );
  }

  printf("acc=%f\n", (double)acc/tot);

  cleanup();
  return 0;
}

void
setup (){
  acc=0;
  tot=0;
  // ready random num generator
  gsl_rng_env_setup ();
  T = gsl_rng_mt19937;
  r = gsl_rng_alloc (T);
  // fill in table of rand matrices
  for(int i=0;i<25;i++){
    rands[2*i]=su3_rand(r);
    rands[2*i+1]=su3_inv(rands[2*i]);
  }
  //init links to identities
  for(int i=0;i<4*N*N*N*N;i++){
    links[i]=su3_unit();
  }

  output = fopen("out.txt", "w+");
}

void
cleanup (){
  gsl_rng_free(r);
  fclose(output);
}
