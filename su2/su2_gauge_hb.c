// Ref:
// G.P. Lepage, LQCD for Novices
#include "su2_utils.c"
#include <stdlib.h>
#include <stdio.h>

// globals
const gsl_rng_type *T;
gsl_rng *r; // random generator
su2_matrix links[4*N*N*N*N]; // link variables
su2_matrix rands[50]; // random unitary matrices

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


  for(int n=0;n<Ncf;n++){
    // Decorrelate
    for(int k=0;k<Ncor;k++){
      hb_update(links,r,rands);
    }

    // Measurement
    //fprintf(output, "%f \n",plaquette(links,0,1,O));
    fprintf(output, "%f \n", wloop_mean(4,8,links));
    //fprintf(output, "%f %f %f %f %f %f %f %f \n",wloop_mean(2,N-1,links),wloop_mean(2,N,links),wloop_mean(3,N-1,links),wloop_mean(3,N,links),wloop_mean(4,N-1,links),wloop_mean(4,N,links),wloop_mean(5,N-1,links),wloop_mean(5,N,links));
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
    rands[2*i]=su2_rand(r);
    rands[2*i+1]=su2_inv(rands[2*i]);
  }
  //init links to identities
  for(int i=0;i<4*N*N*N*N;i++){
    links[i]=su2_unit();
  }

  output = fopen("out.txt", "w+");
}

void
cleanup (){
  gsl_rng_free(r);
  fclose(output);
}
