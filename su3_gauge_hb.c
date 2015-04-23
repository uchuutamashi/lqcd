// Ref:
// G.P. Lepage, LQCD for Novices
#include "su3_utils.c"
#include <stdlib.h>
#include <stdio.h>

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

  /* test matrix methods
  su3_matrix A= su3_rand(r);
  su3_matrix B= su3_rand(r);
  su3_print(A);
  su3_print(B);
  su3_print(su3_mul(A,B));
  su3_print(su3_inv(A));


  return 0;
  // end test */


  int O[4]={0,0,0,0};

  for(int n=0;n<Ntherm;n++) {
    //fprintf(output, "%f \n",plaquette(links,0,1,O));
    hb_update(links,r,rands); //thermalize
  }
  printf("Thermalized! \n");

  for(int n=0;n<Ncf;n++){
    // Decorrelate
    for(int k=0;k<Ncor;k++){
      hb_update(links,r,rands);
    }

    if(n%20==0){printf("%d \n",n);}
    // Measurement
    //fprintf(output, "%f %f \n", wloop_mean(1,1,links),wloop_mean(1,2,links));
    fprintf(output, "%f %f %f %f %f \n",stat_pot_mean(1,links),stat_pot_mean(2,links),stat_pot_mean(3,links),stat_pot_mean(4,links),stat_pot_mean(5,links));
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
