/* =======================================
 * Code for 8.06 Term Paper (Spring 2015)
 *
 * Author : To Chin Yu (ytc@mit.edu)
 * Purpose: 1D lattice simulation
 * V(x) = alpha*x^2 + lambda*x^4
 * =======================================
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define sq(x) pow((x),2.0)
#define cb(x) pow((x),3.0)
#define qd(x) pow((x),4.0)

// GSL Random Number Generator (https://www.gnu.org/software/gsl/)
#include <gsl/gsl_rng.h>
const gsl_rng_type * T; gsl_rng * r;

// Parameters
#define N_lat   5000  // Lattice size
#define a       0.1   // Lattice spacing
#define N_therm 500   // Thermalization
#define N_meas  1000  // Number of configuration to measure
#define N_skip  100   // Decorrelation
#define alpha   0     // Strength of potential (harmonic)
#define lambda  1     // Strength of potential (anharmonic)
#define epsilon 0.2   // Size of fluctuation

#define N_bin 1000
#define dx 0.01
#define x_range 5.0 // x in [-range, range]

// Globals
double x[N_lat];      // Lattice
int acc,tot;          // Acceptance rate = acc/tot
double E0;            // Ground state energy
double C[N_lat];      // Correlator
double psi[N_bin];    // Wavefunction

// Output Files
FILE *fE0, *fCorRatio, *fpsi;

// Prototypes
void setup();
void cleanup();
void update();

int main(void){
  setup();

  // Thermalize
  for(int n=0;n<N_therm;n++){
    update();

    // Plot E0 during thermalization
    if(0){
      E0=0;
      for(int i=0;i<N_lat;i++){
        E0 += ( 2*x[i]*lambda*cb(x[i]) + lambda*qd(x[i]) )/N_lat;
      }
      fprintf(fE0,"%f\n", E0);
    }
  }

  // Do updates and measurements
  for(int n=0;n<N_meas;n++){
    // Decorrelate
    for(int m=0;m<N_skip;m++){ update(); }
    // Measurements:
    // 1. E0

    // Plot E0
    E0=0;
    for(int i=0;i<N_lat;i++){
      E0 += ( 2*x[i]*lambda*cb(x[i]) + lambda*qd(x[i]) )/N_lat;
    }
    fprintf(fE0,"%f\n", E0);

    // 2. Correlator
    for(int i=0;i<N_lat;i++){
      for(int T=0;T<N_lat;T++){
        C[T] += x[i]*x[(i+T)%N_lat] / (N_lat*N_meas);
      }
    }

    // 3. Wavefunction
    for(int i=0;i<N_lat;i++){
      if(abs(x[i])<x_range){
        psi[(int)floor((x[i]+x_range)/dx)]+=1;
      }
    }

  }

  // Plot correlator ratio
  for(int T=0;T<N_lat;T++){
    fprintf(fCorRatio,"%f \n", C[T]/C[(T+1)%N_lat]); // for large T, this gives a(E1-E0)
  }

  // Plot wavefunction
  for(int i=0;i<N_bin;i++){
    fprintf(fpsi,"%f, %f\n", i*dx-x_range, sqrt(psi[i]/(dx*N_meas*N_lat)));
  }

  printf("Acceptance rate: %f \n", (double)acc/tot);

  cleanup();
  return 0;
}

// Calculate difference in action
double dS(int i, double x_old){
  return (x[i]*(x[i]-x[(i+1)%N_lat]-x[(i-1+N_lat)%N_lat])/a
          + alpha*a*sq(x[i]) + lambda*a*qd(x[i]))
          - (x_old*(x_old-x[(i+1)%N_lat]-x[(i-1+N_lat)%N_lat])/a
               + alpha*a*sq(x_old) + lambda*a*qd(x_old));
}

void update(){
  double x_old;
  double d; // Difference in action
  for(int i=0;i<N_lat;i++){
    x_old=x[i];
    x[i]=x[i] + epsilon * (gsl_rng_uniform(r)-0.5)*2;
    acc++; tot++;
    d=dS(i,x_old);
    if(d>0){
      if(exp(-d)<gsl_rng_uniform(r)){ x[i]=x_old; acc--; } // Revert to old value
    }
  }
}

void setup(){
  // Cold start
  for(int i=0;i<N_lat;i++){ x[i]=0; C[i]=0; }

  // Initialize variables
  acc = 0;
  tot = 0;
  E0 = 0;
  for(int i=0;i<N_bin;i++){ psi[i]=0; }

  // Create Files
  fE0=fopen("E0.txt", "w+");
  fCorRatio=fopen("CorRatio.txt", "w+");
  fpsi=fopen("psi.txt", "w+");

  // Setup random environment
  gsl_rng_env_setup();
  T = gsl_rng_mt19937;
  r = gsl_rng_alloc (T);
}

void cleanup(){
  gsl_rng_free (r);

  // Close Files
  fclose(fE0);
  fclose(fCorRatio);
  fclose(fpsi);
}
