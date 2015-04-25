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
#define N_lat   2000  // Lattice size
#define a       0.1   // Lattice spacing
#define N_therm 500   // Thermalization
#define N_meas  1000  // Number of configuration to measure
#define N_skip  50    // Decorrelation
#define alpha   0     // Strength of potential (harmonic)
#define lambda  1     // Strength of potential (anharmonic)
#define epsilon 0.2   // Size of fluctuation

// something not right with the space lattice...
#define N_bin 10000
#define dx 0.0006
#define x_range 3 // x in [-range, range]

// Globals
double x[N_lat];      // Lattice
int acc,tot;          // Acceptance rate = acc/tot
double E0;            // Ground state energy
double C[N_lat];      // Correlator
double psi[N_bin];    // Wavefunction
double psip[N_bin];   // 1st derivative of wavefunction
double psipp[N_bin];  // 2nd derivative of wavefunction
double H_susy[N_bin]; // SUSY action


// Output Files
FILE *fE0, *fpsi0, *fpsi;

// Prototypes
void setup();
void cleanup();
void update();
int xtoi(double x);

int main(void){
  setup();

  // Thermalize
  for(int n=0;n<N_therm;n++){
    update();
  }

  // Calculate Wavefunction
  for(int n=0;n<N_meas;n++){
    // Decorrelate
    for(int m=0;m<N_skip;m++){ update(); }

    for(int i=0;i<N_lat;i++){
      if(abs(x[i])<x_range){
        psi[xtoi(x[i])]+=1;
      }
    }
  }

  printf("1st round acceptance rate: %f \n", (double)acc/tot);

  // Smooth out the wavefunction
  // Contributing paths are non-differentiable
  double psi_smooth[N_bin];
  for(int i=0;i<N_bin;i++){
    psi_smooth[i]=0;
    for(int k=-100;k<100;k++){
      psi_smooth[i] += psi[(i+k+N_bin)%N_bin];
    }
    psi_smooth[i]=psi_smooth[i]/200;
  }
  for(int i=0;i<N_bin;i++){
    psi[i]=psi_smooth[i];
  }




  // Calculate derivatives
  for(int i=0;i<N_bin;i++){
    // O(dx^4) approximation
    // Ref: http://www.geometrictools.com/Documentation/FiniteDifferences.pdf
    psip[i]= (-psi[(i+2)%N_bin]+8*psi[(i+1)%N_bin]-8*psi[(i-1+N_bin)%N_bin]+psi[(i-2+N_bin)%N_bin])/sq(dx);
    psipp[i]=(-psi[(i+2)%N_bin]+16*psi[(i+1)%N_bin]-30*psi[i]+16*psi[(i-1+N_bin)%N_bin]-psi[(i-2+N_bin)%N_bin])/(12*sq(dx));

    // Calculate shift in action (Ref:E Cooper et al. /Physics Reports 251 (1995) 267-385)
    // Note that dH = + sqrt(2) W'  (Old= W^2 - W'/sqrt(2), New= W^2 + W'/sqrt(2))
    if(psi[i]!=0){
      H_susy[i] = sq(psip[i]/psi[i])-psipp[i]/psi[i];
    }else{
      H_susy[i]=9999999;
    }
  }

  // Plot wavefunction
  for(int i=0;i<N_bin;i++){
    fprintf(fpsi0,"%f %f %f\n", i*dx-x_range, psip[i], psipp[i]);
  }


  // Reinitialize
  for(int i=0;i<N_bin;i++){psi[i]=0;}
  for(int i=0;i<N_lat;i++){x[i]=0;}

  // Thermalize
  for(int n=0;n<N_therm;n++){
    update();
  }

  // Calculate Wavefunction
  for(int n=0;n<N_meas;n++){
    // Decorrelate
    for(int m=0;m<N_skip;m++){ update(); }
      for(int i=0;i<N_lat;i++){
      if(abs(x[i])<x_range){
        psi[xtoi(x[i])]+=1;
      }
    }
  }


  printf("2nd round acceptance rate: %f \n", (double)acc/tot);

  // Plot wavefunction
  for(int i=0;i<N_bin;i++){
    fprintf(fpsi,"%f, %f\n", i*dx-x_range, sqrt(psi[i]/(dx*N_meas*N_lat)));
  }

  cleanup();
  return 0;
}

// Convert x into bin index
int xtoi(double x){
  return (int)((x+x_range)/dx);
}

// Calculate difference in action
double dS(int i, double x_old){
  return (x[i]*(x[i]-x[(i+1)%N_lat]-x[(i-1+N_lat)%N_lat])/a
          + alpha*a*sq(x[i]) + lambda*a*qd(x[i]) + a*H_susy[xtoi(x[i])])
          - (x_old*(x_old-x[(i+1)%N_lat]-x[(i-1+N_lat)%N_lat])/a
               + alpha*a*sq(x_old) + lambda*a*qd(x_old) + a*H_susy[xtoi(x_old)]);
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
  for(int i=0;i<N_bin;i++){ psi[i]=0; psip[i]=0; psipp[i]=0; H_susy[i]=0;}

  // Create Files
  fE0=fopen("E0.txt", "w+");
  fpsi0=fopen("psi0.txt", "w+");
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
  fclose(fpsi0);
  fclose(fpsi);
}
