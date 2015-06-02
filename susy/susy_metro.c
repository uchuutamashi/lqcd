/* =======================================
 * SUSY
 *
 * TODO:
 * - still giving qualitatively wrong answer
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

// GSL Spline (for smoothing)
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
gsl_interp_accel *spline_accel;
gsl_spline *spline;
#define N_spl 5       // Factor of N_bin +1 (end points inclusive)

// Parameters
#define N_lat   2000  // Lattice size
#define a       0.1   // Lattice spacing
#define N_therm 500   // Thermalization
#define N_meas  1000  // Number of configuration to measure
#define N_skip  50    // Decorrelation
#define alpha   0     // Strength of potential (harmonic)
#define lambda  1     // Strength of potential (anharmonic)
#define epsilon 0.2   // Size of fluctuation

// SUSY Parameters
#define epsilon_susy 0.4   // Size of fluctuation


#define N_bin 100
#define dx 0.04   // x_range*2 / bin
#define x_range 2 // x in [-range, range]

// Globals
double x[N_lat];           // Lattice
int N_acc,N_tot;           // Acceptance rate = acc/tot
double E0;                 // Ground state energy
double C[N_lat];           // Correlator
double psi[N_bin];         // Wavefunction
double psip[N_bin];        // 1st derivative of wavefunction
double psipp[N_bin];       // 2nd derivative of wavefunction
double H_susy[N_bin];      // SUSY action


// Output Files
FILE *fE0, *fpsi0, *fpsi;

// Prototypes
void setup();
void cleanup();
void update();
void update_susy();
int xtoi(double x);
double itox(int i);

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

  printf("1st round acceptance rate: %f \n", (double)N_acc/N_tot);

  // Normalize wavefunction
  for(int i=0;i<N_bin;i++){
    psi[i] = sqrt(psi[i]/(dx*N_meas*N_lat));
  }

  // Pick points for smoothing spline
  double x_arr[N_spl];
  double y_arr[N_spl];
  int skip = nearbyint(N_bin/(N_spl-1));
  for(int i=0;i<N_bin;i++){
    if(i%skip==0){
      x_arr[(i/skip)]=itox(i);
      y_arr[(i/skip)]=psi[i];
    }
  }
  // manual add end point
  x_arr[N_spl-1]=x_range;
  y_arr[N_spl-1]=psi[N_bin-1];

  // Smoothing spline
  gsl_spline_init(spline, x_arr, y_arr, N_spl);

  for(int i=0;i<N_bin;i++){
    psi[i]  = gsl_spline_eval(spline,itox(i),spline_accel);
    if(psi[i]<0){psi[i]=0;}
    psip[i] = gsl_spline_eval_deriv(spline, itox(i), spline_accel);
    psipp[i] = gsl_spline_eval_deriv2(spline, itox(i), spline_accel);
  }

  // Smoothing second derivative
  for(int i=0;i<N_bin;i++){
    if(i%skip==0){
      y_arr[(i/skip)]=psipp[i];
    }
  }
  // manual add end point
  y_arr[N_spl-1]=psipp[N_bin-1];
  // Smoothing spline
  gsl_spline_init(spline, x_arr, y_arr, N_spl);

  for(int i=0;i<N_bin;i++){
    psipp[i]  = gsl_spline_eval(spline,itox(i),spline_accel);

    if(psi[i]>0){
      H_susy[i] =  - (psipp[i]/psi[i]-sq(psip[i]/psi[i]))/sqrt(2); // TODO: This is technically wrong...
    }else{
      H_susy[i]=9999999;
    }
  }

  // Plot wavefunction and derivatives
  for(int i=0;i<N_bin;i++){
    fprintf(fpsi0,"%f %f %f %f\n", itox(i), psi[i], psip[i], psipp[i]);
    fprintf(fE0,"%f %f\n", itox(i), H_susy[i]);
  }


  // Reinitialize
  for(int i=0;i<N_bin;i++){psi[i]=0;}
  for(int i=0;i<N_lat;i++){x[i]=0;}
  N_acc=0;N_tot=0;

  // Thermalize
  for(int n=0;n<N_therm*10;n++){
    update_susy();
  }

  // Calculate Wavefunction
  for(int n=0;n<N_meas;n++){
    // Decorrelate
    for(int m=0;m<N_skip*2;m++){ update_susy(); }
      for(int i=0;i<N_lat;i++){
      if(abs(x[i])<x_range){
        psi[xtoi(x[i])]+=1;
      }
    }
  }


  printf("2nd round acceptance rate: %f \n", (double)N_acc/N_tot);

  // Plot wavefunction
  for(int i=0;i<N_bin;i++){
    fprintf(fpsi,"%f, %f\n", itox(i), sqrt(psi[i]/(dx*N_meas*N_lat)));
  }

  cleanup();
  return 0;
}

// Convert x into bin index
int xtoi(double x){
  return nearbyint((x+x_range)/dx);
}

// Convert bin index to x
double itox(int i){
  return i*dx-x_range;
}

// Calculate difference in action
double dS(int i, double x_old){
  return (x[i]*(x[i]-x[(i+1)%N_lat]-x[(i-1+N_lat)%N_lat])/a
          + alpha*a*sq(x[i]) + lambda*a*qd(x[i]))
          - (x_old*(x_old-x[(i+1)%N_lat]-x[(i-1+N_lat)%N_lat])/a
               + alpha*a*sq(x_old) + lambda*a*qd(x_old));
}

double dS_susy(int i, double x_old){
  return (x[i]*(x[i]-x[(i+1)%N_lat]-x[(i-1+N_lat)%N_lat])/a
          + a*H_susy[xtoi(x[i])])
          - (x_old*(x_old-x[(i+1)%N_lat]-x[(i-1+N_lat)%N_lat])/a
               + a*H_susy[xtoi(x_old)]);
}

void update(){
  double x_old;
  double d; // Difference in action
  for(int i=0;i<N_lat;i++){
    x_old=x[i];
    x[i]=x[i] + epsilon * (gsl_rng_uniform(r)-0.5)*2;
    N_acc++; N_tot++;
    d=dS(i,x_old);
    if(d>0){
      if(exp(-d)<gsl_rng_uniform(r)){ x[i]=x_old; N_acc--; } // Revert to old value
    }
  }
}


void update_susy(){
  double x_old;
  double d; // Difference in action
  for(int i=0;i<N_lat;i++){
    x_old=x[i];
    x[i]=x[i] + epsilon_susy * (gsl_rng_uniform(r)-0.5)*2;
    N_acc++; N_tot++;
    d=dS_susy(i,x_old);
    if(d>0){
      if(exp(-d)<gsl_rng_uniform(r)){ x[i]=x_old; N_acc--; } // Revert to old value
    }
  }
}

void setup(){
  // Cold start
  for(int i=0;i<N_lat;i++){ x[i]=0; C[i]=0; }

  // Initialize variables
  N_acc = 0;
  N_tot = 0;
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

  // Setup spline
  spline_accel = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_cspline, N_spl);
}

void cleanup(){
  gsl_rng_free (r);

  gsl_spline_free(spline);
  gsl_interp_accel_free(spline_accel);

  // Close Files
  fclose(fE0);
  fclose(fpsi0);
  fclose(fpsi);
}
