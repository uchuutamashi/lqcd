#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>

#ifndef gauge_config
#include "config.h"
#define gauge_config 1

#define sq(x) ((x)*(x))
#define cb(x) ((x)*(x)*(x))
#endif


// SU(2) Matrix
// only the first row carries useful information (Lang p.80)
// note that for the first row must be normalized, i.e. |a|^2+ |b|^2 = 1
// use su2_normalize() for this
typedef struct{
  gsl_complex a;
  gsl_complex b;
} su2_matrix;

// SU(2) matrix methods --------------------------------------

// print matrix
void su2_print(const su2_matrix A){
  printf("%f + i*%f \t", GSL_REAL(A.a), GSL_IMAG(A.a));
  printf("%f + i*%f \n", GSL_REAL(A.b), GSL_IMAG(A.b));
}

// the identity matrix
su2_matrix su2_unit(){
  su2_matrix u;
  u.a=GSL_COMPLEX_ONE;
  u.b=GSL_COMPLEX_ZERO;
  return u;
}

// normalize matrix
su2_matrix su2_norm(const su2_matrix A){
  su2_matrix result;
  double norm = sqrt(gsl_complex_abs2(A.a)+gsl_complex_abs2(A.b));
  result.a = gsl_complex_div(A.a,gsl_complex_rect(norm,0));
  result.b = gsl_complex_div(A.b,gsl_complex_rect(norm,0));
  return result;
}

// multiply matrices
su2_matrix su2_mul(const su2_matrix A, const su2_matrix B){
  su2_matrix result;
  result.a = gsl_complex_sub(gsl_complex_mul(A.a,B.a), gsl_complex_mul(A.b,gsl_complex_conjugate(B.b))); // ac-bd*
  result.b = gsl_complex_add(gsl_complex_mul(A.a,B.b), gsl_complex_mul(A.b,gsl_complex_conjugate(B.a))); // ad+bc*
  return result;
}

// inverse matrix
su2_matrix su2_inv(const su2_matrix A){
  su2_matrix result;
  result.a = gsl_complex_conjugate(A.a);
  result.b = gsl_complex_sub(GSL_COMPLEX_ZERO,A.b);
  return result;
}

// generate random matrix
su2_matrix su2_rand(const gsl_rng* r){
  su2_matrix result;
  result.a = gsl_complex_rect(1,eps*(gsl_rng_uniform(r)*2-1));
  result.b = gsl_complex_rect(0,eps*(gsl_rng_uniform(r)*2-1));
  return su2_norm(result);
}

double su2_trace(const su2_matrix A){
  return 2*GSL_REAL(A.a);
}

// end of SU(2) matrix methods ---------------------------------------------
