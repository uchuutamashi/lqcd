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


// SU(3) Matrix
// only the first two rows carries useful information (Lang p.80)
// note that the first two rows must be orthonormal
// use su3_normalize() for this
typedef struct{
  gsl_complex a1;
  gsl_complex a2;
  gsl_complex a3;
  gsl_complex b1;
  gsl_complex b2;
  gsl_complex b3;
} su3_matrix;

// SU(3) matrix methods --------------------------------------

// print matrix
void su3_print(const su3_matrix A){
  printf("%f + i*%f \t", GSL_REAL(A.a1), GSL_IMAG(A.a1));
  printf("%f + i*%f \t", GSL_REAL(A.a2), GSL_IMAG(A.a2));
  printf("%f + i*%f \n", GSL_REAL(A.a3), GSL_IMAG(A.a3));
  printf("%f + i*%f \t", GSL_REAL(A.b1), GSL_IMAG(A.b1));
  printf("%f + i*%f \t", GSL_REAL(A.b2), GSL_IMAG(A.b2));
  printf("%f + i*%f \n", GSL_REAL(A.b3), GSL_IMAG(A.b3));
}

// the identity matrix
su3_matrix su3_unit(){
  su3_matrix u;
  u.a1=GSL_COMPLEX_ONE;
  u.a2=GSL_COMPLEX_ZERO;
  u.a3=GSL_COMPLEX_ZERO;
  u.b1=GSL_COMPLEX_ZERO;
  u.b2=GSL_COMPLEX_ONE;
  u.b3=GSL_COMPLEX_ZERO;
  return u;
}

// normalize matrix
su3_matrix su3_norm(const su3_matrix A){
  su3_matrix result;
  // normalize first row
  double norm = sqrt(gsl_complex_abs2(A.a1)+gsl_complex_abs2(A.a2)+gsl_complex_abs2(A.a3));
  result.a1 = gsl_complex_div(A.a1,gsl_complex_rect(norm,0));
  result.a2 = gsl_complex_div(A.a2,gsl_complex_rect(norm,0));
  result.a3 = gsl_complex_div(A.a3,gsl_complex_rect(norm,0));

  // orthogonize second row
  gsl_complex dprod = gsl_complex_add(gsl_complex_add(gsl_complex_mul(A.b1,gsl_complex_conjugate(result.a1)),
                      gsl_complex_mul(A.b2,gsl_complex_conjugate(result.a2))),
                      gsl_complex_mul(A.b3,gsl_complex_conjugate(result.a3)));
  result.b1 = gsl_complex_sub(A.b1,gsl_complex_mul(result.a1,dprod));
  result.b2 = gsl_complex_sub(A.b2,gsl_complex_mul(result.a2,dprod));
  result.b3 = gsl_complex_sub(A.b3,gsl_complex_mul(result.a3,dprod));

  // normalize second row
  norm = sqrt(gsl_complex_abs2(result.b1)+gsl_complex_abs2(result.b2)+gsl_complex_abs2(result.b3));
  result.b1 = gsl_complex_div(result.b1,gsl_complex_rect(norm,0));
  result.b2 = gsl_complex_div(result.b2,gsl_complex_rect(norm,0));
  result.b3 = gsl_complex_div(result.b3,gsl_complex_rect(norm,0));
  return result;
}

// multiply matrices
su3_matrix su3_mul(const su3_matrix A, const su3_matrix B){
  su3_matrix result;
  gsl_complex w1,w2,w3;

  // Third row of B as a cross product (c.f. Lang p.80)
  w1 = gsl_complex_sub( gsl_complex_mul(gsl_complex_conjugate(B.a2),gsl_complex_conjugate(B.b3)), gsl_complex_mul(gsl_complex_conjugate(B.a3),gsl_complex_conjugate(B.b2)) );
  w2 = gsl_complex_sub( gsl_complex_mul(gsl_complex_conjugate(B.a3),gsl_complex_conjugate(B.b1)), gsl_complex_mul(gsl_complex_conjugate(B.a1),gsl_complex_conjugate(B.b3)) );
  w3 = gsl_complex_sub( gsl_complex_mul(gsl_complex_conjugate(B.a1),gsl_complex_conjugate(B.b2)), gsl_complex_mul(gsl_complex_conjugate(B.a2),gsl_complex_conjugate(B.b1)) );

  result.a1 = gsl_complex_add(gsl_complex_add(gsl_complex_mul(A.a1,B.a1),gsl_complex_mul(A.a2,B.b1)), gsl_complex_mul(A.a3,w1));
  result.a2 = gsl_complex_add(gsl_complex_add(gsl_complex_mul(A.a1,B.a2),gsl_complex_mul(A.a2,B.b2)), gsl_complex_mul(A.a3,w2));
  result.a3 = gsl_complex_add(gsl_complex_add(gsl_complex_mul(A.a1,B.a3),gsl_complex_mul(A.a2,B.b3)), gsl_complex_mul(A.a3,w3));

  result.b1 = gsl_complex_add(gsl_complex_add(gsl_complex_mul(A.b1,B.a1),gsl_complex_mul(A.b2,B.b1)), gsl_complex_mul(A.b3,w1));
  result.b2 = gsl_complex_add(gsl_complex_add(gsl_complex_mul(A.b1,B.a2),gsl_complex_mul(A.b2,B.b2)), gsl_complex_mul(A.b3,w2));
  result.b3 = gsl_complex_add(gsl_complex_add(gsl_complex_mul(A.b1,B.a3),gsl_complex_mul(A.b2,B.b3)), gsl_complex_mul(A.b3,w3));

  return result;
}

// inverse matrix
su3_matrix su3_inv(const su3_matrix A){
  su3_matrix result;
  gsl_complex w1,w2,w3;

  // Third row of B as a cross product (c.f. Lang p.80)
  w1 = gsl_complex_sub( gsl_complex_mul(gsl_complex_conjugate(A.a2),gsl_complex_conjugate(A.b3)), gsl_complex_mul(gsl_complex_conjugate(A.a3),gsl_complex_conjugate(A.b2)) );
  w2 = gsl_complex_sub( gsl_complex_mul(gsl_complex_conjugate(A.a3),gsl_complex_conjugate(A.b1)), gsl_complex_mul(gsl_complex_conjugate(A.a1),gsl_complex_conjugate(A.b3)) );
  w3 = gsl_complex_sub( gsl_complex_mul(gsl_complex_conjugate(A.a1),gsl_complex_conjugate(A.b2)), gsl_complex_mul(gsl_complex_conjugate(A.a2),gsl_complex_conjugate(A.b1)) );

  result.a1 = gsl_complex_sub( gsl_complex_mul(A.b2,w3), gsl_complex_mul( A.b3,w2) );
  result.a2 = gsl_complex_sub( gsl_complex_mul(A.a3,w2), gsl_complex_mul( A.a2,w3) );
  result.a3 = gsl_complex_sub( gsl_complex_mul(A.a2,A.b3), gsl_complex_mul( A.a3,A.b2) );

  result.b1 = gsl_complex_sub( gsl_complex_mul(A.b3,w1), gsl_complex_mul( A.b1,w3) );
  result.b2 = gsl_complex_sub( gsl_complex_mul(A.a1,w3), gsl_complex_mul( A.a3,w1) );
  result.b3 = gsl_complex_sub( gsl_complex_mul(A.a3,A.b1), gsl_complex_mul( A.a1,A.b3) );

  return result;
}

// generate random matrix
su3_matrix su3_rand(const gsl_rng* r){
  su3_matrix result;
  result.a1 = gsl_complex_rect(1,eps*(gsl_rng_uniform(r)*2-1));
  result.a2 = gsl_complex_rect(0,eps*(gsl_rng_uniform(r)*2-1));
  result.a3 = gsl_complex_rect(0,eps*(gsl_rng_uniform(r)*2-1));

  result.b1 = gsl_complex_rect(0,eps*(gsl_rng_uniform(r)*2-1));
  result.b2 = gsl_complex_rect(1,eps*(gsl_rng_uniform(r)*2-1));
  result.b3 = gsl_complex_rect(0,eps*(gsl_rng_uniform(r)*2-1));
  return su3_norm(result);
}

double su3_trace(const su3_matrix A){
  gsl_complex c3;
  c3 = gsl_complex_sub( gsl_complex_mul(gsl_complex_conjugate(A.a1),gsl_complex_conjugate(A.b2)), gsl_complex_mul(gsl_complex_conjugate(A.a2),gsl_complex_conjugate(A.b1)) );

  return GSL_REAL(gsl_complex_add(gsl_complex_add(A.a1,A.b2),c3));
}

// end of SU(3) matrix methods ---------------------------------------------
