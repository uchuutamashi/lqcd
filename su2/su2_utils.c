#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>

#include "su2_matrix.c"

#ifndef gauge_config
#include "config.h"
#define gauge_config 1

#define sq(x) ((x)*(x))
#define cb(x) ((x)*(x)*(x))
#endif

int acc, tot; // acceptance rate

#define FORWARD 1 //true
#define BACKWARD 0 //false

// Access pattern methods --------------------------------------------------

// index of link site=(t,x,y,z) in direction dir
int index(const int* site, int dir){
  return dir+site[0]*4+site[1]*N*4+site[2]*N*N*4+site[3]*N*N*N*4;
}

// shift site array by one to direction dir
void shift(int* site, int dir, int forward){
  if(forward == 1){
    site[dir] = (site[dir]+1)%N;
  }else{
    site[dir] = (site[dir]-1+N)%N;
  }
}

// end of access pattern methods -------------------------------------------

// Plaquette methods -------------------------------------------------------

// plaquette operator at site (t,x,y,z) in plane (mu,nu)
// normally mu>nu
double plaquette(const su2_matrix* links, int mu, int nu, int* site){
  su2_matrix A,B,C,D;
  A = links[index(site,mu)];
  shift(site,mu,FORWARD);
  B = links[index(site,nu)];
  shift(site,mu,BACKWARD);
  shift(site,nu,FORWARD);
  C = su2_inv(links[index(site,mu)]);
  shift(site,nu,BACKWARD); // restore site position
  D = su2_inv(links[index(site,nu)]);

  A = su2_mul(A,B);
  A = su2_mul(A,C);
  A = su2_mul(A,D);

  return su2_trace(A)/2.0;
}

double plaquette_mean(const su2_matrix* links, int mu, int nu){
  double mean = 0;
  int site[4];
  for(int z=0;z<N;z++){
    for(int y=0;y<N;y++){
      for(int x=0;x<N;x++){
        for(int t=0;t<N;t++){
            site[0]=t; site[1]=x; site[2]=y; site[3]=z;
            mean += plaquette(links,mu,nu,site);
        }
      }
    }
  }
  return mean / (N*N*N*N);
}


// end of plaquette methods ------------------------------------------------

// Metropolis update --------------------------------------------------------

void do_hb_update(su2_matrix* links, const gsl_rng* r, const su2_matrix* rands, int* site, int dir);
double action(su2_matrix* links, int* site, int dir);

// the big loop
void hb_update(su2_matrix* links, const gsl_rng* r, const su2_matrix* rands){
  int site[4];
  for(int z=0;z<N;z++){
    for(int y=0;y<N;y++){
      for(int x=0;x<N;x++){
        for(int t=0;t<N;t++){
          for(int dir=0;dir<4;dir++){
            site[0]=t; site[1]=x; site[2]=y; site[3]=z;
            do_hb_update(links,r,rands,site,dir);
          }
        }
      }
    }
  }
}

// update
void do_hb_update(su2_matrix* links, const gsl_rng* r, const su2_matrix* rands, int* site, int dir){
  if(verbose) {printf("Updating {%d, %d, %d, %d} direction %d",site[0],site[1],site[2],site[3],dir);getchar();}
  double dS,action_old;
  su2_matrix link_old;
  for(int times=0; times<Nk; times++){
    action_old = action(links,site,dir);
    link_old = links[index(site,dir)];
    links[index(site,dir)]=su2_mul(links[index(site,dir)],rands[gsl_rng_uniform_int(r, 50)]);

    dS = action(links,site,dir)-action_old;

    if(verbose) {printf("Old action = %f, ",action_old); su2_print(link_old); getchar();}
    if(verbose) {printf("New action = %f, ",action_old+dS); su2_print(links[index(site,dir)]); getchar();}
    acc++;tot++;
    if(dS>0){
      if(exp(-dS)<gsl_rng_uniform(r)){
        acc--;
        links[index(site,dir)]=link_old; // revert update
        if(verbose) {printf("Reverted"); getchar();}
      }
    }
  }

}

// calculating the action
// for each direction we need to find plaquette at 2 sites in each of the 3 planes
// mu > nu
double action(su2_matrix* links, int* site, int dir){
  int mu,nu;
  double result=0;
  for(int k=0; k<4; k++){
    if(k==dir) continue; // loop over planes
    if(k>dir){mu=k; nu=dir;} // mu=larger number
    if(k<dir){mu=dir; nu=k;}

    result += plaquette(links, mu,nu,site);
    shift(site,k,BACKWARD); // shift to the perpendicular direction backwards
    result += plaquette(links, mu,nu,site);
    shift(site,k,FORWARD);  // shift it back
  }

  return -beta*result;
}

// end of heat bath update -------------------------------------------------

// Wilson loop -------------------------------------------------------------

double wloop(int r, int t, int mu,const su2_matrix* links, int* site){
  su2_matrix A[(r+t)*2];

  // bottom side
  for(int i=0;i<t;i++){
    A[i]=links[index(site,0)];
    shift(site,0,FORWARD);
  }

  // right side
  for(int i=0;i<r;i++){
    A[t+i]=links[index(site,mu)];
    shift(site,mu,FORWARD);
  }

  // top side
  for(int i=0;i<t;i++){
    shift(site,0,BACKWARD);
    A[t+r+i]=su2_inv(links[index(site,0)]);
  }

  // left side
  for(int i=0;i<r;i++){
    shift(site,mu,BACKWARD);
    A[t+r+t+i]=su2_inv(links[index(site,mu)]);
  }

  // multiply
  for(int i=1;i<2*(r+t);i++){
    A[0]=su2_mul(A[0],A[i]);
  }

  return su2_trace(A[0])/2.0;
}

double wloop_mean(int wR, int wT, const su2_matrix* links){
  double mean = 0;
  int site[4];
  for(int z=0;z<N;z++){
    for(int y=0;y<N;y++){
      for(int x=0;x<N;x++){
        for(int t=0;t<N;t++){
            site[0]=t; site[1]=x; site[2]=y; site[3]=z;
            for(int mu=1;mu<4;mu++){
              mean += wloop(wR,wT,mu,links,site);
            }
        }
      }
    }
  }
  return mean/ (3*N*N*N*N); // 3 directions & N^3 sites
}

// end of Wilson loop ------------------------------------------------------
