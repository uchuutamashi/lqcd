#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>

#include "su3_matrix.c"

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
double plaquette(const su3_matrix* links, int mu, int nu, int* site){
  su3_matrix A,B,C,D;
  A = links[index(site,mu)];
  shift(site,mu,FORWARD);
  B = links[index(site,nu)];
  shift(site,mu,BACKWARD);
  shift(site,nu,FORWARD);
  C = su3_inv(links[index(site,mu)]);
  shift(site,nu,BACKWARD); // restore site position
  D = su3_inv(links[index(site,nu)]);

  A = su3_mul(A,B);
  A = su3_mul(A,C);
  A = su3_mul(A,D);

  return su3_trace(A)/3.0;
}

// plaquette spatial mean
double plaq_mean_s(int t, const su3_matrix* links){
  double mean = 0;
  int site[4];
  for(int z=0;z<N;z++){
    for(int y=0;y<N;y++){
      for(int x=0;x<N;x++){
          site[0]=t; site[1]=x; site[2]=y; site[3]=z;
          for(int mu=1;mu<4;mu++){
            for(int nu=mu+1;nu<4;nu++){
              mean += plaquette(links,mu,nu,site);
            }
          }
      }
    }
  }
  return mean/ (3*N*N*N); // 3 directions & N^3 sites
}


// end of plaquette methods ------------------------------------------------

// Metropolis update --------------------------------------------------------

void do_hb_update(su3_matrix* links, const gsl_rng* r, const su3_matrix* rands, int* site, int dir);
double action(su3_matrix* links, int* site, int dir);

// the big loop
void hb_update(su3_matrix* links, const gsl_rng* r, const su3_matrix* rands){
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
void do_hb_update(su3_matrix* links, const gsl_rng* r, const su3_matrix* rands, int* site, int dir){
  if(verbose) {printf("Updating {%d, %d, %d, %d} direction %d",site[0],site[1],site[2],site[3],dir);getchar();}
  double dS,action_old;
  su3_matrix link_old;
  for(int times=0; times<Nk; times++){
    action_old = action(links,site,dir);
    link_old = links[index(site,dir)];
    links[index(site,dir)]=su3_mul(links[index(site,dir)],rands[gsl_rng_uniform_int(r, 50)]);

    dS = action(links,site,dir)-action_old;

    if(verbose) {printf("Old action = %f, ",action_old); su3_print(link_old); getchar();}
    if(verbose) {printf("New action = %f, ",action_old+dS); su3_print(links[index(site,dir)]); getchar();}
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
double action(su3_matrix* links, int* site, int dir){
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

// Loops -------------------------------------------------------------------
// Generic loop ------------------------------------------------------------

double loop(int r, int t, int mu, int nu, const su3_matrix* links, int* site){
  su3_matrix A[(r+t)*2];

  // bottom side
  for(int i=0;i<t;i++){
    A[i]=links[index(site,nu)];
    shift(site,nu,FORWARD);
  }

  // right side
  for(int i=0;i<r;i++){
    A[t+i]=links[index(site,mu)];
    shift(site,mu,FORWARD);
  }

  // top side
  for(int i=0;i<t;i++){
    shift(site,nu,BACKWARD);
    A[t+r+i]=su3_inv(links[index(site,nu)]);
  }

  // left side
  for(int i=0;i<r;i++){
    shift(site,mu,BACKWARD);
    A[t+r+t+i]=su3_inv(links[index(site,mu)]);
  }

  // multiply
  for(int i=1;i<2*(r+t);i++){
    A[0]=su3_mul(A[0],A[i]);
  }

  return su3_trace(A[0])/3.0;
}

// end of Loop -------------------------------------------------------------
// Wilson loop -------------------------------------------------------------

double wloop(int r, int t, int mu,const su3_matrix* links, int* site){
  return loop(r, t, mu, 0, links, site);
}

double wloop_mean(int wR, int wT, const su3_matrix* links){
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
// Polyakov loop -----------------------------------------------------------
double polyakov(const su3_matrix* links, int* site){
  su3_matrix A[N];

  for(int i=0;i<N;i++){
    A[i]=links[index(site,0)];
    shift(site,0,FORWARD);
  }

  // multiply
  for(int i=1;i<N;i++){
    A[0]=su3_mul(A[0],A[i]);
  }

  return su3_trace(A[0])/3.0;
}

double polyakov_mean(const su3_matrix* links){
  double mean = 0;
  int site[4];
  site[0]=0;
  for(int z=0;z<N;z++){
    for(int y=0;y<N;y++){
      for(int x=0;x<N;x++){
        site[1]=x; site[2]=y; site[3]=z;
        mean += polyakov(links,site);
      }
    }
  }
  return mean/ (N*N*N);
}
// end of Polyakov loop ----------------------------------------------------

// Twisted Wilson loop -----------------------------------------------------
// Ref: Lepage p.33

// r=sqrt(2)
double twloop2(int t, int mu, int nu ,const su3_matrix* links, int* site){
  su3_matrix A[(2+t)*2];

  // bottom side
  for(int i=0;i<t;i++){
    A[i]=links[index(site,0)];
    shift(site,0,FORWARD);
  }

  // right side, twisted
  A[t+0]=links[index(site,mu)];
  shift(site,mu,FORWARD);
  A[t+1]=links[index(site,nu)];
  shift(site,nu,FORWARD);


  // top side
  for(int i=0;i<t;i++){
    shift(site,0,BACKWARD);
    A[t+2+i]=su3_inv(links[index(site,0)]);
  }

  // left side, twisted
  shift(site,mu,BACKWARD);
  A[t+2+t+0]=su3_inv(links[index(site,mu)]);
  shift(site,nu,BACKWARD);
  A[t+2+t+1]=su3_inv(links[index(site,nu)]);

  // multiply
  for(int i=1;i<2*(2+t);i++){
    A[0]=su3_mul(A[0],A[i]);
  }

  return su3_trace(A[0])/3.0;
}

double twloop2_mean(int wT, const su3_matrix* links){
  double mean = 0;
  int site[4];
  for(int z=0;z<N;z++){
    for(int y=0;y<N;y++){
      for(int x=0;x<N;x++){
        for(int t=0;t<N;t++){
            site[0]=t; site[1]=x; site[2]=y; site[3]=z;
            for(int mu=1;mu<4;mu++){
              for(int nu=1;nu<4;nu++){
                if(mu!=nu){
                  mean += twloop2(wT,mu,nu,links,site);
                }
              }
            }
        }
      }
    }
  }
  return mean/ (6*N*N*N*N); // 3 x 2 directions & N^3 sites
}

// r=sqrt(3)
double twloop3(int t, const su3_matrix* links, int* site){
  su3_matrix A[(3+t)*2];

  // bottom side
  for(int i=0;i<t;i++){
    A[i]=links[index(site,0)];
    shift(site,0,FORWARD);
  }

  // right side, twisted
  A[t+0]=links[index(site,1)];
  shift(site,1,FORWARD);
  A[t+1]=links[index(site,2)];
  shift(site,2,FORWARD);
  A[t+2]=links[index(site,3)];
  shift(site,3,FORWARD);


  // top side
  for(int i=0;i<t;i++){
    shift(site,0,BACKWARD);
    A[t+3+i]=su3_inv(links[index(site,0)]);
  }

  // left side, twisted
  shift(site,1,BACKWARD);
  A[t+3+t+0]=su3_inv(links[index(site,1)]);
  shift(site,2,BACKWARD);
  A[t+3+t+1]=su3_inv(links[index(site,2)]);
  shift(site,3,BACKWARD);
  A[t+3+t+2]=su3_inv(links[index(site,3)]);

  // multiply
  for(int i=1;i<2*(3+t);i++){
    A[0]=su3_mul(A[0],A[i]);
  }

  return su3_trace(A[0])/3.0;
}

double twloop3_mean(int wT, const su3_matrix* links){
  double mean = 0;
  int site[4];
  for(int z=0;z<N;z++){
    for(int y=0;y<N;y++){
      for(int x=0;x<N;x++){
        for(int t=0;t<N;t++){
            site[0]=t; site[1]=x; site[2]=y; site[3]=z;
            mean += twloop3(wT,links,site);
        }
      }
    }
  }
  return mean/ (N*N*N*N); // N^3 sites
}
// end of twisted Wilson loop --------------------------------------------

// Plaquette correlator
double plaq_corr(int T, const su3_matrix* links){
  double mean=0;
  for(int i=0;i<N;i++){
    mean+=plaq_mean_s(i,links)*plaq_mean_s(i+T,links);
  }
  return mean/N;
}
// end of plaquette correlator

// TODO: vvvvvvvv


// Smeared Wilson loop ---------------------------------------------------


// end of smeared Wilson loop --------------------------------------------
