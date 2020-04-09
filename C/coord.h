/*
 * This file declares function evolve and two inlined functions.
 *
 *  Nbody	  Number of particles
 *  Npair	  Number of particle pairs
 *  pos		  Position of the particles
 *  r         distance of partice from central mass 
 *  vel		  velocity of the particles
 *  f		  Forces acting on each particle
 *  vis       viscosity coefficient for each particle
 *  mass	  mass of each particle
 *  delta_pos separation vector for each particle pair
 *  delta_r	  separation for each particle pair
 */
#include <math.h>

#ifdef DECL
#define DEF
#else
#define DEF extern
#endif
#define Nbody 4*1024
//#define  Npair ((Nbody*(Nbody-1))/2)

enum{ Xcoord=0, Ycoord, Zcoord, Ndim };
DEF int collisions;
#define PADDING 64
#define G 2.0
#define M_central 1000.0
#define cM_multi_G 2000.0

void evolve(int Nstep, double dt, double pos[Nbody][Ndim],double velo[Nbody][Ndim],
        double * restrict vis,double* restrict radius,double* restrict mass,double* restrict wind);

inline void outside_force(int N,double *f, double vis, double *velo, double *wind)
{
    int j;
    #pragma ivdep
    //#pragma omp simd
    #pragma vector aligned
    for(j=0;j<N;j++){
        //f[i] = -vis[i] * velo[i];
        f[j] = -vis * (velo[j] + wind[j]);
    }
}

inline double force(double W, double delta, double r){
    return W*delta/(pow(r,3.0));
}
