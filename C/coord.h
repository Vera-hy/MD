/*
 * This file defines static arrays that contains the primary coordinates
 * of the particles,
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
// Declare arrays
/*
DEF double pos[Nbody][Ndim]__attribute__((aligned(64)));
DEF double velo[Nbody][Ndim]__attribute__((aligned(64)));
DEF double f[Nbody][Ndim]__attribute__((aligned(64)));
DEF double vis[Nbody]__attribute__((aligned(64)));
DEF double mass[Nbody]__attribute__((aligned(64)));
DEF double radius[Nbody]__attribute__((aligned(64)));
DEF double delta_pos[Npair][Ndim]__attribute__((aligned(64)));
DEF double r[Nbody]__attribute__((aligned(64)));
DEF double delta_r[Npair]__attribute__((aligned(64)));
DEF double wind[Ndim]__attribute__((aligned(64)));
DEF int collisions;*/
#define PADDING 64
DEF double pos[Nbody][Ndim] __attribute__((aligned(64)));
DEF double velo[Nbody][Ndim] __attribute__((aligned(64)));
//DEF double f[Nbody][Ndim] __attribute__((aligned(64)));
DEF double vis[Nbody + PADDING] __attribute__((aligned(64)));
DEF double radius[Nbody + PADDING] __attribute__((aligned(64)));
DEF double mass[Nbody + PADDING] __attribute__((aligned(64)));
//DEF double delta_pos[Npair][Ndim] __attribute__((aligned(64)));
//DEF double r[Nbody + PADDING] __attribute__((aligned(64)));
//DEF double delta_r[Npair + PADDING] __attribute__((aligned(64)));
DEF double wind[Ndim + PADDING] __attribute__((aligned(64)));
DEF int collisions;

#define G 2.0
#define M_central 1000.0
#define cM_multi_G 2000.0

void evolve(int Nstep, double dt);

inline void visc_wind_force(int N,double *f, double vis, double *velo, double *wind)
{
    int j;
    for(j=0;j<N;j++){
        //f[i] = -vis[i] * velo[i];
        f[j] = -vis * (velo[j] + wind[j]);
    }
}
/*inline double add_norm(int N, double *delta)
{
    int k;
    double dist_sum;
    dist_sum = 0;
    for(k=0;k<N;k++){
        dist_sum += (delta[k] * delta[k]);
    }
   // double r = sqrt(r_sum);
    return dist_sum;
}*/

inline double force(double W, double delta, double r){
    return W*delta/(pow(r,3.0));
}
