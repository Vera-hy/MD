/*
 *  Simple molecular dynamics code.
 *
 */
#include <stdio.h>
#include <math.h>
#include "coord.h"

/* In C, use file scope globals in preference to externals if possible */
static double f[Nbody][Ndim] __attribute__((aligned(64)));
static double delta_pos[Ndim + PADDING] __attribute__((aligned(64)));

void evolve(int count,double dt,double pos[Nbody][Ndim],double velo[Nbody][Ndim],
        double * restrict vis,double* restrict radius,double* restrict mass,double* restrict wind) {
    int step;
    int i, j, k, l;
    double Size;
/*
 * Loop over timesteps.
 */
    for (step = 1; step <= count; step++) {
        printf("timestep %d\n", step);
        printf("collisions %d\n", collisions);

        for (k = 0; k < Nbody; k++) {
/* set the viscosity term and wind term in the force calculation */
            outside_force(Ndim, f[k], vis[k], velo[k], wind);
/* calculate distance from central mass */
            double r = 0;
#pragma ivdep
#pragma omp simd
#pragma vector aligned
            for (l = 0; l < Ndim; l++) {
                r += pos[k][l] * pos[k][l];
            }
            r = sqrt(r);
#pragma ivdep
#pragma omp simd
#pragma vector aligned
/* add central force */
            for (l = Ndim -1 ; l >= 0 ; l--) {
                f[k][l] -= force(mass[k] * cM_multi_G, pos[k][l], r);
            }
        }

/* add pairwise forces */
        for(i=0; i<Nbody; i++) {
            for(j=i+1; j<Nbody; j++) {
                double G_multi_mSquare = G*mass[i]*mass[j];
                double delta_r = 0;
#pragma ivdep
#pragma omp simd
#pragma vector aligned
                /* calculate pairwise separation of particles */
                for(l=0; l<Ndim; l++) {
                    delta_pos[l] = pos[i][l] - pos[j][l];
                    delta_r += delta_pos[l] * delta_pos[l];
                }
                delta_r = sqrt(delta_r);


                Size = radius[i] + radius[j];

                if( delta_r >= Size ) {
#pragma ivdep
#pragma omp simd
#pragma vector aligned
                    for(l=0; l<Ndim; l++) {
                        double calc_force = force(G_multi_mSquare,delta_pos[l],delta_r);
                        f[i][l] -= calc_force;
                        f[j][l] += calc_force;
                    }
                }
                    /* if two particles are too close, they will collide */
                else {
#pragma ivdep
#pragma omp simd
#pragma vector aligned
                    for(l=0; l<Ndim; l++) {
                        double calc_force = force(G_multi_mSquare,delta_pos[l],delta_r);
                        f[i][l] += calc_force;
                        f[j][l] -= calc_force;
                    }
                    collisions++;
                }
            }
        }

/* update positions and velocities */
        for(i=Nbody-1; i>=0; i--) {
#pragma ivdep
//#pragma omp simd aligned(pos:64,velo:64,f:64,mass:64)
#pragma vector aligned
#pragma omp simd
            for (j = 0; j < Ndim; j++) {
                pos[i][j] += +dt * velo[i][j];
                velo[i][j] += +dt * (f[i][j] / mass[i]);
            }
        }
    }

}







