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

void evolve(int count,double dt) {
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
            visc_wind_force(Ndim, f[k], vis[k], velo[k], wind);
/* calculate distance from central mass */
            double r = 0;
            for (l = 0; l < Ndim; l++) {
                r += pos[k][l] * pos[k][l];
            }
            r = sqrt(r);

            //r = sqrt(add_norm(Ndim, pos[k]));
/* calculate central force */
           //for (l = 0; l < Ndim; l++) {
            for (l = Ndim -1 ; l >= 0 ; l--) {
                f[k][l] -= force(mass[k] * cM_multi_G, pos[k][l], r);
            }
        }

/*
    * add pairwise forces.
    */
        for(i=0; i<Nbody; i++) {
            for(j=i+1; j<Nbody; j++) {
                double G_multi_mSquare = G*mass[i]*mass[j];
                double delta_r = 0;

                /* calculate pairwise separation of particles */
                for(l=0; l<Ndim; l++) {
                    delta_pos[l] = pos[i][l] - pos[j][l];
                    delta_r += delta_pos[l] * delta_pos[l];
                }
                delta_r = sqrt(delta_r);

//                double delta_r;
//                delta_r = sqrt(add_norm(Ndim, delta_pos));

                Size = radius[i] + radius[j];

                if( delta_r >= Size ) {
                    for(l=0; l<Ndim; l++) {
                   // for(l=Ndim - 1; l >= 0; l--) {
                        double calc_force = force(G_multi_mSquare,delta_pos[l],delta_r);
                        f[i][l] -= calc_force;
                        f[j][l] += calc_force;
                    }
                }
                else {
                    for(l=0; l<Ndim; l++) {
                    //for(l=Ndim - 1; l >= 0; l--) {
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
            for (j = 0; j < Ndim; j++) {
                pos[i][j] += +dt * velo[i][j];
                velo[i][j] += +dt * (f[i][j] / mass[i]);
            }
        }
    }

}







