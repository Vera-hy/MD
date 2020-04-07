/*
 *  Simple molecular dynamics code.
 *
 */
#include <stdio.h>
#include <math.h>
#include "coord.h"
#include "util.h"

//void vis_force(int N,double *f, double *vis, double *vel);
//void vis_force(int N,double *f, double vis, double *vel);
//void add_norm(int N,double *r, double *delta);
//double force(double W, double delta, double r);
//void wind_force(int N,double *f, double *vis, double vel);
//void wind_force(int N,double *f, double vis, double *wind);
//void visc_wind_force(int N,double *f, double vis, double *velo, double *wind);


void evolve(int count,double dt) {
    int step;
    int i, j, k, l;
    //int collided;
    double Size;
/*
 * Loop over timesteps.
 */
    int cM_multi_G = M_central * G;
    for (step = 1; step <= count; step++) {
        printf("timestep %d\n", step);
        printf("collisions %d\n", collisions);
        for (k = 0; k < Nbody; k++) {
/* set the viscosity term and wind term in the force calculation */
            visc_wind_force(Ndim, f[k], vis[k], velo[k], wind);
/* calculate distance from central mass */
            r[k] = 0.0;
            add_norm(Ndim, &r[k], pos[k]);
            r[k] = sqrt(r[k]);
/* calculate central force */
            for (l = 0; l < Ndim; l++) {
                //f[k][l] -= force(G * mass[k] * M_central, pos[k][l], r[k]);
                f[k][l] -= force(mass[k] * cM_multi_G, pos[k][l], r[k]);
            }
        }

/*
    * add pairwise forces.
    */
        k = 0;
        for(i=0; i<Nbody; i++) {
            for(j=i+1; j<Nbody; j++) {
                double G_multi_mSquare = G*mass[i]*mass[j];

                /* calculate pairwise separation of particles */
                for(l=0; l<Ndim; l++) {
                    delta_pos[k][l] = pos[i][l] - pos[j][l];
                }

                /* calculate norm of separation vector */
                delta_r[k] = 0.0;
                add_norm(Ndim, &delta_r[k], delta_pos[k]);
                delta_r[k] = sqrt(delta_r[k]);

                Size = radius[i] + radius[j];
               // collided=0;

                /*  flip force if close in */
                if( delta_r[k] >= Size ) {
                    for(l=0; l<Ndim; l++) {
                        //f[i][l] -= force(G*mass[i]*mass[j],delta_pos[k][l],delta_r[k]);
                        //f[j][l] += force(G*mass[i]*mass[j],delta_pos[k][l],delta_r[k]);
                       // f[i][l] -= force(G_multi_mSquare,delta_pos[k][l],delta_r[k]);
                       // f[j][l] += force(G_multi_mSquare,delta_pos[k][l],delta_r[k]);
                       double calc_force = force(G_multi_mSquare,delta_pos[k][l],delta_r[k]);
                        f[i][l] -= calc_force;
                        f[j][l] += calc_force;
                    }
                }
                else {
                    for(l=0; l<Ndim; l++) {
                        //f[i][l] += force(G*mass[i]*mass[j],delta_pos[k][l],delta_r[k]);
                        //f[j][l] -= force(G*mass[i]*mass[j],delta_pos[k][l],delta_r[k]);
                        //f[i][l] += force(G_multi_mSquare,delta_pos[k][l],delta_r[k]);
                        //f[j][l] -= force(G_multi_mSquare,delta_pos[k][l],delta_r[k]);
                        double calc_force = force(G_multi_mSquare,delta_pos[k][l],delta_r[k]);
                        f[i][l] += calc_force;
                        f[j][l] -= calc_force;
                       // collided=1;
                    }
                    collisions++;
                }

//                if( collided == 1 ){
//                    collisions++;
//                }
                k = k + 1;
            }
        }

/* update positions and velocities */
        //for (i = 0; i < Nbody; i++) {
        for(i=Nbody-1; i>=0; i--) {
            for (j = 0; j < Ndim; j++) {
                pos[i][j] += +dt * velo[i][j];
                velo[i][j] += +dt * (f[i][j] / mass[i]);
            }
        }
    }

}







