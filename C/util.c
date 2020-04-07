#include <math.h>

//void visc_force(int N,double *f, double *vis, double *velo)
/*void visc_force(int N,double *f, double vis, double *velo)
{
    int j;
    for(j=0;j<N;j++){
        //f[i] = -vis[i] * velo[i];
        f[j] = -vis * velo[j];
    }
}*/
//void wind_force(int N,double *f, double *vis, double velo)
/*void wind_force(int N,double *f, double vis, double *wind)
{
    int j;
    for(j=0;j<N;j++){
        //f[i] = f[i] -vis[i] * velo;
        f[j] = f[j] -vis * wind[j];
    }
}*/
void visc_wind_force(int N,double *f, double vis, double *velo, double *wind)
{
    int j;
    for(j=0;j<N;j++){
        //f[i] = -vis[i] * velo[i];
        f[j] = -vis * (velo[j] + wind[j]);
    }
}
void add_norm(int N,double *r, double *delta)
{
    int k;
    for(k=0;k<N;k++){
        *r += (delta[k] * delta[k]);
    }
}

double force(double W, double delta, double r){
    return W*delta/(pow(r,3.0));
}



