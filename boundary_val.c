#include "boundary_val.h"
#include <stdio.h>

void boundaryvalues(int imax,int jmax,double** U,double** V){
    for (int j=1;j<imax;j++){
        u[0][j] =0;
        u[imax][j]=0;
        v[0][j] = -v[1][j];
        v[imax+1][j] = -v[imax][j];
        p[0][j] = p[1][j];
        p[imax+1][j] = p[imax][j];
        F[0][j] = u[0][j];
        F[imax][j] = umax[i][j];
    }
    for (int i=1;i<imax;i++){
        v[i][0]=0;
        v[i][jmax]=0;
        u[i][0]=-u[i][1];
        u[i][jmax+1] = -u[i][jmax];
        p[i][0] = p[i][1];
        p[i][jmax+1]=p[i][jmax];
        G[i][0] = v[i][0];
        G[i][jmax]=v[i][jmax];
    }
    
}
