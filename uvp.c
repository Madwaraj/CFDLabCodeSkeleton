#include "uvp.h"
#include "init.h"

#include <math.h>
#include <stdio.h>

void calculate_fg(
double Re,
double GX,
double GY,
double alpha,
double dt,
double dx,
double dy,
int imax,
int jmax,
double **U,
double **V,
double **F,
double **G
){
    double a, b,c,d,du2x2,du2y2,du2dx,duvy,dv2y2,dv2x2,dv2dy,duvx;
    for (int j=1; j<=jmax;j++){
        F[0][j]=U[0][j];
        F[imax][j]=U[imax][j];
    }
    for (int i=1; i<=imax; i++){
        G[i][0]=V[i][0];
        G[i][jmax]=V[i][jmax];
    }
    for (int i=1;i<imax;i++){
        for (int j=1;j<=jmax;j++){
            du2x2= (U[i+1][j]-2*U[i][j]+U[i-1][j])/(dx*dx);
            du2y2= (U[i][j+1]-2*U[i][j]+U[i][j-1])/(dy*dy);
            a=(U[i][j]+U[i+1][j])/2;
            b=(U[i-1][j]+U[i][j])/2;
            du2dx=(a*a-b*b+ alpha*(fabs(a)*((U[i][j]-U[i+1][j])/2)-fabs(b)*((U[i-1][j]-U[i][j])/2)))/dx;
            duvy=((V[i][j]+V[i+1][j])*(U[i][j]+U[i][j+1])-(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]+U[i][j])+alpha*(fabs(V[i][j]+V[i+1][j])*(U[i][j]-U[i][j+1])-fabs(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]-U[i][j])))/(4*dy);
            F[i][j]=U[i][j]+dt*((du2x2+du2y2)*(1/Re)-du2dx-duvy+GX);
        }
    }
    for (int i=1;i<=imax;i++){
        for (int j=1;j<jmax;j++){
            dv2y2= (V[i][j+1]-2*V[i][j]+V[i][j-1])/(dy*dy);
            dv2x2= (V[i+1][j]-2*V[i][j]+V[i-1][j])/(dx*dx);
            c=(V[i][j]+V[i][j+1])/2;
            d=(V[i][j-1]+V[i][j])/2;
            dv2dy=(c*c-d*d+ alpha*(fabs(c)*((V[i][j]-V[i][j+1])/2)-fabs(d)*((V[i][j-1]-V[i][j])/2)))/dy;
            duvx=((U[i][j]+U[i][j+1])*(V[i][j]+V[i+1][j])-(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]+V[i][j])+alpha*(fabs(U[i][j]+U[i][j+1])*(V[i][j]-V[i+1][j])-fabs(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]-V[i][j])))/(4*dy);
            G[i][j]=V[i][j]+dt*((dv2x2+dv2y2)*(1/Re)-dv2dy-duvx+GY);
        }
    }
   
    return;
}

void calculate_rs(
double dt,
double dx,
double dy,
int imax,
int jmax,
double **F,
double **G,
double **RS
){
    int i,j;
    for (i=1;i<=imax;i++){
        for (j=1;j<=jmax;j++){
           RS[i][j] = ((F[i][j] - F[i-1][j])/dx + (G[i][j] - G[i][j-1])/dy)/dt;
        }
    }
    return;
}

void calculate_dt(
double Re,
double tau,
double *dt,
double dx,
double dy,
int imax,
int jmax,
double **U,
double **V
){
    double dt1,dt2,dt3;
    double U1=fabs(U[0][0]);
    double V1=fabs(V[0][0]);
    
    for(int c=0 ; c <=imax ; c++ ){
        for(int d= 0 ; d <=jmax ; d++ ){
            if ( fabs(U[c][d]) > fabs(U1) )
                U1= U[c][d];
            if ( fabs(V[c][d]) > fabs(V1) )
                V1 = V[c][d];
        }
    }
    dt1 = 0.5*Re/(1/(dx*dx) + 1/(dy*dy));
    dt2 = dx/fabs(U1);
    dt3 = dy/fabs(V1);
    
    *dt = dt1;
    if (dt2 < dt1){
        *dt =dt2;
        if (dt3 < dt2){
            *dt = dt3;
        }
    }
    else if (dt3 < dt1){
        *dt = dt3;
    }
    *dt = *dt*tau;
    return;
}

void calculate_uv(
double dt,
double dx,
double dy,
int imax,
int jmax,
double **U,
double **V,
double **F,
double **G,
double **P){
    for (int i=1;i<imax;i++){
        for (int j=1;j<=jmax;j++){
            U[i][j] = F[i][j] - dt*(P[i+1][j] - P[i][j])/dx;
        }
    }
    for (int i=1;i<=imax;i++){
        for (int j=1;j<jmax;j++){
            V[i][j] = G[i][j] - dt*(P[i][j+1] - P[i][j])/dy;
        }
    }
    return;
}
