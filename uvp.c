#include "uvp.h"
#include "init.h"

#include <math.h>
#include <stdio.h>

void calculate_temp(
double dt,
double **U,
double **V,
double dx,
double dy,
double Re,
double Pr,
double gama,
double **T
){
    double d2t_dx2, d2t_dy2, dut_dx, dvt_dy;
    for (int i=1; i<imax; i++){
        for (int j=1; j<=jmax; j++){
            d2t_dx2[i][j] = (T[i+1][j] - 2*T[i][j] + T[i-1][j])/(dx*dx);
            dut_dx[i][j] = 0.5*dx*(U[i][j]*(T[i][j]+T[i+1][j]) - U[i-1][j]*(T[i-1][j]+T[i][j])) + 0.5*gama*(U[i][j]*(T[i][j]-T[i+1][j]) - U[i-1][j]*(T[i-1][j]-T[i][j]))/dx;
       }
    }
    for (int i=1; i<=imax; i++){
        for (int j=1; j<jmax; j++){
            dvt_dy[i][j] = 0.5*dy*(V[i][j]*(T[i][j]+T[i][j+1]) - V[i-1][j]*(T[i][j-1]+T[i][j])) + 0.5*gama*(V[i][j]*(T[i][j]-T[i][j+1]) - V[i][j-1]*(T[i][j-1]-T[i][j]))/dx;
            d2t_dy2[i][j] = (T[i][j+1] - 2*T[i][j] + T[i][j-1])/(dy*dy);
        }
    }
    for (int i=1;i<imax;i++){
        for (int j=1;j<jmax;j++){
            T[i][j] = T[i][j] + dt*(-dut_dx[i][j] - dvt_dy[i][j] + (d2t_dy2[i][j] + d2t_dx2[i][j])/(Re*Pr));
        }
    }
}

void calculate_fg(
double Re,
double GX,
double GY,
double beta,
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
 
    for (int i=1;i<imax;i++){
        for (int j=1;j<=jmax;j++){
         F[i][j] = F[i][j] - 0.5*beta*dt*(T[i+1][j]+T[i][j]);
        }
    }
    for (int i=1;i<=imax;i++){
        for (int j=1;j<jmax;j++){
          G[i][j] = G[i][j] - 0.5*beta*dt*(T[i][j+1]+T[i][j]);
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
    dt1 = 0.5*Re*Pr/(1/(dx*dx) + 1/(dy*dy));
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
