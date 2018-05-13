#include "uvp.h"
#include"helper.h"
#include"boundary_val.h"

#include <math.h>
#include <stdio.h>
#include <string.h>


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
                  double **G,
                  int **flag,
                  double beta,
                  double **T,
                  int include_T
                  ){
    double a, b,c,d,du2dx,duvy,dv2dy,duvx, d2u_dx2, d2u_dy2, d2v_dx2, d2v_dy2;
    for(int i = 0; i<imax; i++)
    {
        for(int j = 0; j<jmax; j++)
        {
            if ( B_N(flag[i][j]) )  G[i][j] = V[i][j];
            
            if ( B_S(flag[i][j]) )  G[i][j-1] = V[i][j-1];
            
            if ( B_O(flag[i][j]) )  F[i][j] = U[i][j];
            
            if ( B_W(flag[i][j]) )  F[i-1][j] = U[i-1][j];
            
            if ( B_NO(flag[i][j]) ) { F[i][j] = U[i][j]; G[i][j] = V[i][j]; }
            
            if ( B_NW(flag[i][j]) ) { F[i-1][j] = U[i-1][j]; G[i][j] = V[i][j]; }
            
            if ( B_SO(flag[i][j]) ) { F[i][j] = U[i][j]; G[i][j-1] = V[i][j-1]; }
            
            if ( B_SW(flag[i][j]) ) { F[i-1][j] = U[i-1][j]; G[i][j-1] = V[i][j-1]; }
            
            if (flag[i][j]&(1<<4) ) F[i][j] = U[i][j];
        }
    }
    
    for (int i=1;i<imax-1;i++){
        for (int j=1;j<jmax;j++){
            if((flag[i][j]&(1<<0))&flag[i+1][j])
            {
                if(include_T){
                    d2u_dx2 = (U[i-1][j]-2*U[i][j]+U[i+1][j])/(dx*dx);
                    d2u_dy2 = (U[i][j-1]-2*U[i][j]+U[i][j+1])/(dy*dy);
                    a=(U[i][j]+U[i+1][j])/2;
                    b=(U[i-1][j]+U[i][j])/2;
                    du2dx=(a*a-b*b+ alpha*(fabs(a)*((U[i][j]-U[i+1][j])/2)-fabs(b)*((U[i-1][j]-U[i][j])/2)))/dx;
                     duvy=((V[i][j]+V[i+1][j])*(U[i][j]+U[i][j+1])-(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]+U[i][j])+alpha*(fabs(V[i][j]+V[i+1][j])*(U[i][j]-U[i][j+1])-fabs(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]-U[i][j])))/(4*dy);
                    F[i][j]=U[i][j]+dt*((d2u_dx2+d2u_dy2)*(1/Re)-du2dx-duvy+GX);
                    F[i][j] = F[i][j]-GX*0.5*((beta*dt)*(T[i][j]+T[i+1][j]));
                }
                
                else{
                    d2u_dx2 = (U[i-1][j]-2*U[i][j]+U[i+1][j])/(dx*dx);
                    d2u_dy2 = (U[i][j-1]-2*U[i][j]+U[i][j+1])/(dy*dy);
                    a=(U[i][j]+U[i+1][j])/2;
                    b=(U[i-1][j]+U[i][j])/2;
                    du2dx=(a*a-b*b+ alpha*(fabs(a)*((U[i][j]-U[i+1][j])/2)-fabs(b)*((U[i-1][j]-U[i][j])/2)))/dx;
                    duvy=((V[i][j]+V[i+1][j])*(U[i][j]+U[i][j+1])-(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]+U[i][j])+alpha*(fabs(V[i][j]+V[i+1][j])*(U[i][j]-U[i][j+1])-fabs(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]-U[i][j])))/(4*dy);
                    F[i][j]=U[i][j]+dt*((d2u_dx2+d2u_dy2)*(1/Re)-du2dx-duvy+GX);

                }
            }
        }
    }
    for (int i=1;i<imax;i++){
        for (int j=1;j<jmax-1;j++){
            if((flag[i][j]&(1<<0))&flag[i+1][j])
            {
                if(include_T){
                    d2v_dy2= (V[i][j+1]-2*V[i][j]+V[i][j-1])/(dy*dy);
                    d2v_dx2= (V[i+1][j]-2*V[i][j]+V[i-1][j])/(dx*dx);
                    c=(V[i][j]+V[i][j+1])/2;
                    d=(V[i][j-1]+V[i][j])/2;
                    dv2dy=(c*c-d*d+ alpha*(fabs(c)*((V[i][j]-V[i][j+1])/2)-fabs(d)*((V[i][j-1]-V[i][j])/2)))/dy;
                    duvx=((U[i][j]+U[i][j+1])*(V[i][j]+V[i+1][j])-(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]+V[i][j])+alpha*(fabs(U[i][j]+U[i][j+1])*(V[i][j]-V[i+1][j])-fabs(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]-V[i][j])))/(4*dy);
                    G[i][j]=V[i][j]+dt*((d2v_dx2+d2v_dy2)*(1/Re)-dv2dy-duvx+GY);
                    G[i][j] = G[i][j]-GY*0.5*((beta*dt)*(T[i][j]+T[i][j+1]));
        }
                else {
                    d2v_dy2= (V[i][j+1]-2*V[i][j]+V[i][j-1])/(dy*dy);
                    d2v_dx2= (V[i+1][j]-2*V[i][j]+V[i-1][j])/(dx*dx);
                    c=(V[i][j]+V[i][j+1])/2;
                    d=(V[i][j-1]+V[i][j])/2;
                    dv2dy=(c*c-d*d+ alpha*(fabs(c)*((V[i][j]-V[i][j+1])/2)-fabs(d)*((V[i][j-1]-V[i][j])/2)))/dy;
                    duvx=((U[i][j]+U[i][j+1])*(V[i][j]+V[i+1][j])-(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]+V[i][j])+alpha*(fabs(U[i][j]+U[i][j+1])*(V[i][j]-V[i+1][j])-fabs(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]-V[i][j])))/(4*dy);
                    G[i][j]=V[i][j]+dt*((d2v_dx2+d2v_dy2)*(1/Re)-dv2dy-duvx+GY);
                }
    }
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
                  double **RS,
                  int **flag
                  ){
    int i,j;
    for (i=1;i<=imax;i++){
        for (j=1;j<=jmax;j++){
            if(flag[i][j]&(1<<0))
            {
                RS[i][j] = ((F[i][j] - F[i-1][j])/dx + (G[i][j] - G[i][j-1])/dy)/dt;
            }
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
                  double **V,
                  int include_T,
                  double Pr
                  ){
    double dt1,dt2,dt3,dt4,dtmin;
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
    
    dtmin = fmin(dt1,dt2);
    dtmin = fmin(dtmin,dt3);
    
    if(include_T){
    dt4 = 0.5*Re*Pr/(1/(dx*dx) + 1/(dy*dy));
        dtmin = fmin(dt4,dtmin);
    }
    if( (tau > 0) && (tau < 1))
    {
        *dt = dtmin*tau;
    }
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
                  double **P,
                  int **flag){
    for (int i=1;i<imax;i++){
        for (int j=1;j<=jmax;j++){
            if((flag[i][j]&(1<<0))&flag[i+1][j])
            {
                U[i][j] = F[i][j] - dt*(P[i+1][j] - P[i][j])/dx;
            }
        }
    }
    for (int i=1;i<=imax;i++){
        for (int j=1;j<jmax;j++){
            if((flag[i][j]&(1<<0))&flag[i][j+1])
            {
                V[i][j] = G[i][j] - dt*(P[i][j+1] - P[i][j])/dy;
            }
        }
    }
    return;
}

void calculate_temp(
                    double **T,
                    double Pr,
                    double Re,
                    int imax,
                    int jmax,
                    double dx,
                    double dy,
                    double dt,
                    double alpha,
                    double **U,
                    double **V,
                    int **flag,
                    double TI,
                    double TH,
                    double TC,
                    char* problem
                    ){
    for(int i = 0; i<imax; ++i)
    {
        for(int j = 0; j<jmax; ++j)
        {
            if ( B_N(flag[i][j]) )  T[i][j] = T[i][j+1];
            
            if ( B_S(flag[i][j]) )  T[i][j] = T[i][j-1];
            
            if ( B_O(flag[i][j]) )  T[i][j] = T[i+1][j];
            
            if ( B_W(flag[i][j]) )  T[i][j] = T[i-1][j];
            
            if ( B_NO(flag[i][j]) ) T[i][j] = (T[i][j+1] + T[i+1][j])/2;
            
            if ( B_NW(flag[i][j]) ) T[i][j] = (T[i][j+1] + T[i-1][j])/2;
            
            if ( B_SO(flag[i][j]) ) T[i][j] = (T[i][j-1] + T[i+1][j])/2;
            
            if ( B_SW(flag[i][j]) ) T[i][j] = (T[i][j-1] + T[i-1][j])/2;
            
            if (flag[i][j]&(1<<3) ) T[i][j] = T[i-1][j];
            
            if (flag[i][j]&(1<<4) ) T[i][j] = TI;
        }
    }
    
    if( strcmp(problem,"natural_convection") )
    {
        for(int j=1; j<jmax; j++)
        {
            T[0][j] = 2*TH - T[1][j];
            T[imax-1][j] = 2*TC - T[imax-2][j];
        }
    }
    
    if( strcmp(problem,"fluid_trap") )
    {
        for(int j=1; j<=jmax; j++)
        {
            T[0][j] = 2*TH - T[1][j];
            T[imax+1][j] = 2*TC - T[imax][j];
        }
        for(int i=1; i<=imax; i++)
        {
            T[i][0] = 2*TC - T[i][1];
            T[i][jmax+1] = 2*TC - T[i][jmax];
        }
    }
    
    if( strcmp(problem,"rb_convection") )
    {
        for(int i=0; i<imax; i++)
        {
            T[i][0] = 2*TH - T[i][1];
            T[i][jmax-1] = 2*TC - T[i][jmax-2];
        }
        
    }
    double d2t_dx2, d2t_dy2, dut_dx, dvt_dy;
    for (int i=1; i<imax; i++){
        for (int j=1; j<=jmax; j++){
            if(flag[i][j]&(1<<0)){
                d2t_dx2 = (T[i+1][j] - 2*T[i][j] + T[i-1][j])/(dx*dx);
                dut_dx = 0.5*dx*(U[i][j]*(T[i][j]+T[i+1][j]) - U[i-1][j]*(T[i-1][j]+T[i][j])) + 0.5*alpha*(U[i][j]*(T[i][j]-T[i+1][j]) - U[i-1][j]*(T[i-1][j]-T[i][j]))/dx;
                
                dvt_dy = 0.5*dy*(V[i][j]*(T[i][j]+T[i][j+1]) - V[i-1][j]*(T[i][j-1]+T[i][j])) + 0.5*alpha*(V[i][j]*(T[i][j]-T[i][j+1]) - V[i][j-1]*(T[i][j-1]-T[i][j]))/dy;
                d2t_dy2 = (T[i][j+1] - 2*T[i][j] + T[i][j-1])/(dy*dy);
                
                T[i][j] = T[i][j] + dt*(-dut_dx - dvt_dy + (d2t_dy2 + d2t_dx2)/(Re*Pr));
            }
        }
    }
}

