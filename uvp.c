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
                    double alpha,
                    double **T,
                    double TI,
                    double T_h,
                    double T_c,
                    const char* problem
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
        for(int j=0; j<jmax; j++)
        {
            T[0][j] = 2*TH - T[1][j];
            T[imax-1][j] = 2*TC - T[imax-2][j];
        }
    }
    
    if( strcmp(problem,"fluid_trap") )
    {
        for(int j=0; j<jmax; j++)
        {
            T[0][j] = 2*TH - T[1][j];
            T[imax-1][j] = 2*TC - T[imax-2][j];
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
                  double beta,
                  int **flag,
                  double **T,
                  int include_T
                  ){
    double d2u_dx2, d2u_dy2, d2v_dx2, d2v_dy2;
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
                    d2v_dy2 = (U[i][j-1]-2*U[i][j]+U[i][j+1])/(dy*dy);
                    
                    F[i][j] = U[i][j] + dt*(((d2u_dx2 + d2u_dy2)/Re) - 0.25*(1/dx)*((pow((U[i+1][j]+U[i][j]),2.0) - pow((U[i-1][j]+U[i][j]),2.0))
                                                                                    +alpha*(fabs(U[i+1][j]+U[i][j])*(U[i][j]-U[i+1][j])-fabs(U[i-1][j]+U[i][j])*(U[i-1][j]-U[i][j]))
                                                                                    )/dx -0.25*(1/dy)*(((V[i][j]+V[i+1][j])*(U[i][j]+U[i][j+1])- (V[i][j-1]+V[i+1][j-1])*(U[i][j-1]+U[i][j]))+alpha*(fabs(V[i][j]+V[i+1][j])*(U[i][j]-U[i][j+1])-fabs(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]-U[i][j])))+GX*0.5*((beta*dt)*(T[i][j]+T[i+1][j])));;
                    
                }
                
                else{
                    F[i][j] = U[i][j] + dt*(((d2u_dx2 + d2u_dy2)/Re) - 0.25*(1/dx)*((pow((U[i+1][j]+U[i][j]),2.0) - pow((U[i-1][j]+U[i][j]),2.0))
                                                                                    +alpha*(fabs(U[i+1][j]+U[i][j])*(U[i][j]-U[i+1][j])-fabs(U[i-1][j]+U[i][j])*(U[i-1][j]-U[i][j]))
                                                                                    )/dx -0.25*(1/dy)*(((V[i][j]+V[i+1][j])*(U[i][j]+U[i][j+1])- (V[i][j-1]+V[i+1][j-1])*(U[i][j-1]+U[i][j]))+alpha*(fabs(V[i][j]+V[i+1][j])*(U[i][j]-U[i][j+1])-fabs(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]-U[i][j])))+GX);
                }
            }
        }
    }
    for (int i=1;i<imax;i++){
        for (int j=1;j<jmax-1;j++){
            if((flag[i][j]&(1<<0))&flag[i+1][j])
            {
                if(include_T){
            G[i][j] = V[i][j] + dt*(((d2v_dx2 + d2v_dy2)/Re) -(1/dx)*0.25*(
                                                                           ((U[i][j]+U[i][j+1])*(V[i][j]+V[i+1][j])- (U[i-1][j]+U[i-1][j+1])*(V[i-1][j]+V[i][j]))
                                                                           +alpha*(fabs(U[i][j]+U[i][j+1])*(V[i][j]-V[i+1][j])-fabs(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]-V[i][j]))
                                                                           )
                                    
                                    -(1/dy)*0.25*(
                                                  (pow((V[i][j]+V[i][j+1]),2.0) - pow((V[i][j-1]+V[i][j]),2.0))
                                                  +alpha*(fabs(V[i][j]+V[i][j+1])*(V[i][j]-V[i][j+1])-fabs(V[i][j-1]+V[i][j])*(V[i][j-1]-V[i][j]))
                                                  )
                                    +GY*0.5*(((beta*dt)*(T[i][j]+T[i][j+1]))));;
        }
                else {
                    G[i][j] = V[i][j] + dt*(((d2v_dx2 + d2v_dy2)/Re) -(1/dx)*0.25*(
                                                                                   ((U[i][j]+U[i][j+1])*(V[i][j]+V[i+1][j])- (U[i-1][j]+U[i-1][j+1])*(V[i-1][j]+V[i][j]))
                                                                                   +alpha*(fabs(U[i][j]+U[i][j+1])*(V[i][j]-V[i+1][j])-fabs(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]-V[i][j]))
                                                                                   )
                                            
                                            -(1/dy)*0.25*(
                                                          (pow((V[i][j]+V[i][j+1]),2.0) - pow((V[i][j-1]+V[i][j]),2.0))
                                                          +alpha*(fabs(V[i][j]+V[i][j+1])*(V[i][j]-V[i][j+1])-fabs(V[i][j-1]+V[i][j])*(V[i][j-1]-V[i][j]))
                                                          )
                                            +GY);;
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
                  double **RS
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
                  int include_T
                  ){
    double dt1,dt2,dt3,dt4;
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
    
    *dt = fmin(dt1,dt2);
    *dt = fmin(dt,dt3);
    
    if(include_T){
    dt4 = 0.5*Re*Pr/(1/(dx*dx) + 1/(dy*dy));
        *dt = fmin(dt4,dt);
    }
    if( (tau > 0) && (tau < 1))
    {
        *dt = *dt*tau;
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
