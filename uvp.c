#include "uvp.h"
#include "init.h"

#include <stdio.h>

void calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, U, V, F, G){
    for (i=1;i<imax;i++){
        for (j=1;j<imax;j++){
    F(i,j) = u(i,j) + dt*(((u(i+1,j)-2*u(i,j)+u(i-1,j))/(dx*dx)+(u(i,j+1)-2*u(i,j)+u(i,j-1))/(dy*dy))/Re - 0.25*(((u(i,j)+u(i+1,j))^2-(u(i,j)+u(i-1,j))^2)/dx + ((u(i,j)+u(i,j+1))*(v(i,j)+v(i,j+1))-(u(i,j)+u(i,j-1))*(v(i,j)*v(i,j-1)))/dy )+ GX);
    
    G(i,j) = v(i,j) + dt*(((v(i+1,j)-2*v(i,j)+v(i-1,j))/(dx*dx)+(v(i,j+1)-2*v(i,j)+v(i,j-1))/(dy*dy))/Re - 0.25*(((v(i,j)+v(i+1,j))^2-(v(i,j)+v(i-1,j))^2)/dy + ((u(i,j)+u(i,j+1))*(v(i,j)+v(i,j+1))-(u(i,j)+u(i,j-1))*(v(i,j)*v(i,j-1)))/dx )+ GY);
        }
    }
    for (j=1; j<jmax;j++){
        F(0,j) = U(0,j);
        F(imax,j) = U(imax,j);
    }
    for (i=1; i<imax; i++){
        G(i,0) = V(i,0);
        G(i,jmax) = V(i,jmax);
    }
    return;
}

void calculate_rs(dt,dx, dy, imax, jmax, F, G, RS){
    for (i=1;i<imax;i++){
        for (j=1;j,jmax;j++){
            RS = (((F(i,j)-F(i-1,j))/dx + ((G(i,j)-G(i,j-1))/dy))/dt;
        }
    }
    return;
}

void calculate_dt(){
    dt1 = 0.5*Re/(1/(dx*dx) + 1/(dy*dy));
    dt2 = dx/abs(umax);
    dt3 = dy/abs(vmax);
    
    dt = dt1;
    if (dt2 < dt1){
        dt =dt2;
        if (dt3 < dt2){
            dt = dt3;
        }
    }
    else if (dt3 < dt1){
        dt = dt3;
    }
    dt = dt*0.5;
    return dt;
}

void calculate_uv(dt, dx, dy, imax, jmax, F, G, RS){
    for (i=1;i<imax;i++){
        for (j=1;j<jmax;j++){
            U(i,j) = F(i,j) - dt*(p(i+1,j)-p(i,j))/dx;
            V(i,j) = G(i,j) - dt*(p(i,j+1)-p(i,j))/dy;
        }
    }
    
    for (i=1;i<imax;i++){
        U(i,0) = -U(i,1);
        U(i,jmax+1) = -U(i,jmax);
    }
    for (j=1;j<jmax;j++){
        V(0,j) = -V(1,j);
        V(imax+1,j) = -V(imax,j);
    }
    return;
}
