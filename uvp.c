#include <math.h>
#include "uvp.h"
#include"helper.h"
#include"boundary_val.h"
#include <stdio.h>

/////////////////////////////////////////////////////////////////////
void calculate_dt(double Re,
                  double tau,
                  double *dt,
                  double dx,
                  double dy,
                  int imax,
                  int jmax,
                  double **U,
                  double **V,double Pr, int include_T)
{

 double Umax = fabs(U[0][0]);
 
   for( int i = 0 ; i < imax ; i++ )
   {
      for( int j = 0 ; j < jmax ; j++ )
      {
         if ( fabs(U[i][j]) > Umax )
            Umax = fabs(U[i][j]);
      }
   }
   
  double Vmax = fabs(V[0][0]);
 
   for( int i = 0 ; i < imax ; i++ )
   {
      for( int j = 0 ; j < jmax ; j++ )
      {
         if ( fabs(V[i][j]) > Vmax )
            Vmax = fabs(V[i][j]);
      }
   }
 
 double dt1, dt2, dt3, dt4, tmin;
 
    dt1 = 0.5*Re/(1/(dx*dx) + 1/(dy*dy));

 dt2 = dx/(fabs(Umax));

 dt3 = dy/(fabs(Vmax));
 
 tmin = fmin(dt1,dt2);
 tmin = fmin(tmin,dt3);

 if(include_T){
     dt4 = (0.5*Re*Pr)/(1/(dx*dx) + 1/(dy*dy));
    tmin = fmin(tmin,dt4);
 }

 
 if( (tau > 0) && (tau < 1))
 {
   *dt = tau * tmin;
 }
}

/////////////////////////////////////////////////////////////////////////////////
void calculate_fg(double Re,
		 double GX, double GY,
		 double alpha,
		 double dt,
		 double dx, double dy,
		 int imax, int jmax,
		 double** U, double** V,
		 double** F, double** G,int **flag,
		 double beta, double** T, int include_T)
{

for(int i = 0; i<imax; ++i)
{
	for(int j = 0; j<jmax; ++j)
	{
		if ( B_O(flag[i][j]) )  F[i][j] = U[i][j];

		if ( B_W(flag[i][j]) )  F[i-1][j] = U[i-1][j];

		if ( B_N(flag[i][j]) )  G[i][j] = V[i][j];

		if ( B_S(flag[i][j]) )  G[i][j-1] = V[i][j-1];

		if ( B_NO(flag[i][j]) ) { F[i][j] = U[i][j]; G[i][j] = V[i][j]; }

		if ( B_NW(flag[i][j]) ) { F[i-1][j] = U[i-1][j]; G[i][j] = V[i][j]; }

		if ( B_SO(flag[i][j]) ) { F[i][j] = U[i][j]; G[i][j-1] = V[i][j-1]; }

		if ( B_SW(flag[i][j]) ) { F[i-1][j] = U[i-1][j]; G[i][j-1] = V[i][j-1]; }

		if (flag[i][j]&(1<<4) ) F[i][j] = U[i][j];
	}
}
    double a, b,c,d,du2x2,du2y2,du2dx,duvy,dv2y2,dv2x2,dv2dy,duvx;


    for(int i=0; i<imax-1; i++)
    {
        for(int j=0; j<jmax; j++)
	{
	if( ((flag[i][j]&(1<<0))&flag[i+1][j]) || ( (flag[i+1][j] & (1<<3)) && (flag[i][j]&(1<<0))) )
	{
        if(include_T)
        {
            du2x2= (U[i+1][j]-2*U[i][j]+U[i-1][j])/(dx*dx);
            du2y2= (U[i][j+1]-2*U[i][j]+U[i][j-1])/(dy*dy);
            a=(U[i][j]+U[i+1][j])/2;
            b=(U[i-1][j]+U[i][j])/2;
            du2dx=(a*a-b*b+ alpha*(fabs(a)*((U[i][j]-U[i+1][j])/2)-fabs(b)*((U[i-1][j]-U[i][j])/2)))/dx;
            duvy=((V[i][j]+V[i+1][j])*(U[i][j]+U[i][j+1])-(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]+U[i][j])+alpha*(fabs(V[i][j]+V[i+1][j])*(U[i][j]-U[i][j+1])-fabs(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]-U[i][j])))/(4*dy);
            F[i][j]=U[i][j]+dt*((du2x2+du2y2)*(1/Re)-du2dx-duvy+GX*((beta*dt)*(T[i][j]+T[i+1][j]))/2);
            
        }
        else
        {
            du2x2= (U[i+1][j]-2*U[i][j]+U[i-1][j])/(dx*dx);
            du2y2= (U[i][j+1]-2*U[i][j]+U[i][j-1])/(dy*dy);
            a=(U[i][j]+U[i+1][j])/2;
            b=(U[i-1][j]+U[i][j])/2;
            du2dx=(a*a-b*b+ alpha*(fabs(a)*((U[i][j]-U[i+1][j])/2)-fabs(b)*((U[i-1][j]-U[i][j])/2)))/dx;
            duvy=((V[i][j]+V[i+1][j])*(U[i][j]+U[i][j+1])-(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]+U[i][j])+alpha*(fabs(V[i][j]+V[i+1][j])*(U[i][j]-U[i][j+1])-fabs(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]-U[i][j])))/(4*dy);
            F[i][j]=U[i][j]+dt*((du2x2+du2y2)*(1/Re)-du2dx-duvy+GX);
	}

	}
	}
    }


    for(int i=0; i<imax; i++)
	{
        for(int j=0; j<jmax-1; j++)
	{
	if((flag[i][j]&(1<<0))&flag[i][j+1])
	{
        if(include_T)
        {
            dv2y2= (V[i][j+1]-2*V[i][j]+V[i][j-1])/(dy*dy);
            dv2x2= (V[i+1][j]-2*V[i][j]+V[i-1][j])/(dx*dx);
            c=(V[i][j]+V[i][j+1])/2;
            d=(V[i][j-1]+V[i][j])/2;
            dv2dy=(c*c-d*d+ alpha*(fabs(c)*((V[i][j]-V[i][j+1])/2)-fabs(d)*((V[i][j-1]-V[i][j])/2)))/dy;
            duvx=((U[i][j]+U[i][j+1])*(V[i][j]+V[i+1][j])-(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]+V[i][j])+alpha*(fabs(U[i][j]+U[i][j+1])*(V[i][j]-V[i+1][j])-fabs(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]-V[i][j])))/(4*dy);
            G[i][j]=V[i][j]+dt*((dv2x2+dv2y2)*(1/Re)-dv2dy-duvx+GY*(((beta*dt)*(T[i][j]+T[i][j+1]))/2));
        }
        else
        {
            dv2y2= (V[i][j+1]-2*V[i][j]+V[i][j-1])/(dy*dy);
            dv2x2= (V[i+1][j]-2*V[i][j]+V[i-1][j])/(dx*dx);
            c=(V[i][j]+V[i][j+1])/2;
            d=(V[i][j-1]+V[i][j])/2;
            dv2dy=(c*c-d*d+ alpha*(fabs(c)*((V[i][j]-V[i][j+1])/2)-fabs(d)*((V[i][j-1]-V[i][j])/2)))/dy;
            duvx=((U[i][j]+U[i][j+1])*(V[i][j]+V[i+1][j])-(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]+V[i][j])+alpha*(fabs(U[i][j]+U[i][j+1])*(V[i][j]-V[i+1][j])-fabs(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]-V[i][j])))/(4*dy);
            G[i][j]=V[i][j]+dt*((dv2x2+dv2y2)*(1/Re)-dv2dy-duvx+GY);
	
	}
	}
        }

    }

}


////////////////////////////////////////////////////////////////////////////////
void calculate_uv(double dt,double dx,double dy,int imax, int jmax,
		 double**U, double**V,double**F,double**G,double **P,int **flag)
{
	for (int i = 0; i< imax-1;i++)
	{
		for (int j = 0; j<jmax;j++)
		{
			if(((flag[i][j]&(1<<0))&flag[i+1][j]) || ( (flag[i+1][j] & (1<<3)) && (flag[i][j]&(1<<0))))
			U[i][j] = F[i][j] - (dt/dx)*(P[i+1][j]-P[i][j]);
		}
	}
	
	for (int i = 0; i< imax;i++)
	{
		for (int j = 0; j<jmax-1;j++)
		{
			if((flag[i][j]&(1<<0))&flag[i][j+1])
			V[i][j] = G[i][j] - (dt/dy)*(P[i][j+1] -P[i][j]);
		}
	}

}

/////////////////////////////////////////////////////////////////
void calculate_rs(double dt,
		  double dx,
		  double dy,
		  int imax,
		  int jmax,
		  double **F,
		  double **G,
		  double **RS,int **flag)
{

	for(int i=0; i<imax; i++)
   	{
        	for(int j=0; j<jmax; j++)
		{
			if(flag[i][j]&(1<<0))
			RS[i][j] =  (1/dt)*( (F[i][j]-F[i-1][j])/dx + (G[i][j]-G[i][j-1])/dy );
		}
	}
}

void calculate_temp(double **T, double **T1, double Pr, double Re, int imax,int jmax,double dx, double dy,
		double dt, double alpha,double **U,double **V,int **flag, double TI, double T_h, double T_c, int select)
{   
    for(int i = 0; i<imax; ++i)
    {
  	for(int j = 0; j<jmax; ++j)
  	{
  		if ( B_O(flag[i][j]) )  T[i][j] = T[i+1][j];

  		if ( B_W(flag[i][j]) )  T[i][j] = T[i-1][j];

  		if ( B_N(flag[i][j]) )  T[i][j] = T[i][j+1];

  		if ( B_S(flag[i][j]) )  T[i][j] = T[i][j-1];

  		if ( B_NO(flag[i][j]) ) T[i][j] = (T[i][j+1] + T[i+1][j])/2;

  		if ( B_NW(flag[i][j]) ) T[i][j] = (T[i][j+1] + T[i-1][j])/2;

  		if ( B_SO(flag[i][j]) ) T[i][j] = (T[i][j-1] + T[i+1][j])/2;

  		if ( B_SW(flag[i][j]) ) T[i][j] = (T[i][j-1] + T[i-1][j])/2;

  		if (flag[i][j]&(1<<3) ) T[i][j] = T[i-1][j];

  		if (flag[i][j]&(1<<4) ) T[i][j] = TI;
  	}
  }

switch(select)
{
	case 3:
	for(int j=0; j<jmax; j++)
	{
		T[0][j] = 2*T_h - T[1][j];
		T[imax-1][j] = 2*T_c - T[imax-2][j];
	}
	break;
	
	case 4:
	for(int j=0; j<jmax; j++)
	{
		T[0][j] = 2*T_h - T[1][j];
		T[imax-1][j] = 2*T_c - T[imax-2][j];
	}
	break;

	case 5:
	for(int i=0; i<imax; i++)
	{
		T[i][0] = 2*T_h - T[i][1];
		T[i][jmax-1] = 2*T_c - T[i][jmax-2];
	}

}

double dut_dx;
double dvt_dy;
double dt2_dx2;
double dt2_dy2;
double Z;
for (int i = 0; i< imax;i++)
{
    for (int j=0; j<jmax;j++)
	{

	if(flag[i][j]&((1<<0)|(1<<3)|(1<<4)))
	{

  		dut_dx = (1/dx)*( (U[i][j]*(T[i][j]+T[i+1][j])*0.5)-(U[i-1][j]*(T[i-1][j]+T[i][j])*0.5))+(alpha/dx)*(
				(fabs(U[i][j])*(T[i][j]-T[i+1][j])*0.5) - (fabs(U[i-1][j])*(T[i-1][j]-T[i][j])*0.5)
				);

  		dvt_dy = (1/dy)*( (V[i][j]*(T[i][j]+T[i][j+1])*0.5)-(V[i][j-1]*(T[i][j-1]+T[i][j])*0.5) )+(alpha/dy)*(
				(fabs(V[i][j])*(T[i][j]-T[i][j+1])*0.5) - (fabs(V[i][j-1])*(T[i][j-1]-T[i][j])*0.5)
				);

  		dt2_dx2 = (T[i+1][j] - 2*T[i][j] + T[i-1][j])/(dx*dx);

  		dt2_dy2 = (T[i][j+1] - 2*T[i][j] +T[i][j-1])/(dy*dy);

  		Z = (1/(Re*Pr))*(dt2_dx2+dt2_dy2) - dut_dx - dvt_dy;

      	T1[i][j] = T[i][j]+ (dt*Z);
    }
  }
	
}

  for (int i = 0; i< imax;i++){
    for (int j=0;j<jmax;j++){

		if(flag[i][j]&((1<<0)|(1<<3)|(1<<4))){

		T[i][j] = T1[i][j];
    }
  }
}

}

///////////////////////////////////////////////////////////////////////////////

void normal_boundary(double **U, double **V, double **P, int **flag, int imax, int jmax)
{

	for (int i = 0; i< imax;i++)
	{
	    for (int j=0;j<jmax; j++)
		{
			if(flag[i][j]&( (1<<1)|(1<<2)) )
			{
				U[i][j] = 0;
				V[i][j] = 0;
				P[i][j] = 0;
			}
			
		}
	}

}


void normal_boundary_T(double **U, double **V, double **P, double **T, int **flag, int imax, int jmax)
{
	for (int i = 0; i< imax;i++)
	{
	    for (int j=0;j<jmax; j++)
		{
			if(flag[i][j]&( (1<<1)|(1<<2)) )
			{
				U[i][j] = 0;
				V[i][j] = 0;
				P[i][j] = 0;
				T[i][j] = 0;
			}
		}
	}


}

