#include "boundary_val.h"
#include <stdio.h>
#include "mpi.h"

void boundaryvalues(int imax,int jmax,double** U,double** V, int il, int ir, int jb, int jt, int rank_l, int rank_r, int rank_b, int rank_t){
    if ( MPI_PROC_NULL == rank_r )
    {
        for (int j=jb-1; j<=jt+1;j++){
  //  for (int j=1;j<=jmax;j++){
     //   U[0][j] =0;
        U[imax][j]=0;
       // V[0][j] = -V[1][j];
        V[imax+1][j] = -V[imax][j];
        
    }
    }
    if ( MPI_PROC_NULL == rank_l )
    {
        for (int j=jb-1; j<=jt+1;j++){
            U[imax][j]=0;
            V[imax+1][j] = -V[imax][j];
        }
    }
    
    //for (int i=1;i<=imax;i++){
    if ( MPI_PROC_NULL == rank_b )
    {
        for (int i=il-1; i<=ir+1; i++){
            
        V[i][0]=0;
     //   V[i][jmax]=0;
        U[i][0]=-U[i][1];
      //  U[i][jmax+1] = 2-U[i][jmax];
        }
    }
    if ( MPI_PROC_NULL == rank_t )
    {
        for (int i=il-1; i<=ir+1; i++){
            U[i][0]=-U[i][1];
            V[i][0]=0;
       
    }
    }
}
