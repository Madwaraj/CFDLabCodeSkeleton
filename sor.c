#include "sor.h"
#include <math.h>
#include "mpi.h"

void sor(double omg, double dx, double dy, int imax, int jmax, double **P,
         double **RS, double *res, int il, int ir, int jb, int jt, int rank_l,
         int rank_r, int rank_b, int rank_t, double *bufSend, double *bufRecv
         
         ) {
    int i, j;
    double rloc;
    MPI_Status *status;
    int chunk = 0;
    int locRank;
    double coeff = omg / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy)));
    
    MPI_Comm_rank(MPI_COMM_WORLD, &locRank);
    
    /* SOR iteration */
    for (i = il; i <= imax; i++) {
        for (j = 1; j <= jmax; j++) {
            P[i][j] = (1.0 - omg) * P[i][j]
            + coeff
            * ((P[i + 1][j] + P[i - 1][j]) / (dx * dx)
               + (P[i][j + 1] + P[i][j - 1]) / (dy * dy)
               - RS[i][j]);
        }
    }
    // Exchange Pressures
    pressure_comm(P, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, bufSend,
                  bufRecv, status, chunk);
    
    /* compute the residual */
    rloc = 0;
    for (i = il; i < ir + 1; i++) {
        for (j = jb; j < jt + 1; j++) {
            rloc += ((P[i + 1][j] - 2.0 * P[i][j] + P[i - 1][j]) / (dx * dx)
                     + (P[i][j + 1] - 2.0 * P[i][j] + P[i][j - 1]) / (dy * dy)
                     - RS[i][j])
            * ((P[i + 1][j] - 2.0 * P[i][j] + P[i - 1][j]) / (dx * dx)
               + (P[i][j + 1] - 2.0 * P[i][j] + P[i][j - 1])
               / (dy * dy) - RS[i][j]);
        }
    }
    
    double glRes;
    /*Send local residual sum ain process and set residual*/
    MPI_Reduce(&rloc, &glRes, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    double norm;
    
    if (locRank == 0) {
        norm = glRes / (imax * jmax);
        norm = sqrt(norm);
        MPI_Bcast(&norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    *res = norm;
    
    /* set boundary values */
    
    for (i = il; i < ir+2; i++) {
        P[i][0] = P[i][1];
        P[i][jt + 1] = P[i][jt];
    }
    for (j = jb; j < jt+2; j++) {
        P[0][j] = P[1][j];
        P[ir + 1][j] = P[ir][j];
    }
}

if ( MPI_PROC_NULL == rank_r )
{
    for (int j=jb-1; j<=jt+1;j++){
        P[ir+1][j]=P[ir][j];
    }
}
if ( MPI_PROC_NULL == rank_l ){
    for (int j=jb-1; j<=jt+1;j++){
        P[0][j] = P[1][j];
    }
}
if ( MPI_PROC_NULL == rank_b)
{
    for (int i=il-1; i<=ir+1; i++){
        P[i][0]=P[i][1];
        //     G[i][jmax]=V[i][jmax];
    }
}
if ( MPI_PROC_NULL == rank_t )
{
    for (int i=il-1; i<=ir+1; i++){
        //   G[i][0]=V[i][0];
        P[i][jt+1]=P[i][jt];
    }
}
