#include "boundary_val.h"
#include <stdio.h>
#include <mpi.h>
#include "parallel.h"

void boundaryvalues(int imax, int jmax, double** U, double** V, int iMaxUF,
		int jMaxUF, int iMaxVG, int jMaxVG, int rank_l, int rank_r, int rank_b,
		int rank_t) {
	int iMaxU = iMaxUF - 1;
	int jMaxU = jMaxUF - 1;
	int iMaxV = iMaxVG - 1;
	int jMaxV = jMaxVG - 1;
	//printf("iMaxU=%d, jMaxU=%d, iMaxV=%d, jMaxV=%d",iMaxU,jMaxU,iMaxV,jMaxV);
	//Programm_Sync("Got here, boundaryvalues.");
	//Left Boundary
	if ( MPI_PROC_NULL == rank_l) {
		for (int ctrj = 0; ctrj < jMaxU + 1; ctrj++) {
			U[0][ctrj] = 0;
			U[1][ctrj] = 0;
		}
		for (int ctrj = 0; ctrj < jMaxV + 1; ctrj++) {
			V[0][ctrj] = -V[1][ctrj];
		}
	}

	//Right Boundary
	if ( MPI_PROC_NULL == rank_r) {
		for (int ctrj = 0; ctrj < jMaxU + 1; ctrj++) {
			U[iMaxU][ctrj] = 0;
		}
		for (int ctrj = 0; ctrj < jMaxV + 1; ctrj++) {
			V[iMaxV][ctrj] = -V[iMaxV - 1][ctrj];
		}
	}

	//Bottom Boundary
	if ( MPI_PROC_NULL == rank_b) {
		for (int ctri = 0; ctri < iMaxV + 1; ctri++) {
			V[ctri][0] = 0;
			V[ctri][1] = 0;
		}
		for (int ctri = 0; ctri < iMaxU + 1; ctri++) {
			U[ctri][0] = -U[ctri][1];
		}
	}

	//Top Boundary
	if ( MPI_PROC_NULL == rank_t) {
		for (int ctri = 0; ctri < iMaxV + 1; ctri++) {
			V[ctri][jMaxV] = 0;
		}
		for (int ctri = 0; ctri < iMaxU + 1; ctri++) {
			U[ctri][jMaxU] = 2 - U[ctri][jMaxU - 1];
		}
	}
}

