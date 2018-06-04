#include "boundary_val.h"
#include <stdio.h>
#include <mpi.h>

void boundaryvalues(int imax, int jmax, double** U, double** V, int iMaxUF,
		int jMaxUF, int iMaxVG, int jMaxVG, int rank_l, int rank_r, int rank_b,
		int rank_t) {
	iMaxUF = iMaxUF - 1;
	jMaxUF = jMaxUF - 1;
	iMaxVG = iMaxVG - 1;
	jMaxVG = jMaxVG - 1;
	printf("Got here, boundaryvalues.");
	//Left Boundary
	if ( MPI_PROC_NULL == rank_l) {
		for (int ctrj = 0; ctrj < jMaxUF + 1; ctrj++) {
			U[0][ctrj] = 0;
			U[1][ctrj] = 0;
		}
		for (int ctrj = 0; ctrj < jMaxVG + 1; ctrj++) {
			V[0][ctrj] = -V[1][ctrj];
		}
	}

	//Right Boundary
	if ( MPI_PROC_NULL == rank_r) {
		for (int ctrj = 0; ctrj < jMaxUF + 1; ctrj++) {
			U[iMaxUF][ctrj] = 0;
		}
		for (int ctrj = 0; ctrj < jMaxVG + 1; ctrj++) {
			V[iMaxVG][ctrj] = -V[iMaxVG - 1][ctrj];
		}
	}

	//Bottom Boundary
	if ( MPI_PROC_NULL == rank_b) {
		for (int ctri = 0; ctri < iMaxVG + 1; ctri++) {
			V[ctri][0] = 0;
			V[ctri][1] = 0;
		}
		for (int ctri = 0; ctri < iMaxUF + 1; ctri++) {
			U[ctri][0] = -U[ctri][1];
		}
	}

	//Top Boundary
	if ( MPI_PROC_NULL == rank_t) {
		for (int ctri = 0; ctri < iMaxVG + 1; ctri++) {
			V[ctri][jMaxVG] = 0;
		}
		for (int ctri = 0; ctri < iMaxUF + 1; ctri++) {
			U[ctri][jMaxUF] = 2 - U[ctri][jMaxUF - 1];
		}
	}
}

