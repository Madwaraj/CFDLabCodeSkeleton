#include "boundary_val.h"
#include <stdio.h>
#include <mpi.h>
#include "parallel.h"

void boundaryvalues(int imax, int jmax, double** U, double** V, int iTotElsUF,
		int jTotElsUF, int iTotElsVG, int jTotElsVG, int rank_l, int rank_r, int rank_b,
		int rank_t) {
	int iMaxPlus1U = iTotElsUF - 1;
	int jMaxPlus1U = jTotElsUF - 1;
	int iMaxPlus1V = iTotElsVG - 1;
	int jMaxPlus1V = jTotElsVG - 1;

	//Left Boundary
	if (rank_l == MPI_PROC_NULL) {
		for (int ctrj = 1; ctrj < jMaxPlus1U; ctrj++) {
			U[0][ctrj] = 0;
			U[1][ctrj] = 0;
		}
		for (int ctrj = 1; ctrj < jMaxPlus1V; ctrj++) {
			V[0][ctrj] = -V[1][ctrj];
		}
	}

	//Right Boundary
	if (rank_r == MPI_PROC_NULL) {
		for (int ctrj = 1; ctrj < jMaxPlus1U; ctrj++) {
			U[iMaxPlus1U-1][ctrj] = 0;
		}
		for (int ctrj = 1; ctrj < jMaxPlus1V; ctrj++) {
			V[iMaxPlus1V][ctrj] = -V[iMaxPlus1V - 1][ctrj];
		}
	}

	//Bottom Boundary
	if (rank_b == MPI_PROC_NULL) {
		for (int ctri = 1; ctri < iMaxPlus1V; ctri++) {
			V[ctri][0] = 0;
			V[ctri][1] = 0;
		}
		for (int ctri = 1; ctri < iMaxPlus1U; ctri++) {
			U[ctri][0] = -U[ctri][1];
		}
	}

	//Top Boundary
	if (rank_t == MPI_PROC_NULL) {
		for (int ctri = 1; ctri < iMaxPlus1V; ctri++) {
			V[ctri][jMaxPlus1V-1] = 0;
		}
		for (int ctri = 1; ctri < iMaxPlus1U; ctri++) {
			U[ctri][jMaxPlus1U] = 2 - U[ctri][jMaxPlus1U - 1];
		}
	}
}

