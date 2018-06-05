#include "sor.h"
#include <math.h>
#include <mpi.h>
#include "parallel.h"

void sor(double omg, double dx, double dy, int imax, int jmax, double **P,
		double **RS, double *res, int il, int ir, int jb, int jt, int iTotElsUF,
		int jTotElsUF, int iTotElsVG, int jTotElsVG, int rank_l, int rank_r, int rank_b,
		int rank_t, double *bufSend, double *bufRecv, int chunk) {
	int iRS, jRS;
	//int iMaxU = iMaxUF - 1;
	int iMaxPlus1P = iTotElsVG - 1;
	int jMaxPlus1P = jTotElsUF - 1;
	//int jMaxV = jMaxVG - 1;
	double rloc;
	MPI_Status status;
	int locRank;
	double coeff = omg / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy)));

	/* set pressure boundary values */

	if (rank_l == MPI_PROC_NULL) {
		for (int j = 1; j < jMaxPlus1P; j++) {
			P[0][j] = P[1][j];
		}
	}
	if (rank_r == MPI_PROC_NULL) {
		for (int j = 1; j < jMaxPlus1P; j++) {
			P[iMaxPlus1P][j] = P[iMaxPlus1P - 1][j];
		}
	}
	if (rank_b == MPI_PROC_NULL) {
		for (int i = 1; i < iMaxPlus1P; i++) {
			P[i][0] = P[i][1];
		}
	}
	if (rank_t == MPI_PROC_NULL) {
		for (int i = 1; i < iMaxPlus1P; i++) {
			P[i][jMaxPlus1P] = P[i][jMaxPlus1P - 1];
		}
	}

	MPI_Comm_rank(MPI_COMM_WORLD, &locRank);

	/* SOR iteration */
	for (int i = 1; i < iMaxPlus1P; i++) {
		for (int j = 1; j < jMaxPlus1P; j++) {
			iRS = i - 1;
			jRS = j - 1;
			P[i][j] = (1.0 - omg) * P[i][j]	+ coeff	* ((P[i + 1][j] + P[i - 1][j]) / (dx * dx)	+ (P[i][j + 1] + P[i][j - 1]) / (dy * dy)- RS[iRS][jRS]);
		}
	}

// Exchange Pressures
	pressure_comm(P, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, bufSend,
			bufRecv, &status, chunk);


	/* compute the residual */
	rloc = 0;
	for (int i = 1; i < iMaxPlus1P; i++) {
		for (int j = 1; j < jMaxPlus1P; j++) {
			iRS = i - 1;
			jRS = j - 1;
			rloc += ((P[i + 1][j] - 2.0 * P[i][j] + P[i - 1][j]) / (dx * dx) + (P[i][j + 1] - 2.0 * P[i][j] + P[i][j - 1]) / (dy * dy)- RS[iRS][jRS]) * ((P[i + 1][j] - 2.0 * P[i][j] + P[i - 1][j]) / (dx * dx) + (P[i][j + 1] - 2.0 * P[i][j] + P[i][j - 1])/ (dy * dy) - RS[iRS][jRS]);
		}
	}

	double glRes, norm;

	/*Send local residual sum to main process and set residual*/
	MPI_Reduce(&rloc, &glRes, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (locRank == 0) {

		norm = glRes / (imax * jmax);
		norm = sqrt(norm);
	}
	MPI_Bcast(&norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	*res = norm; //Updated Total Residual
	MPI_Barrier(MPI_COMM_WORLD); /* synchronize*/
}
