#include "sor.h"
#include <math.h>
#include <mpi.h>
#include "parallel.h"

void sor(double omg, double dx, double dy, int imax, int jmax, double **P,
		double **RS, double *res, int il, int ir, int jb, int jt, int iMaxUF,
		int jMaxUF, int iMaxVG, int jMaxVG, int rank_l, int rank_r, int rank_b,
		int rank_t, double *bufSend, double *bufRecv, int chunk) {
	int i, j, iRS, jRS;
	//int iMaxU = iMaxUF - 1;
	int jMaxU = jMaxUF - 1;
	int iMaxV = iMaxVG - 1;
	//int jMaxV = jMaxVG - 1;
	double rloc;
	MPI_Status status;
	int locRank;
	double coeff = omg / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy)));

	//printf("Inside SOR: Variables defined \n \n");

	/* set pressure boundary values */

	if ( MPI_PROC_NULL == rank_l) {
		for (int j = 0; j < jMaxU + 1; j++) {
			P[0][j] = P[1][j];
		}
	}
	if ( MPI_PROC_NULL == rank_r) {
		for (int j = 0; j < jMaxU + 1; j++) {
			P[iMaxV][j] = P[iMaxV - 1][j];
		}
	}
	if ( MPI_PROC_NULL == rank_b) {
		for (int i = 0; i < iMaxV + 1; i++) {
			P[i][0] = P[i][1];
		}
	}
	if ( MPI_PROC_NULL == rank_t) {
		for (int i = 0; i < iMaxV + 1; i++) {
			P[i][jMaxU] = P[i][jMaxU - 1];
		}
	}

	MPI_Comm_rank(MPI_COMM_WORLD, &locRank);

	/* SOR iteration */
	for (i = 1; i < iMaxV; i++) {
		for (j = 1; j < jMaxU; j++) {
			iRS = i - 1;
			jRS = j - 1;
			P[i][j] = (1.0 - omg) * P[i][j]
					+ coeff
							* ((P[i + 1][j] + P[i - 1][j]) / (dx * dx)
									+ (P[i][j + 1] + P[i][j - 1]) / (dy * dy)
									- RS[iRS][jRS]);
		}
	}

	printf("Inside SOR: SOR Iteration done \n \n");

	// Exchange Pressures
	pressure_comm(P, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, bufSend,
			bufRecv, &status, chunk);

	printf("Inside SOR: Pressures Exchanged \n \n");

	/* compute the residual */
	rloc = 0;
	for (i = 1; i < iMaxV; i++) {
		for (j = 1; j < jMaxU; j++) {
			iRS = i - 1;
			jRS = j - 1;
			rloc += ((P[i + 1][j] - 2.0 * P[i][j] + P[i - 1][j]) / (dx * dx)
					+ (P[i][j + 1] - 2.0 * P[i][j] + P[i][j - 1]) / (dy * dy)
					- RS[iRS][jRS])
					* ((P[i + 1][j] - 2.0 * P[i][j] + P[i - 1][j]) / (dx * dx)
							+ (P[i][j + 1] - 2.0 * P[i][j] + P[i][j - 1])
									/ (dy * dy) - RS[iRS][jRS]);
		}
	}

	printf("Inside SOR: Local Residual calculated \n \n");

	double glRes;

	/*Send local residual sum ain process and set residual*/

	printf("P%d: ", locRank);
	Programm_Sync("Sync for eps reduction");
	MPI_Reduce(&rloc, &glRes, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	printf("Inside SOR: MPI_Reduce called \n \n");

	double norm;

	if (locRank == 0) {
		norm = glRes / (imax * jmax);
		norm = sqrt(norm);
		MPI_Bcast(&norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	*res = norm;
	Programm_Sync("Global eps Broadcasted");
	printf("Inside SOR: Global Residual Norm calculated \n \n");

}
