#include "uvp.h"
#include "init.h"
#include "parallel.h"
#include <math.h>
#include <stdio.h>
#include <mpi.h>

void calculate_fg(double Re, double GX, double GY, double alpha, double dt,
		double dx, double dy, int imax, int jmax, double **U, double **V,
		double **F, double **G, int iMaxUF, int jMaxUF, int iMaxVG, int jMaxVG,
		int rank_l, int rank_r, int rank_b, int rank_t) {
	double a, b, c, d, du2x2, du2y2, du2dx, duvy, dv2y2, dv2x2, dv2dy, duvx;
	int iMaxU = iMaxUF - 1;
	int jMaxU = jMaxUF - 1;
	int iMaxV = iMaxVG - 1;
	int jMaxV = jMaxVG - 1;
	//printf("P:%d P:%d jMaxVG=%d", rank_l, rank_t, jMaxV);
	if ( MPI_PROC_NULL == rank_l) {
		for (int j = 2; j < jMaxU; j++) {
			F[0][j] = U[1][j];
			F[1][j] = U[1][j]; //Actually used in calcs
		}
	}
	if ( MPI_PROC_NULL == rank_r) {
		for (int j = 2; j < jMaxU; j++) {
			F[iMaxU - 1][j] = U[iMaxU - 1][j];
		}
	}
	if ( MPI_PROC_NULL == rank_b) {
		for (int i = 1; i < iMaxV; i++) {
			G[i][0] = V[i][1];
			G[i][1] = V[i][1]; //Actually used in calcs
			//     G[i][jmax]=V[i][jmax];
		}
	}
	if ( MPI_PROC_NULL == rank_t) {
		for (int i = 1; i < iMaxV; i++) {
			//   G[i][0]=V[i][0];
			G[i][jMaxV - 1] = V[i][jMaxV - 1];
		}
	}

	int iVforFCalc, jVforFCalc;
	for (int i = 2; i < iMaxU - 1; i++) {
		for (int j = 1; j < jMaxU; j++) {
			iVforFCalc = i - 1;
			jVforFCalc = j + 1;
			du2x2 = (U[i + 1][j] - 2 * U[i][j] + U[i - 1][j]) / (dx * dx);
			du2y2 = (U[i][j + 1] - 2 * U[i][j] + U[i][j - 1]) / (dy * dy);
			a = (U[i][j] + U[i + 1][j]) / 2;
			b = (U[i - 1][j] + U[i][j]) / 2;
			du2dx = (a * a - b * b
					+ alpha
							* (fabs(a) * ((U[i][j] - U[i + 1][j]) / 2)
									- fabs(b) * ((U[i - 1][j] - U[i][j]) / 2)))
					/ dx;
			duvy =
					((V[iVforFCalc][jVforFCalc] + V[iVforFCalc + 1][jVforFCalc])
							* (U[i][j] + U[i][j + 1])
							- (V[iVforFCalc][jVforFCalc - 1]
									+ V[iVforFCalc + 1][jVforFCalc - 1])
									* (U[i][j - 1] + U[i][j])
							+ alpha
									* (fabs(
											V[iVforFCalc][jVforFCalc]
													+ V[iVforFCalc + 1][jVforFCalc])
											* (U[i][j] - U[i][j + 1])
											- fabs(
													V[iVforFCalc][jVforFCalc - 1]
															+ V[iVforFCalc + 1][jVforFCalc
																	- 1])
													* (U[i][j - 1] - U[i][j])))
							/ (4 * dy);
			F[i][j] = U[i][j]
					+ dt * ((du2x2 + du2y2) * (1 / Re) - du2dx - duvy + GX);
		}
	}
	int iUforGCalc, jUforGCalc;
	for (int i = 1; i < iMaxV; i++) {
		for (int j = 2; j < jMaxV - 1; j++) {
			iUforGCalc = i + 1;
			jUforGCalc = j - 1;
			dv2y2 = (V[i][j + 1] - 2 * V[i][j] + V[i][j - 1]) / (dy * dy);
			dv2x2 = (V[i + 1][j] - 2 * V[i][j] + V[i - 1][j]) / (dx * dx);
			c = (V[i][j] + V[i][j + 1]) / 2;
			d = (V[i][j - 1] + V[i][j]) / 2;
			dv2dy = (c * c - d * d
					+ alpha
							* (fabs(c) * ((V[i][j] - V[i][j + 1]) / 2)
									- fabs(d) * ((V[i][j - 1] - V[i][j]) / 2)))
					/ dy;
			duvx =
					((U[iUforGCalc][jUforGCalc] + U[iUforGCalc][jUforGCalc + 1])
							* (V[i][j] + V[i + 1][j])
							- (U[iUforGCalc - 1][jUforGCalc]
									+ U[iUforGCalc - 1][jUforGCalc + 1])
									* (V[i - 1][j] + V[i][j])
							+ alpha
									* (fabs(
											U[iUforGCalc][jUforGCalc]
													+ U[iUforGCalc][jUforGCalc
															+ 1])
											* (V[i][j] - V[i + 1][j])
											- fabs(
													U[iUforGCalc - 1][jUforGCalc]
															+ U[iUforGCalc - 1][jUforGCalc
																	+ 1])
													* (V[i - 1][j] - V[i][j])))
							/ (4 * dy);
			G[i][j] = V[i][j]
					+ dt * ((dv2x2 + dv2y2) * (1 / Re) - dv2dy - duvx + GY);
		}
	}

	return;
}

void calculate_rs(double dt, double dx, double dy, int imax, int jmax,
		double **F, double **G, double **RS, int iMaxUF, int jMaxUF, int iMaxVG,
		int jMaxVG, int iMaxRS, int jMaxRS) {
	int i, j, iF, iG, jF, jG;
	for (i = 0; i < iMaxRS; i++) {
		for (j = 0; j < jMaxRS; j++) {
			iF = i + 2;
			iG = i + 1;
			jF = j + 1;
			jG = j + 2;
			RS[i][j] = ((F[iF][jF] - F[iF - 1][jF]) / dx
					+ (G[iG][jG] - G[iG][jG - 1]) / dy) / dt;
		}
	}
	return;
}

void calculate_dt(double Re, double tau, double *dt, double dx, double dy,
		int imax, int jmax, double **U, double **V, int iMaxUF, int jMaxUF,
		int iMaxVG, int jMaxVG) {
	double dt1, dt2, dt3, Umax, Vmax;
	double U1 = fabs(U[2][1]);
	double V1 = fabs(V[1][2]);
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	for (int c = 2; c < iMaxUF - 1; c++) {
		for (int d = 1; d < jMaxUF - 1; d++) {
			if (fabs(U[c][d]) > fabs(U1))
				U1 = U[c][d];
		}
	}
	for (int c = 1; c < iMaxVG - 1; c++) {
		for (int d = 2; d < jMaxVG - 1; d++) {
			if (fabs(V[c][d]) > fabs(V1))
				V1 = V[c][d];
		}
	}
	MPI_Reduce(&U1, &Umax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&V1, &Vmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	if (myrank == 0) {

		dt1 = 0.5 * Re / (1 / (dx * dx) + 1 / (dy * dy));
		dt2 = dx / fabs(Umax);
		dt3 = dy / fabs(Vmax);

		*dt = dt1;
		if (dt2 < dt1) {
			*dt = dt2;
			if (dt3 < dt2) {
				*dt = dt3;
			}
		} else if (dt3 < dt1) {
			*dt = dt3;
		}
		*dt = *dt * tau;
		MPI_Bcast(dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	Programm_Sync("New dt Value Received");
	return;
}

void calculate_uv(double dt, double dx, double dy, int imax, int jmax,
		double **U, double **V, double **F, double **G, double **P, int iMaxUF,
		int jMaxUF, int iMaxVG, int jMaxVG) {
	int iPU, jPV;
	int iMaxU = iMaxUF - 1;
	int jMaxU = jMaxUF - 1;
	int iMaxV = iMaxVG - 1;
	int jMaxV = jMaxVG - 1;
	for (int i = 2; i < iMaxU - 1; i++) {
		for (int j = 1; j < jMaxU; j++) {
			iPU = i - 1;
			U[i][j] = F[i][j] - dt * (P[iPU + 1][j] - P[iPU][j]) / dx;
		}
	}
	for (int i = 1; i < iMaxV; i++) {
		for (int j = 2; j < jMaxV - 1; j++) {
			jPV = j - 1;
			V[i][j] = G[i][j] - dt * (P[i][jPV + 1] - P[i][jPV]) / dy;
		}
	}
	return;
}
