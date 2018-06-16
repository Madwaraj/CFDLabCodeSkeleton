#include "precice_adapter.h"

#include "adapters/c/SolverInterfaceC.h"
#include <stdlib.h>

int *precice_set_interface_vertices(int imax, int jmax, double dx, double dy,
		double x_origin, double y_origin, int num_coupling_cells, int meshID,
		int **flag) {
	int dimension = 3;
	int coupledcell = 0;
	double* vertices = (double*)malloc(sizeof(double) * num_coupling_cells * dimension);
	int* vertexIDs = (int*)malloc(num_coupling_cells * sizeof(int));

	for (int j = 1; j < jmax + 1; j++) {
		if (~(flag[imax + 1][j] & (0 << 0 | 0 << 1 | 0 << 2 | 0 << 3 | 0 << 4))
				&& (flag[imax + 1][j] & 1 << 8)) {      // Right boundary
			vertices[dimension * coupledcell] = x_origin + (imax * dx);
			vertices[dimension * coupledcell + 1] = y_origin + (j - 0.5) * dy;
			vertices[dimension * coupledcell + 2] = 0;
			coupledcell++;
		}
	}
	for (int j = 1; j < jmax + 1; j++) {
		if (~(flag[0][j] & (0 << 0 | 0 << 1 | 0 << 2 | 0 << 3 | 0 << 4))
				&& (flag[0][j] & 1 << 7)) {     // Left boundary
			vertices[dimension * coupledcell] = 0;
			vertices[dimension * coupledcell + 1] = y_origin + (j - 0.5) * dy;
			vertices[dimension * coupledcell + 2] = 0;
			coupledcell++;
		}
	}
	for (int i = 1; i < imax + 1; i++) {
		if (~(flag[i][0] & (0 << 0 | 0 << 1 | 0 << 2 | 0 << 3 | 0 << 4))
				&& (flag[i][0] & 1 << 6)) {    // Bottom boundary
			vertices[dimension * coupledcell] = x_origin + (i - 0.5) * dx;
			vertices[dimension * coupledcell + 1] = 0;
			vertices[dimension * coupledcell + 2] = 0;
			coupledcell++;
		}
	}
	for (int i = 1; i < imax + 1; i++) {
		if (~(flag[i][jmax + 1] & (0 << 0 | 0 << 1 | 0 << 2 | 0 << 3 | 0 << 4))
				&& (flag[i][jmax + 1] & 1 << 5)) {   // Top boundary
			vertices[dimension * coupledcell] = x_origin + (i - 0.5) * dx;
			vertices[dimension * coupledcell + 1] = y_origin + (jmax * dy);
			vertices[dimension * coupledcell + 2] = 0;

			coupledcell++;
		}
	}
	precicec_setMeshVertices(meshID, num_coupling_cells, vertices, vertexIDs);
	return vertexIDs;
}

void precice_write_temperature(int imax, int jmax, int num_coupling_cells,
		double *temperature, int *vertexIDs, int temperatureID, double **T,
		int **flag) {
	int count = 0;

	for (int j = 1; j < jmax + 1; j++) {
		if (~(flag[imax + 1][j] & (0 << 0 | 0 << 1 | 0 << 2 | 0 << 3 | 0 << 4))
				&& (flag[imax + 1][j] & 1 << 8)) { // Right boundary
			temperature[count] = T[imax][j];
			count++;
		}
	}
	for (int j = 1; j < jmax + 1; j++) {
		if (~(flag[0][j] & (0 << 0 | 0 << 1 | 0 << 2 | 0 << 3 | 0 << 4))
				&& (flag[0][j] & 1 << 7)) { // Left boundary
			temperature[count] = T[1][j];
			count++;
		}
	}

	for (int i = 1; i < imax + 1; i++) {
		if (~(flag[i][0] & (0 << 0 | 0 << 1 | 0 << 2 | 0 << 3 | 0 << 4))
				&& (flag[i][0] & 1 << 6)) { // Bottom boundary
			temperature[count] = T[i][1];
			count++;
		}
	}

	for (int i = 1; i < imax + 1; i++) {
		if (~(flag[i][jmax + 1] & (0 << 0 | 0 << 1 | 0 << 2 | 0 << 3 | 0 << 4))
				&& (flag[i][jmax + 1] & 1 << 5)) { // Top boundary
			temperature[count] = T[i][jmax];
			count++;
		}
	}
	precicec_writeBlockScalarData(temperatureID, num_coupling_cells, vertexIDs,
			temperature);
}

void write_checkpoint(double time, double **U, double **V, double **T,
		double *time_cp, double **U_cp, double **V_cp, double **T_cp, int imax,
		int jmax) {
	time_cp = &time;

	for (int i = 1; i <= imax; i++) {
		for (int j = 1; j <= jmax; j++) {
			T_cp[i][j] = T[i][j];
			U_cp[i][j] = U[i][j];
			V_cp[i][j] = V[i][j];
		}
	}
}

void restore_checkpoint(double *time, double **U, double **V, double **T,
		double time_cp, double **U_cp, double **V_cp, double **T_cp, int imax,
		int jmax) {

	time = &time_cp;

	for (int i = 1; i <= imax; i++) {
		for (int j = 1; j <= jmax; j++) {
			T[i][j] = T_cp[i][j];
			U[i][j] = U_cp[i][j];
			V[i][j] = V_cp[i][j];
		}
	}
}

void set_coupling_boundary(int imax, int jmax, double dx, double dy,
		double *heatflux, double **T, int **flag) {
	int count = 0;

	for (int j = 1; j < jmax + 1; j++) {
		if (~(flag[imax + 1][j] & (0 << 0 | 0 << 1 | 0 << 2 | 0 << 3 | 0 << 4))
				&& (flag[imax + 1][j] & 1 << 8)) { // Right boundary
			T[imax+1][j] = T[imax][j]+(dx*heatflux[count]);
			count++;
		}
	}

	for (int j = 1; j < jmax + 1; j++) {
		if (~(flag[0][j] & (0 << 0 | 0 << 1 | 0 << 2 | 0 << 3 | 0 << 4))
				&& (flag[0][j] & 1 << 7)) { // Left boundary
			T[0][j] = T[1][j]+(dx*heatflux[count]);
			count++;
		}
	}

	for (int i = 1; i < imax + 1; i++) {
		if (~(flag[i][0] & (0 << 0 | 0 << 1 | 0 << 2 | 0 << 3 | 0 << 4))
				&& (flag[i][0] & 1 << 6)) { // Bottom boundary
			T[i][0] = T[i][1]+(dy*heatflux[count]);
			count++;
		}
	}

	for (int i = 1; i < imax + 1; i++) {
		if (~(flag[i][jmax + 1] & (0 << 0 | 0 << 1 | 0 << 2 | 0 << 3 | 0 << 4))
				&& (flag[i][jmax + 1] & 1 << 5)) { // Top boundary
			T[i][jmax+1] = T[i][jmax]+(dy*heatflux[count]);
			count++;
		}
	}
}
