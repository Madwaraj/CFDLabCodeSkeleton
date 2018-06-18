#include "helper.h"
#include "init.h"
#include <stdio.h>

void read_parameters(const char *szFileName, /* name of the file */
int *imax, /* number of cells x-direction*/
int *jmax, /* number of cells y-direction*/
double *xlength, /* length of the domain x-dir.*/
double *ylength, /* length of the domain y-dir.*/
double *dt, /* time step */
double *t_end, /* end time */
double *tau, /* safety factor for time step*/
double *dt_value, /* time for output */
double *eps, /* accuracy bound for pressure*/
double *omg, /* relaxation factor */
double *alpha, /* uppwind differencing factor*/
int *itermax, /* max. number of iterations  */
double *GX, /* gravitation x-direction */
double *GY, /* gravitation y-direction */
double *Re, /* reynolds number   */
double *Pr, double *UI, /* velocity x-direction */
double *VI, /* velocity y-direction */
double *PI, /* pressure */
double *TI, double *beta, double *dx, /* length of a cell x-dir. */
double *dy, /* length of a cell y-dir. */
char *problem, char *geometry,char *precice_config, char *participant_name, char *mesh_name, char *read_data_name, char *write_data_name
) {
	printf("PROGRESS: Reading .dat file... \n");
	//READ_STRING( szFileName, *problem );
	//READ_STRING( szFileName, geometry );

	READ_INT(szFileName, *imax);
	READ_INT(szFileName, *jmax);

	READ_DOUBLE(szFileName, *xlength);
	READ_DOUBLE(szFileName, *ylength);

	READ_DOUBLE(szFileName, *dt);
	READ_DOUBLE(szFileName, *t_end);
	READ_DOUBLE(szFileName, *tau);
	READ_DOUBLE(szFileName, *dt_value);

	READ_DOUBLE(szFileName, *eps);
	READ_DOUBLE(szFileName, *omg);
	READ_DOUBLE(szFileName, *alpha);
	READ_INT(szFileName, *itermax);

	READ_DOUBLE(szFileName, *GX);
	READ_DOUBLE(szFileName, *GY);
	READ_DOUBLE(szFileName, *Re);
	READ_DOUBLE(szFileName, *Pr);

	READ_DOUBLE(szFileName, *UI);
	READ_DOUBLE(szFileName, *VI);
	READ_DOUBLE(szFileName, *PI);
	READ_DOUBLE(szFileName, *TI);
	READ_DOUBLE(szFileName, *beta);

	/*READ_STRING(szFileName, problem);
	READ_STRING(szFileName, geometry);*/

	*dx = *xlength / (double) (*imax);
	*dy = *ylength / (double) (*jmax);

	printf("PROGRESS: .dat file read... \n \n");

}

void init_uvp(double UI, double VI, double PI, int imax, int jmax, double** U,
		double** V, double** P, int** flag) {
	printf("PROGRESS: Starting initialization of U,V,P ... \n");

	for (int i = 0; i < imax + 2; i++) {
		for (int j = 0; j < jmax + 2; j++) {
			if (flag[i][j] & (1 << 0)) {

				U[i][j] = UI;
				V[i][j] = VI;
				P[i][j] = PI;
			}
		}
	}
	printf("PROGRESS: U,V,P matrices initialized... \n \n");
}

void init_uvpt(double UI, double VI, double PI, double TI, int imax, int jmax,
		double** U, double** V, double** P, double** T, int** flag) {
	printf("PROGRESS: Starting initialization of U,V,P,T ... \n");

	for (int i = 0; i < imax + 2; i++) {
		for (int j = 0; j < jmax + 2; j++) {
			if (flag[i][j] & (1 << 0)) {

				U[i][j] = UI;
				V[i][j] = VI;
				P[i][j] = PI;
				T[i][j] = TI;
			}
		}
	}
	printf("PROGRESS: U,V,P,T matrices initialized... \n \n");
}

int isfluid(int pic) {
	if ((pic == 2) || (pic == 3) || (pic == 6)) {
		return 1;
	} else {
		return 0;
	}
}

void call_assert_error() {
	char szBuff[80];
	sprintf(szBuff, "Current geometry is forbidden. Modify .pgm file. \n");
	ERROR(szBuff);
}

// fluid on opposite sides Left and right
int forbidden_LR(int **pic, int i, int j) {
	if ((pic[i - 1][j] == 6) && (pic[i + 1][j] == 6)) {
		return 1;
	} else {
		return 0;
	}
}

// fluid on opposite sides top and bottom
int forbidden_TB(int **pic, int i, int j) {
	if ((pic[i][j + 1] == 6) && (pic[i][j - 1] == 6)) {
		return 1;
	} else {
		return 0;
	}
}

//Avoids any forbidden configuration
void forbid_assert(int imax, int jmax, int **pic) {
	//inner obstacles
	for (int i = 1; i < imax + 1; i++) {
		for (int j = 1; j < jmax + 1; j++) {
			if (pic[i][j] != 6) {
				if (forbidden_LR(pic, i, j) || forbidden_TB(pic, i, j)) {
					call_assert_error();
				}
			}
		}
	}
	//left boundary
	for (int j = 1; j < jmax + 1; j++) {
		if (pic[0][j] != 6) {
			if (forbidden_TB(pic, 0, j)) {
				call_assert_error();
			}
		}
	}
	//right boundary
	for (int j = 1; j < jmax + 1; j++) {
		if (pic[imax - 1][j] != 6) {
			if (forbidden_TB(pic, imax - 1, j)) {
				call_assert_error();
			}
		}
	}
	//top boundary
	for (int i = 1; i < imax + 1; i++) {
		if (pic[i][jmax - 1] != 6) {
			if (forbidden_LR(pic, i, jmax - 1)) {
				call_assert_error();
			}
		}
	}
	//bottom boundary
	for (int i = 1; i < imax + 1; i++) {
		if (pic[i][0] != 6) {
			if (forbidden_LR(pic, i, 0)) {
				call_assert_error();
			}
		}
	}

}

void init_flag(char* problem, char* geometry, int imax, int jmax, int **flag,
		int *num_coupling_cells) {
	printf("PROGRESS: Setting flags... \n");
	int **pic = imatrix(0, imax + 1, 0, jmax + 1);
	int NumCoupCells = 0;
	pic = read_pgm(geometry);
	forbid_assert(imax, jmax, pic); //Checks for disallowed geometries
	for (int i = 0; i < imax + 2; i++) {
		for (int j = 0; j < jmax + 2; j++) {

			flag[i][j] = 0;

			switch (pic[i][j]) {
			case 0: //no-slip
				flag[i][j] = 1 << 1;
				break;

			case 1: //free-slip
				flag[i][j] = 1 << 2;
				break;

			case 2: //outflow
				flag[i][j] = 1 << 3;
				break;

			case 3: //inflow
				flag[i][j] = 1 << 4;
				break;

			case 4: //coupling
				//Implies that the first 5 values (from right to left) are all zeros.
				//This condition is taken as a coupling obstacle cell.
				flag[i][j] = 0 << 0;
				NumCoupCells++;
				break;

			case 6: //fluid
				flag[i][j] = 1 << 0;
				break;
			}

			if (!isfluid(pic[i][j])) //set neighbors if not fluid
					{

				if (i < imax + 1 && pic[i + 1][j] == 6) {
					flag[i][j] |= 1 << 8;  //Set B_O
				}
				if (i > 0 && pic[i - 1][j] == 6) {
					flag[i][j] |= 1 << 7; //Set B_W
				}
				if (j < jmax + 1 && pic[i][j + 1] == 6) {
					flag[i][j] |= 1 << 5; //Set B_N
				}
				if (j > 0 && pic[i][j - 1] == 6) {
					flag[i][j] |= 1 << 6; //Set B_S
				}
			}

		}

	}
	*num_coupling_cells = NumCoupCells;
	free_imatrix(pic, 0, imax + 1, 0, jmax + 1);
	printf("PROGRESS: flags set using .pgm file. num_coupling_cells assigned.\n \n");

}

int num_coupling( char* geometry, int imax, int jmax)
{
	int **pic = imatrix(0,imax-1,0,jmax-1);
	pic = read_pgm(geometry);

	int counter = 0;

		for (int i=0; i<imax; i++) {
		for (int j=0; j<jmax; j++) {

		if (pic[i][j] == 9)
		counter++;

					   }	
					   }

return counter;
}
