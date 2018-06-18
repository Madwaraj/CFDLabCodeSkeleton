#include "helper.h"
#include "visual.h"
#include "init.h"
#include"uvp.h"
#include "precice_adapter.h"
#include"boundary_val.h"
#include"sor.h"
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "adapters/c/SolverInterfaceC.h"

/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed using the matrix() command
 * - create the initial setup init_uvp(), init_flag(), output_uvp()
 * - perform the main loop
 * - trailer: destroy memory allocated and do some statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop the following big steps are done (for some of the 
 * operations a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */
int main(int argn, char** args) {
	int select;
	//char* geometry = (char*) (malloc(sizeof(char) * 6));
	char *filename = "configs/heated-plate.dat";

	/*strcpy(filename, problem);
	 strcat(filename, ".dat");*/
	//scanf("%d", &select);
	//select problem
	/*	const char* filename = "0";
	 switch (select) {
	 case 1:
	 filename = "configs/heated-plate.dat";
	 break;

	 case 2:
	 filename = "configs/convection.dat";
	 break;

	 case 3:
	 filename = "configs/F1-heat-exchange.dat";
	 break;

	 case 4:
	 filename = "configs/F2-heat-exchange.dat";
	 break;
	 }
	 */
	//define parameter variables
	double Re; /* reynolds number   */
	double UI; /* velocity x-direction */
	double VI; /* velocity y-direction */
	double PI; /* pressure */
	double GX; /* gravitation x-direction */
	double GY; /* gravitation y-direction */
	double t_end; /* end time */
	double xlength; /* length of the domain x-dir.*/
	double ylength; /* length of the domain y-dir.*/
	double dt; /* time step */
	double dx; /* length of a cell x-dir. */
	double dy; /* length of a cell y-dir. */
	int imax; /* number of cells x-direction*/
	int jmax; /* number of cells y-direction*/
	double alpha; /* uppwind differencing factor*/
	double omg; /* relaxation factor */
	double tau; /* safety factor for time step*/
	int itermax; /* max. number of iterations  */
	/* for pressure per time step */
	double eps; /* accuracy bound for pressure*/
	double dt_value; /* time for output */
	double Pr;
	double TI;
	double beta;
	double x_origin;
	double y_origin;
	double *temperature;
	char *problem = (char*) malloc(100 * sizeof(char));
	char *geometry = (char*) malloc(100 * sizeof(char));
	char *precice_config = (char*) malloc(100 * sizeof(char));
	char *participant_name = (char*) malloc(100 * sizeof(char));
	char *mesh_name = (char*) malloc(100 * sizeof(char));
	char *read_data_name = (char*) malloc(100 * sizeof(char));
	char *write_data_name = (char*) malloc(100 * sizeof(char));
	//Read and assign the parameter values from file
	read_parameters(filename, &imax, &jmax, &xlength, &ylength, &dt, &t_end,
			&tau, &dt_value, &eps, &omg, &alpha, &itermax, &GX, &GY, &Re, &Pr,
			&UI, &VI, &PI, &TI, &beta, &dx, &dy, problem, geometry,
			precice_config, participant_name, mesh_name, read_data_name,
			write_data_name);
	//include_temp =1 => include temperature equations for solving
	int include_temp = 1;

	/*printf("1st checkpoint\n");
	 printf(
	 "geometry=%s\n precice_config=%s\n participant_name=%s\n mesh_name=%s\n read_data_name=%s\n write_data_name=%s\n",
	 geometry, precice_config, participant_name, mesh_name,
	 read_data_name, write_data_name);*/
	//Allocate the matrices for P(pressure), U(velocity_x), V(velocity_y), F, and G on heap
	printf("PROGRESS: Starting matrix allocation... \n");
	double **P = matrix(0, imax - 1, 0, jmax - 1);
	double **U = matrix(0, imax - 1, 0, jmax - 1);
	double **V = matrix(0, imax - 1, 0, jmax - 1);
	double **F = matrix(0, imax - 1, 0, jmax - 1);
	double **G = matrix(0, imax - 1, 0, jmax - 1);
	double **RS = matrix(0, imax - 1, 0, jmax - 1);
	int **flag = imatrix(0, imax - 1, 0, jmax - 1);
	double **T;
	double **T1;

	precicec_createSolverInterface(participant_name, precice_config, 0, 1);
	int dim = precicec_getDimensions();
	int meshID = precicec_getMeshID(mesh_name);

	T = matrix(0, imax - 1, 0, jmax - 1);
	T1 = matrix(0, imax - 1, 0, jmax - 1);

	printf("PROGRESS: Matrices allocated on heap... \n \n");

	//Initilize flags and get count of num_coupling_cells
	int num_coupling_cells;
	init_flag(problem, geometry, imax, jmax, flag, &num_coupling_cells);
	printf("num_coupling_cells=%d \n", num_coupling_cells);
	//Initialise vertices for preCICE with num_coupling_cells from init_flag

	int *vertexIDs = precice_set_interface_vertices(imax, jmax, dx, dy,
			x_origin, y_origin, num_coupling_cells, meshID, flag);
	//Initialize the U, V and P

	int temperatureID = precicec_getDataID(write_data_name, meshID);
	double* temperatureCoupled = (double*) malloc(
			sizeof(double) * num_coupling_cells);

	int heatFluxID = precicec_getDataID(read_data_name, meshID);
	double* heatflux = (double*) malloc(sizeof(double) * num_coupling_cells);

	double precice_dt = precicec_initialize();
	printf("PROGRESS:precicec_initialize... \n \n");
	precice_write_temperature(imax, jmax, num_coupling_cells, temperature,
			vertexIDs, temperatureID, T, flag);
        printf("precice_write_temperature done\n\n");
	precicec_initialize_data();
        printf("precice_initialize_data done\n\n");
	precicec_readBlockScalarData(heatFluxID, num_coupling_cells, vertexIDs,
			heatflux);
        printf("precice_readscalarblockdata done\n\n");
	init_uvpt(UI, VI, PI, TI, imax, jmax, U, V, P, T, flag);

	//Make solution folder
	struct stat st = { 0 };
	char sol_folder[80];
	sprintf(sol_folder, "Solution_%s", problem);
	if (stat(sol_folder, &st) == -1) {
		mkdir(sol_folder, 0700);
	}

	//VTK File Name Prefix
	char sol_directory[80];
	sprintf(sol_directory, "Solution_%s/sol", problem);

	printf("PROGRESS: Starting the flow simulation...\n");
	double t = 0;
	int n = 0;
	int n1 = 0;

	while (precicec_isCouplingOngoing()) {
        printf("Starting to calculate dt\n\n");
		calculate_dt(Re, tau, &dt, dx, dy, imax, jmax, U, V, Pr, include_temp);
	dt = fmin(dt, precice_dt);
	printf("dt calculated\n\n\n");
		printf("t = %lf, dt = %lf, \n\n",t, dt);

		//Used only if inflow BCs are set in PGM
		spec_boundary_val(imax, jmax, U, V, flag);

		boundaryvalues(imax, jmax, U, V, flag);

		set_coupling_boundary(imax, jmax, dx, dy, heatflux, T, flag);
		
		calculate_temp(T, T1, Pr, Re, imax, jmax, dx, dy, dt, alpha, U, V, flag,
				TI, select);
		printf("Calculated Temperature\n\n");
		calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, U, V, F, G,
				flag, beta, T, include_temp);
		printf("Calculated RHS\n\n");
		calculate_rs(dt, dx, dy, imax, jmax, F, G, RS, flag);

		int it = 0;
		double res = 10.0;

		do {
			sor(omg, dx, dy, imax, jmax, P, RS, &res, flag, UI);
			++it;

		} while (it < itermax && res > eps);
		printf("SOR itertions = %d ,residual = %f \n", it - 1, res);
		if ((it == itermax) && (res > eps)) {
			printf("WARNING: Iteration limit reached before convergence. \n");
		}

		calculate_uv(dt, dx, dy, imax, jmax, U, V, F, G, P, flag);

		precice_write_temperature(imax, jmax, num_coupling_cells, temperature,
				vertexIDs, temperatureID, T, flag);
		precice_dt = precicec_advance(dt);
		precicec_readBlockScalarData(heatFluxID, num_coupling_cells, vertexIDs,
				heatflux);

		if ((t >= n1 * dt_value) && (t != 0.0)) {
			write_vtkFile(sol_directory, n, xlength, ylength, imax, jmax, dx,
					dy, U, V, P, T, include_temp);

			printf("writing result at %f seconds \n", n1 * dt_value);
			n1 = n1 + 1;
			continue;
		}
		t = t + dt;
		n = n + 1;
	}
	precicec_finalize();

	printf("PROGRESS: flow simulation completed...\n \n");

	printf("PROGRESS: Freeing allocated memory...\n");
	//Free memory
	free_matrix(P, 0, imax - 1, 0, jmax - 1);
	free_matrix(U, 0, imax - 1, 0, jmax - 1);
	free_matrix(V, 0, imax - 1, 0, jmax - 1);
	free_matrix(F, 0, imax - 1, 0, jmax - 1);
	free_matrix(G, 0, imax - 1, 0, jmax - 1);
	free_matrix(RS, 0, imax - 1, 0, jmax - 1);
	free_imatrix(flag, 0, imax - 1, 0, jmax - 1);
	free_matrix(T, 0, imax - 1, 0, jmax - 1);
	free_matrix(T1, 0, imax - 1, 0, jmax - 1);
	free(geometry);

	printf("PROGRESS: allocated memory released...\n \n");

	printf("PROGRESS: End of Run.\n");
	return -1;

}

/*Things to do:
 read parameter
 */
