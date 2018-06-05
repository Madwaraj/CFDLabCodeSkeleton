#include "helper.h"
#include "visual.h"
#include "init.h"
#include "uvp.h"
#include "boundary_val.h"
#include "sor.h"
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <mpi.h>
#include "parallel.h"
#include <math.h>

#include <unistd.h> //For debugging commands

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

	// Reading the problem data
	const char* problem_data = "cavity100.dat";

	// Geometry Data

	double xlength; /* length of the domain x-dir.*/
	double ylength; /* length of the domain y-dir.*/
	int imax; /* number of cells x-direction*/
	int jmax; /* number of cells y-direction*/
	double dx; /* length of a cell x-dir. */
	double dy; /* length of a cell y-dir. */
	int il; /*Left bound, x-axis*/
	int ir; /*Right bound, x-axis*/
	int jt; /*Upper bound, y-axis*/
	int jb; /*Lower bound, y-axis*/

	// Time Stepping Data

	double t = 0;
	double tau;
	double t_end; /* end time */
	double dt; /* time step */
	double dt_value; /* time for output */
	int n = 0;

	// Pressure Iteration Data

	int itermax; /* max. number of iterations  */
	double eps; /* accuracy bound for pressure*/
	double omg; /* relaxation factor */
	double alpha; /* uppwind differencing factor*/

	// Problem dependent quantities

	double Re; /* reynolds number   */
	double UI; /* velocity x-direction */
	double VI; /* velocity y-direction */
	double PI; /* pressure */
	double GX; /* gravitation x-direction */
	double GY; /* gravitation y-direction */

	//int data;

	char *message = "X";

	// MPI data
	int myrank;
	int sndrank;
	int omg_i;
	int omg_j;
	int iproc;
	int jproc;
	int num_proc;
	int rank_l;
	int rank_r;
	int rank_t;
	int rank_b;
	MPI_Status status;
	int chunk = 0;
	int iChunk;
	int jChunk;
	int lastUsedir;
	int lastUsedjt;
	MPI_Init(&argn, &args);
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	// Extracting parameter values from data file and assigning them to variables
	read_parameters(problem_data, &Re, &UI, &VI, &PI, &GX, &GY, &t_end,
			&xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau,
			&itermax, &eps, &dt_value, &iproc, &jproc);

	printf("P%d\t Parameters Extracted \n \n", myrank);

	if (myrank == 0) {
		/*Debug Code
		 int sleepVar = 0;
		 char hostname[256];
		 gethostname(hostname, sizeof(hostname));
		 printf("PID %d with rank=%d on %s ready for attach\n", getpid(), myrank,
		 hostname);
		 fflush(stdout);
		 while (0 == sleepVar)
		 sleep(5);
		 End of Debug Code*/

		mkdir("Solution", 0777);
		// Segregating the large domain into smaller sub-domains based on number of processes
		if (iproc == 0 && jproc == 0) {
			if (num_proc % 2 == 0) {
				iproc = (int) sqrt(num_proc);
				jproc = num_proc - omg_i;
			} else if (num_proc % 2 == 1) {
				float temp = sqrt(num_proc);
				if (temp * temp == num_proc) // To check if the number is perfect square
						{
					iproc = temp;
					jproc = temp;
				} else {
					iproc = num_proc;
					jproc = 1;
				}
			} else {
				printf("\n Invalid Domain Size, Exiting!!");
				//break;
			}
		}

		iChunk = imax / iproc;
		jChunk = jmax / jproc;

		// Assign values for il, ir, jb, jt
		sndrank = 0;
		//int *bufTemp[6];
		int *bufTemp = malloc(6 * sizeof(int));
		printf("P%d\t Temporary Buffer created \n \n", myrank);
		for (int j = 1; j < jproc + 1; j++) {
			for (int i = 1; i < iproc + 1; i++) {
				bufTemp[0] = i; //omg_i
				bufTemp[1] = j; //omg_j
				if (i != iproc) {
					bufTemp[2] = ((i - 1) * iChunk) + 1; //il
					bufTemp[3] = bufTemp[2] + iChunk - 1; // ir
					lastUsedir = bufTemp[3];
				} else {
					bufTemp[2] = lastUsedir + 1; //il
					bufTemp[3] = imax; //il
				}
				if (j != jproc) {
					bufTemp[4] = ((j - 1) * jChunk) + 1; //jb
					bufTemp[5] = bufTemp[4] + jChunk - 1; //jt
					lastUsedjt = bufTemp[5];
				} else {
					bufTemp[4] = lastUsedjt + 1; //jb
					bufTemp[5] = jmax; //jt
				}
				if (sndrank != 0) {
					MPI_Send(bufTemp, 6, MPI_INT, sndrank, 1,
					MPI_COMM_WORLD);
				} else { // Assign values for Master Process
					omg_i = bufTemp[0];
					omg_j = bufTemp[1];
					il = bufTemp[2];
					ir = bufTemp[3];
					jb = bufTemp[4];
					jt = bufTemp[5];
				}
				sndrank++;
			}
		}

		//printf("P%d\t Bounds of sub-domains assigned \n \n", myrank);

		//Assign Neighbours for Master Process
		rank_l = MPI_PROC_NULL;
		rank_b = MPI_PROC_NULL;
		if (iproc > 1) { //If there are divisions in the horizontal direction
			rank_r = 1;
		} else { // If there are no horizontal divisions
			rank_r = MPI_PROC_NULL;
		}
		if (jproc > 1) { //Ensures there are divisions in the vertical direction
			rank_t = myrank + iproc;
		} else { // If there are no vertical divisions
			rank_t = MPI_PROC_NULL;
		}
		//printf("L=%d R=%d B=%d T=%d",rank_l,rank_r,rank_b,rank_t);
		//Programm_Sync("");

	}
	// End of work for Master Thread

	else { // Only Slave Threads
		init_parallel(iproc, jproc, imax, jmax, &myrank, &il, &ir, &jb, &jt,
				&rank_l, &rank_r, &rank_b, &rank_t, &omg_i, &omg_j, num_proc); //Initialising the parallel processes
	}
	Programm_Sync("");
	//printf("P%d\t Parallel Processes initialized \n \n", myrank);

	//Matrix extents. Includes Ghost cells + extra cells for U,V,F,G
	int iMaxUF = (ir + 1) - (il - 2) + 1;
	//printf("P%d\t Main: iMaxUF=%d initialized \n \n", myrank,iMaxUF);
	int jMaxUF = (jt + 1) - (jb - 1) + 1;
	//printf("P%d\t Main: jMaxUF=%d initialized \n \n", myrank,jMaxUF);
	int iMaxVG = (ir + 1) - (il - 1) + 1;
	//printf("P%d\t Main: iMaxVG=%d initialized \n \n", myrank,iMaxVG);
	int jMaxVG = (jt + 1) - (jb - 2) + 1;
	//printf("P%d\t Main: jMaxVG=%d initialized \n \n", myrank,jMaxVG);
	int iMaxRS = ir - il + 1;
	//printf("P%d\t Main: iMaxRS=%d initialized \n \n", myrank,iMaxRS);
	int jMaxRS = jt - jb + 1;
	//printf("P%d\t Main: jMaxRS=%d initialized \n \n", myrank,jMaxRS);

	chunk = max(max(iMaxUF - 2, jMaxUF - 2), max(iMaxVG, jMaxVG));
	//printf("P%d\t Main: Max Buff =%d initialized \n \n", myrank, maxBuf);
	double *bufSend = malloc(chunk * sizeof(double));
	//printf("P%d\t Main: bufSend initialized \n \n", myrank);
	double *bufRecv = malloc(chunk * sizeof(double));
	//printf("P%d\t Main: bufRecv initialized \n \n", myrank);

	/* Dynamic allocation of matrices for P(pressure), U(velocity_x), V(velocity_y), F, and G on heap*/
	double **P = matrix(0, iMaxVG - 1, 0, jMaxUF - 1);
	//printf("P%d\t Main: **P initialized \n \n", myrank);
	double **U = matrix(0, iMaxUF - 1, 0, jMaxUF - 1);
	//printf("P%d\t Main: **U initialized \n \n", myrank);
	double **V = matrix(0, iMaxVG - 1, 0, jMaxVG - 1);
	//printf("P%d\t Main: **V initialized \n \n", myrank);
	double **F = matrix(0, iMaxUF - 1, 0, jMaxUF - 1);
	//printf("P%d\t Main: **F initialized \n \n", myrank);
	double **G = matrix(0, iMaxVG - 1, 0, jMaxVG - 1);
	//printf("P%d\t Main: **G initialized \n \n", myrank);
	double **RS = matrix(0, iMaxRS - 1, 0, jMaxRS - 1);
	//printf("P%d\t Main: P, U, V, F, G, RS Matrices initialized \n \n", myrank);

	/* Dynamic allocation of matrices for P(pressure), U(velocity_x), V(velocity_y), F, and G on heap
	 double **P = matrix((il - 1), (ir + 1), (jb - 1), (jt + 1));
	 double **U = matrix((il - 2), (ir + 1), (jb - 1), (jt + 1));
	 double **V = matrix((il - 1), (ir + 1), (jb - 2), (jt + 1));
	 double **F = matrix((il - 2), (ir + 1), (jb - 1), (jt + 1));
	 double **G = matrix((il - 1), (ir + 1), (jb - 2), (jt + 1));
	 double **RS = matrix(il, ir, jb, jt);*/

	//Initialize U, V and P
	init_uvp(UI, VI, PI, U, V, P, iMaxUF, jMaxUF, iMaxVG, jMaxVG);

	//printf("P%d\t U V P Initialized \n \n", myrank);

	int n1 = 0;

	while (t < t_end) {


		//printf("P%d\t Setting Domain BCs \n \n", myrank);
		//Programm_Sync("Got here: Main - Entering boundaryvalues");
		boundaryvalues(imax, jmax, U, V, iMaxUF, jMaxUF, iMaxVG, jMaxVG, rank_l,
				rank_r, rank_b, rank_t); // Assigning Domain Boundary Values

		//printf("P%d\t Domain Boundary values set \n \n", myrank);

		calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, U, V, F, G,
				iMaxUF, jMaxUF, iMaxVG, jMaxVG, rank_l, rank_r, rank_b, rank_t); // Computing Fn and Gn

		//printf("P%d\t Fn & Gn Calculated \n \n", myrank);

		calculate_rs(dt, dx, dy, imax, jmax, F, G, RS, iMaxUF, jMaxUF, iMaxVG,
				jMaxVG, iMaxRS, jMaxRS); // Computing the right hand side of the Pressure Eqn

		//printf("P%d\tRHS Calculated \n \n", myrank);

		int it = 0;

		double res = 1.0; // Residual for the SOR

		while (it < itermax && res > eps) {
			sor(omg, dx, dy, imax, jmax, P, RS, &res, il, ir, jb, jt, iMaxUF,
					jMaxUF, iMaxVG, jMaxVG, rank_l, rank_r, rank_b, rank_t,
					bufSend, bufRecv, chunk); // Successive over-realaxation to solve the Pressure Eqn
			it++;
		}

		Program_Message("SOR Converged \n \n");

		calculate_uv(dt, dx, dy, imax, jmax, U, V, F, G, P, iMaxUF, jMaxUF,
				iMaxVG, jMaxVG); // Computing U, V for the next time-step

		Program_Message(" U V for next time step done\n \n");

		uv_comm(U, V, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, bufSend,
				bufRecv, &status, chunk);
		printf("P%d\t U V values exchanged across boundaries \n \n", myrank);
		char output_file[40];
		sprintf(output_file, "Solution/Output_%d_",myrank);
		if (t >= n1 * dt_value *dt) {
			/*write_vtkFile(output_dir, n, xlength, ylength, iMaxVG - 2,
					jMaxUF - 2, dx, dy, U, V, P);*/

			output_uvp(U, V, P, il, ir, jb, jt, omg_i, omg_j, n, dx, dy, output_file);
			printf("%f Time Elapsed \n", n1 * dt_value);
			n1++;
		}
		/*Programm_Sync("new dt found Complete");
		int sleepVar = 0;
		while (0 == sleepVar)
			sleep(5);*/
		calculate_dt(Re, tau, &dt, dx, dy, imax, jmax, U, V, iMaxUF, jMaxUF,
								iMaxVG, jMaxVG);
		t = t + dt;
		n++;
	}

	//Free memory
	free_matrix(P, 0, iMaxVG - 1, 0, jMaxUF - 1);
	free_matrix(U, 0, iMaxUF - 1, 0, jMaxUF - 1);
	free_matrix(V, 0, iMaxVG - 1, 0, jMaxVG - 1);
	free_matrix(F, 0, iMaxUF - 1, 0, jMaxUF - 1);
	free_matrix(G, 0, iMaxVG - 1, 0, jMaxVG - 1);
	free_matrix(RS, 0, iMaxRS - 1, 0, jMaxRS - 1);

	Programm_Sync("Program end reached");

	Programm_Stop(message);

	return -1;
}

