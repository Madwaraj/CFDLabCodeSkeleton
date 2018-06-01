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
#include <parallel.h>
#include <math.h>

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
int main(int argn, char** args){


	// Reading the problem data
	const char* problem_data = "cavity100.dat";

	// Geometry Data	

	double xlength;           /* length of the domain x-dir.*/
    	double ylength;           /* length of the domain y-dir.*/
	int  imax;                /* number of cells x-direction*/
	int  jmax;                /* number of cells y-direction*/
	double dx;                /* length of a cell x-dir. */
    	double dy;                /* length of a cell y-dir. */
	int il;
	int ir;
	int jt;
	int jb;
	omg_i;
	omg_j;

	// Time Stepping Data	

    	double t = 0;
    	double tau;
	double t_end;             /* end time */
	double dt;                /* time step */
	double dt_value;          /* time for output */
	int n = 0;

	// Pressure Iteration Data

	int  itermax;             /* max. number of iterations  */
	double eps;               /* accuracy bound for pressure*/
	double omg;               /* relaxation factor */
	double alpha;             /* uppwind differencing factor*/

	// Problem dependent quantities

	double Re;                /* reynolds number   */
    	double UI;                /* velocity x-direction */
    	double VI;                /* velocity y-direction */
    	double PI;                /* pressure */
    	double GX;                /* gravitation x-direction */
    	double GY;                /* gravitation y-direction */

	int data;

	char message = 'X';

	// MPI data

	int myrank;
	int omg_i; 
	int omg_j; 
	int iproc; 
	int jproc; 
	int num_proc;
	int rank_l;
	int rank_r;
	int rank_t;
	int rank_b;
	double bufSend;
	double bufRecv;
	MPI_Status status;
	int chuk;
	

	MPI_Init(int argn, char** args);
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc); // Not sure if this is required here since init_parallel already has it
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	if (rank == 0 )
	{

		// Extracting parameter values from data file and assigning them to variables
		data = read_parameters(problem_data, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, 
					&xlength, &ylength, &dt, &dx, &dy, &imax, &jmax,
                                	&alpha, &omg, &tau, &itermax, &eps, &dt_value, &iproc, &jproc);
							
		data++;

		/* Dynamic allocation of matrices for P(pressure), U(velocity_x), V(velocity_y), F, and G on heap
		double **P = matrix(0, imax+1, 0, jmax+1);
    		double **U = matrix(0, imax, 0, jmax+1);
    		double **V = matrix(0, imax+1, 0, jmax);
    		double **F = matrix(0, imax, 0, jmax+1);
    		double **G = matrix(0, imax+1, 0, jmax);
    		double **RS = matrix(1, imax, 1, jmax); */ //Not Sure if this is necessary


		// Segregating the large domain into smaller sub-domains based on number of processes
		if (num_proc % 2 == 0)
		{
			omg_i = (int)sqrt(num_proc);
			omg_j = num_proc - omg_i;
		}
		else if (num_proc % 2 == 1)
		{
			float temp = sqrt(num_proc);
			if (temp * temp == num_proc) // To check if the number is perfect square
			{
				omg_i = temp;
				omg_j = temp;
			}
			else
			{
				omg_i = num_proc;
				omg_j = 1
			} 	
		}
		else
		{
			printf("\n Invalid Domain Size, Exiting!!");
			break;
		}
		
		// Broadcasting data to all processes

		MPI_Bcast(&Re, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&t_end, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&xlength, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&ylength, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&dx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&dy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&imax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&jmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&omg, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&tau, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&itermax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&eps, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&dt_value, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


		for (i = 0; i < omg_i; i++)
		{
			for(j = 0; j < omg_j; j++)
			{
				il[i][j] = (imax / omg_i) * i;
				ir[i][j] = (imax / omg_i) * (i + 1);
				jt[i][j] = (jmax / omg_j) * j;
				jb[i][j] = (jmax / omg_j) * (j + 1);
			}
		}
	} // End of work for Master Thread

	parallel_init(int iproc, int jproc, int imax, int jmax, int *myrank,
		int *il, int *ir, int *jb, int *jt, int *rank_l, int *rank_r,
		int *rank_b, int *rank_t, int *omg_i, int *omg_j, int num_proc)	//Initializing the parallel processes	
		


	//Initialize U, V and P	
	init_uvp(UI, VI, PI, imax, jmax, U, V, P);
									
	int n1 = 0;

	mkdir("Output",0777);


	while (t < t_end) 
	{

	    	calculate_dt(Re, tau, &dt, dx, dy, imax, jmax, U, V); // Adaptive time stepping

	    	printf("Time Step is %f \n", t);
					
	    	boundaryvalues(imax, jmax, U, V); // Assigning Boundary Values
														
	    	calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, U, V, F, G); // Computing Fn and Gn
													
	    	calculate_rs(dt, dx, dy, imax, jmax, F, G, RS); // Computing the right hand side of the Pressure Eqn
													
		int it = 0;

		double res = 1.0; // Residual for the SOR 

		while(it < itermax && res > eps)
		{
			sor(omg, dx, dy, imax, jmax, P, RS, &res); // Successive over-realaxation to solve the Pressure Eqn    	
			it++;

			void pressure_comm(double **P, int il, int ir, int jb, int jt, int rank_l,
		int rank_r, int rank_b, int rank_t, double *bufSend, double *bufRecv, MPI_Status *status, int chunk) // Exchange Pressure values
		}

		calculate_uv(dt, dx, dy, imax, jmax, U, V, F, G, P); // Computing U, V for the next time-step

		if (t >= n1*dt_value)
  		{
   			write_vtkFile("Output/Solution", n, xlength, ylength, imax, jmax, dx, dy, U, V, P);
			printf("%f Time Elapsed \n", n1*dt_value);
    			n1++;
    			continue;
  		}	

		t = t + dt;

		n++;
	}
	
	//Free memory
    	free_matrix( P, 0, imax+1, 0, jmax+1);
    	free_matrix( U, 0, imax, 0, jmax+1);
    	free_matrix( V, 0, imax+1, 0, jmax);
    	free_matrix( F, 0, imax, 0, jmax+1);
    	free_matrix( G, 0, imax+1, 0, jmax);
    	free_matrix(RS, 1, imax, 1, jmax);

	Program_Stop(*message);
	
  return -1;
}













































