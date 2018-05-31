#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>

void pressure_MPI_SndRcv(double **P, int opDirn, int lowBound, int upBound,
		int sndID, int sndColIdx, double *bufSend, int rcvID, int rcvColIdx,
		double *bufRecv, int elCnt, MPI_Status *status);
/*Called within init_parallel*/

void init_parallel (int iproc, int jproc, int imax, int jmax,
		int *myrank, int *il, int *ir, int *jb, int *jt,
		int *rank_l, int *rank_r, int *rank_b, int *rank_t,
		int *omg_i, int *omg_j, int num_proc);
/* Initialises Each Sub-domain/Process */


void Program_Message(char *txt);
/* produces a stderr text output  */



void Programm_Sync(char *txt);
/* produces a stderr textoutput and synchronize all processes */



void Programm_Stop(char *txt);
/* all processes will produce a text output, be synchronized and finished */

