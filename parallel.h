#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>

void init_parallel (int iproc, int jproc, int imax, int jmax,
		int *myrank, int *il, int *ir, int *jb, int *jt,
		int *rank_l, int *rank_r, int *rank_b, int *rank_t,
		int *omg_i, int *omg_j, int num_proc);
/* Initialises Each Sub-domain/Process */



void pressure_MPI_SndRcv(double **P, int opDirn, int lowBound, int upBound,
		int sndID, int sndColIdx, double *bufSend, int rcvID, int rcvColIdx,
		double *bufRecv, int elCnt, MPI_Status *status);
/*Called within pressure_comm*/



void pressure_comm(double **P, int il, int ir, int jb, int jt, int rank_l,
		int rank_r, int rank_b, int rank_t, double *bufSend, double *bufRecv,
		MPI_Status *status, int chunk);
/*Exchanges the pressures on sub-domain edges to the respective neighbouring ghost cells*/

void uv_MPI_SndRcv(double **U, double **V,
		int opDirn, //controls direction of send-receive
		int lowBound, int upBound, //Lower and upper bounds for 1VelComp loops
		int lowBound2, int upBound2, //Lower and upper bounds for 2VelComp loops
		int sndID, //Rank of process to which data must be sent
		int sndColIdx, int sndColIdx2, // Indices of fixed row/col of the sending data
		double *bufSend,
		int rcvID, //Rank of process from which data must be received
		int rcvColIdx, int rcvColIdx2, // Indices of fixed row/col of the receiving data
		double *bufRecv, int elCnt, MPI_Status *status);
/*Called within uv_comm*/

void uv_comm(double**U, double**V, int il, int ir, int jb, int jt, int rank_l,
		int rank_r, int rank_b, int rank_t, double *bufSend, double *bufRecv,
		MPI_Status *status, int chunk);
/* Exchanges the U and V values on sub-domain edges to the respective neighbouring ghost cells */

void Program_Message(char *txt);
/* produces a stderr text output  */



void Programm_Sync(char *txt);
/* produces a stderr textoutput and synchronize all processes */



void Programm_Stop(char *txt);
/* all processes will produce a text output, be synchronized and finished */

