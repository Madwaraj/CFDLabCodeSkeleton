#include "parallel.h"

void init_parallel(int iproc, int jproc, int imax, int jmax, int *myrank,
		int *il, int *ir, int *jb, int *jt, int *rank_l, int *rank_r,
		int *rank_b, int *rank_t, int *omg_i, int *omg_j, int num_proc) {
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	if (0 < myrank && (*il) != 1)
		rank_l = myrank - 1;
	else
		rank_l = MPI_PROC_NULL;
}

void pressure_MPI_SndRcv(double **P, int opDirn, int lowBound, int upBound,
		int sndID, int sndColIdx, double *bufSend, int rcvID, int rcvColIdx,
		double *bufRecv, int elCnt, MPI_Status *status) {
	int *sndP_i, *sndP_j, *rcvP_i, *rcvP_j;
	int ctr, ctr2 = 0;

	/*opDirn controls direction of send-receive
	 *opDirn=1 => //Send to Left Neighbour, Receive from Right Neighbour
	 *opDirn=2 => //Send to Right Neighbour, Receive from Left Neighbour
	 *opDirn=3 => //Send to Bottom Neighbour, Receive from Top Neighbour
	 *opDirn=4 => //Send to Top Neighbour, Receive from Bottom Neighbour
	 */

	switch (opDirn) {
	case 1 || 2:
		sndP_i = &sndColIdx;
		sndP_j = &ctr;
		rcvP_i = &rcvColIdx;
		rcvP_j = &ctr;
		break;
	case 3 || 4:
		sndP_i = &ctr;
		sndP_j = &sndColIdx;
		rcvP_i = &ctr;
		rcvP_j = &rcvColIdx;
		break;
	}

	for (ctr = lowBound; ctr < upBound + 1; ctr++) {
		bufSend[ctr2] = P[(*sndP_i)][(*sndP_j)];
		ctr2++;
	}
	ctr2 = 0;
	MPI_Sendrecv(bufSend, elCnt, MPI_DOUBLE, sndID, 1, bufRecv, elCnt,
	MPI_DOUBLE, rcvID, 2, MPI_COMM_WORLD, status);
	if (status->MPI_SOURCE != MPI_NO_OP) {
		for (ctr = lowBound; ctr < upBound + 1; ctr++) {
			P[(*rcvP_i)][(*rcvP_j)] = bufRecv[ctr2];
			ctr2++;
		}
	}
}

void pressure_comm(double **P, int il, int ir, int jb, int jt, int rank_l,
		int rank_r, int rank_b, int rank_t, double *bufSend, double *bufRecv,
		MPI_Status *status, int chunk) {
	int iBdry_cnt = jt - jb + 1; //Number of elements on a left or right facing boundary
	int jBdry_cnt = ir - il + 1; //Number of elements on a top or bottom facing boundary

	pressure_MPI_SndRcv(P, 1, jb, jt, rank_l, il, bufSend, rank_r, (ir + 1),
			bufRecv, iBdry_cnt, status); //Send to Left Neighbour, Receive from Right Neighbour
	pressure_MPI_SndRcv(P, 2, jb, jt, rank_r, ir, bufSend, rank_l, (il - 1),
			bufRecv, iBdry_cnt, status); //Send to Right Neighbour, Receive from Left Neighbour
	pressure_MPI_SndRcv(P, 3, il, ir, rank_b, jb, bufSend, rank_t, (jt + 1),
			bufRecv, jBdry_cnt, status); //Send to Bottom Neighbour, Receive from Top Neighbour
	pressure_MPI_SndRcv(P, 4, il, ir, rank_t, jt, bufSend, rank_b, (jb - 1),
			bufRecv, jBdry_cnt, status); //Send to Top Neighbour, Receive from Neighbour Bottom

	/*
	 * //Old Approach
	 //Send to Left Neighbour, Receive from Right Neighbour
	 for (ctr=jb;ctr<jt+1;ctr++)
	 {
	 bufSend[ctr2]=P[il][ctr];
	 ctr2++;
	 }
	 ctr2=0;

	 MPI_Sendrecv(bufSend,iBdry_cnt,MPI_DOUBLE,rank_l,1,bufRecv,iBdry_cnt,MPI_DOUBLE,rank_r,2,MPI_COMM_WORLD,status);
	 if(status->MPI_SOURCE!=MPI_NO_OP)
	 {
	 for (ctr=jb;ctr<jt+1;ctr++)
	 {
	 P[ir][ctr]=bufRecv[ctr2];
	 ctr2++;
	 }
	 }
	 ctr2=0;

	//Send to Right Neighbour, Receive from Left Neighbour
	 for (ctr=jb;ctr<jt+1;ctr++)
	 {
	 bufSend[ctr2]=P[ir][ctr];
	 ctr2++;
	 }
	 ctr2=0;
	 MPI_Sendrecv(bufSend,iBdry_cnt,MPI_DOUBLE,rank_r,2,bufRecv,iBdry_cnt,MPI_DOUBLE,rank_l,1,MPI_COMM_WORLD,status);

	//Send to Bottom Neighbour, Receive from Top Neighbour
	 for (ctr=il;ctr<ir+1;ctr++)
	 {
	 bufSend[ctr2]=P[ctr][jb];
	 ctr2++;
	 }
	 ctr2=0;
	 MPI_Sendrecv(bufSend,jBdry_cnt,MPI_DOUBLE,rank_b,3,bufRecv,jBdry_cnt,MPI_DOUBLE,rank_t,4,MPI_COMM_WORLD,status);

	//Send to Top Neighbour, Receive from Neighbour Bottom
	 for (ctr=il;ctr<ir+1;ctr++)
	 {
	 bufSend[ctr2]=P[ctr][jt];
	 ctr2++;
	 }
	 ctr2=0;
	 MPI_Sendrecv(bufSend,jBdry_cnt,MPI_DOUBLE,rank_t,4,bufRecv,jBdry_cnt,MPI_DOUBLE,rank_b,3,MPI_COMM_WORLD,status);*/

}

void Program_Message(char *txt)
/* produces a stderr text output  */

{
	int myrank;

	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	fprintf(stderr, "-MESSAGE- P:%2d : %s\n", myrank, txt);
	fflush(stdout);
	fflush(stderr);
}

void Programm_Sync(char *txt)
/* produces a stderr textoutput and synchronize all processes */

{
	int myrank;

	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Barrier(MPI_COMM_WORLD); /* synchronize output */
	fprintf(stderr, "-MESSAGE- P:%2d : %s\n", myrank, txt);
	fflush(stdout);
	fflush(stderr);
	MPI_Barrier(MPI_COMM_WORLD);
}

void Programm_Stop(char *txt)
/* all processes will produce a text output, be synchronized and finished */

{
	int myrank;

	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Barrier(MPI_COMM_WORLD); /* synchronize output */
	fprintf(stderr, "-STOP- P:%2d : %s\n", myrank, txt);
	fflush(stdout);
	fflush(stderr);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	exit(1);
}
