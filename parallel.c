#include "parallel.h"

#include <unistd.h> //For debugging commands

void init_parallel(int iproc, int jproc, int imax, int jmax, int *myrank,
		int *il, int *ir, int *jb, int *jt, int *rank_l, int *rank_r,
		int *rank_b, int *rank_t, int *omg_i, int *omg_j, int num_proc) {
	//MPI_Comm_rank(MPI_COMM_WORLD, myrank);
	//printf("P%d\t Init_Parallel \n \n", *myrank);
	int bufRcvT[6];
	//printf("P%d\t Init_Parallel: Enter Recv\n \n", *myrank);
	MPI_Recv(bufRcvT, 6, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD,
	MPI_STATUS_IGNORE);
	//printf("P%d\t Init_Parallel: Finished Recv\n \n", *myrank);
	*omg_i = bufRcvT[0];
	/*printf("P%d\t Init_Parallel: Assigned omg_i=%d\n \n", *myrank, *omg_i);*/
	*omg_j = bufRcvT[1];
	/*printf("P%d\t Init_Parallel: Assigned omg_j=%d\n \n", *myrank, *omg_j);*/
	*il = bufRcvT[2];
	/*printf("P%d\t Init_Parallel: Assigned il=%d\n \n", *myrank, *il);*/
	*ir = bufRcvT[3];
	/*	printf("P%d\t Init_Parallel: Assigned ir=%d\n \n", *myrank, *ir);*/
	*jb = bufRcvT[4];
	/*printf("P%d\t Init_Parallel: Assigned jb=%d\n \n", *myrank, *jb);*/
	*jt = bufRcvT[5];
	/*printf("P%d\t Init_Parallel: Assigned jt=%d\n \n", *myrank, *jt);*/

//Assign the Neighbours
	if (iproc > 1) { //If there are divisions in the horizontal direction
		//Left Neighbour
		if (*omg_i > 1) {
			*rank_l = *myrank - 1;
			/*	printf("P%d\t Init_Parallel: Assigned rank_l=%d\n \n", *myrank,
			 *rank_l);*/
		} else {
			*rank_l = MPI_PROC_NULL;
			/*printf("P%d\t Init_Parallel: Assigned rank_l=%d NULL\n \n", *myrank,
			 *rank_l);*/
		}
		//Right Neighbour
		if (*omg_i < iproc) {
			*rank_r = *myrank + 1;
			/*printf("P%d\t Init_Parallel: Assigned rank_r=%d\n \n", *myrank,
			 *rank_r);*/
		} else {
			*rank_r = MPI_PROC_NULL;
			/*printf("P%d\t Init_Parallel: Assigned rank_r=%d NULL\n \n", *myrank,
			 *rank_r);*/
		}
	} else { // If there are no horizontal divisions
		*rank_r = MPI_PROC_NULL;
		*rank_l = MPI_PROC_NULL;
		/*printf(
		 "P%d\t Init_Parallel: Assigned rank_l=%d NULL rank_r=%d NULL\n \n",
		 *myrank, *rank_l, *rank_r);*/
	}
	if (jproc > 1) { //Ensures there are divisions in the vertical direction
		//Bottom Neighbour
		if (*omg_j > 1) {
			*rank_b = *myrank - iproc;
			/*printf("P%d\t Init_Parallel: Assigned rank_b=%d\n \n", *myrank,
			 *rank_b);*/
		} else {
			*rank_b = MPI_PROC_NULL;
			/*printf("P%d\t Init_Parallel: Assigned rank_b=%d NULL\n \n", *myrank,
			 *rank_b);*/
		}
		//Top Neighbour
		if (*omg_j < jproc) {
			*rank_t = *myrank + iproc;
			/*printf("P%d\t Init_Parallel: Assigned rank_t=%d\n \n", *myrank,
			 *rank_t);*/
		} else {
			*rank_t = MPI_PROC_NULL;
			/*printf("P%d\t Init_Parallel: Assigned rank_t=%d NULL\n \n", *myrank,
			 *rank_t);*/
		}
	} else { // If there are no vertical divisions
		*rank_b = MPI_PROC_NULL;
		*rank_t = MPI_PROC_NULL;
		/*printf("P%d\t Init_Parallel: Assigned rank_b=%d NULL rank_t=%d NULL\n \n",
		 *myrank, *rank_b, *rank_t);*/
	}
}

void pressure_MPI_SndRcv(double **P, int opDirn, int lowBound, int upBound,
		int sndID, int sndColIdx, double *bufSend, int rcvID, int rcvColIdx,
		double *bufRecv, int elCnt, MPI_Status *status) {
	int *sndP_i, *sndP_j, *rcvP_i, *rcvP_j;
	int ctr, ctr2 = 0;
	int locRank;

	MPI_Comm_rank(MPI_COMM_WORLD, &locRank);

	/*opDirn controls direction of send-receive
	 *opDirn=1 => //Send to Left Neighbour, Receive from Right Neighbour
	 *opDirn=2 => //Send to Right Neighbour, Receive from Left Neighbour
	 *opDirn=3 => //Send to Bottom Neighbour, Receive from Top Neighbour
	 *opDirn=4 => //Send to Top Neighbour, Receive from Bottom Neighbour
	 */
	//printf("P%d\t Entered SndRcv %d time\n", locRank, opDirn);
	switch (opDirn) {
	case 1:
		sndP_i = &sndColIdx;
		sndP_j = &ctr;
		rcvP_i = &rcvColIdx;
		rcvP_j = &ctr;
		break;

	case 2:
		sndP_i = &sndColIdx;
		sndP_j = &ctr;
		rcvP_i = &rcvColIdx;
		rcvP_j = &ctr;
		break;

	case 3:
		sndP_i = &ctr;
		sndP_j = &sndColIdx;
		rcvP_i = &ctr;
		rcvP_j = &rcvColIdx;
		break;

	case 4:
		sndP_i = &ctr;
		sndP_j = &sndColIdx;
		rcvP_i = &ctr;
		rcvP_j = &rcvColIdx;
		break;
	}
	//printf("P%d\t Checked for the cases\n", locRank);
	for (ctr = lowBound; ctr < upBound; ctr++) {
		bufSend[ctr2] = P[(*sndP_i)][(*sndP_j)];
		ctr2++;
	}
	/*printf("P%d\t ctr2 = %d\t lowBound=%d\t upBound=%d\n", locRank, ctr,
	 lowBound, upBound);*/
	ctr2 = 0;
	// Programm_Sync("Start");
	//printf("P%d: Starting SndRcv %d time\n", locRank, opDirn);
	MPI_Sendrecv(bufSend, elCnt, MPI_DOUBLE, sndID, 1, bufRecv,
			elCnt, MPI_DOUBLE, rcvID, 1, MPI_COMM_WORLD, status);
	//Programm_Sync("Finish ");
	/*printf("P%d\t Exit SndRcv %d time status->MPI_SOURCE=%d errCode=%d\n",
			locRank, opDirn, status->MPI_SOURCE, errCode);*/
	if (status->MPI_SOURCE != MPI_PROC_NULL) {
		//if (MPI_PROC_NULL != rcvID) {
		for (ctr = lowBound; ctr < upBound; ctr++) {
			P[(*rcvP_i)][(*rcvP_j)] = bufRecv[ctr2];
			ctr2++;
		}
	}
}

void pressure_comm(double **P, int il, int ir, int jb, int jt, int rank_l,
		int rank_r, int rank_b, int rank_t, double *bufSend, double *bufRecv,
		MPI_Status *status, int chunk) {
	int iSize = ir - il + 3; //Number of elements on a top or bottom facing boundary
	int jSize = jt - jb + 3; //Number of elements on a left or right facing boundary
	int lRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &lRank);
	//printf("Ilow = %d\t %d \t %d\t %d\t %d\t %d\n", il, ir, jb, jt, jSize, iSize);
	//printf("Entering SndRcv 3rd time\n");
	Programm_Sync("Start Pressure Swap");
	//printf("Entering SndRcv 1st time\n");
	pressure_MPI_SndRcv(P, 1, 1, (jSize - 2), rank_l, 1, bufSend, rank_r,
			(iSize - 1), bufRecv, chunk, status); //Send to Left Neighbour, Receive from Right Neighbour
	//Programm_Sync("Sent Pressures Left");
	//printf("Entering SndRcv 2nd time\n");
	pressure_MPI_SndRcv(P, 2, 1, (jSize - 2), rank_r, (iSize - 2), bufSend,
			rank_l, 0, bufRecv, chunk, status); //Send to Right Neighbour, Receive from Left Neighbour
	//Programm_Sync("Sent Pressures Right");
	pressure_MPI_SndRcv(P, 3, 1, (iSize - 2), rank_b, 1, bufSend, rank_t,
			(jSize - 1), bufRecv, iSize - 2, status); //Send to Bottom Neighbour, Receive from Top Neighbour
	//Programm_Sync("Sent Pressures Top");
	//printf("Entering SndRcv 4th time\n");
	pressure_MPI_SndRcv(P, 4, 1, (iSize - 2), rank_t, (jSize - 2), bufSend,
			rank_b, 0, bufRecv, iSize - 2, status); //Send to Top Neighbour, Receive from Neighbour Bottom


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

void uv_MPI_SndRcv(double **U, double **V, int opDirn, //controls direction of send-receive
		int lowBound, int upBound, //Lower and upper bounds for 1VelComp loops
		int lowBound2, int upBound2, //Lower and upper bounds for 2VelComp loops
		int sndID, //Rank of process to which data must be sent
		int sndColIdx, int sndColIdx2, // Indices of fixed row/col of the sending data
		double *bufSend, int rcvID, //Rank of process from which data must be received
		int rcvColIdx, int rcvColIdx2, // Indices of fixed row/col of the receiving data
		double *bufRecv, int elCnt, MPI_Status *status) {

	int *snd_i, *snd_j, *rcv_i, *rcv_j, *snd_i2, *snd_j2, *rcv_i2, *rcv_j2;
	double ***vel1, ***vel2;
	int ctr, ctr2 = 0;

	/*opDirn controls direction of send-receive
	 *opDirn=1 => //Send to Left Neighbour, Receive from Right Neighbour
	 *opDirn=2 => //Send to Right Neighbour, Receive from Left Neighbour
	 *opDirn=3 => //Send to Bottom Neighbour, Receive from Top Neighbour
	 *opDirn=4 => //Send to Top Neighbour, Receive from Bottom Neighbour
	 */

	switch (opDirn) {
	case 1:
		vel1 = &U;
		vel2 = &V;
		snd_i = &sndColIdx;
		snd_j = &ctr;
		snd_i2 = &sndColIdx2;
		snd_j2 = &ctr;
		rcv_i = &rcvColIdx;
		rcv_j = &ctr;
		rcv_i2 = &rcvColIdx2;
		rcv_j2 = &ctr;
		break;
	case 2:
		vel1 = &U;
		vel2 = &V;
		snd_i = &sndColIdx;
		snd_j = &ctr;
		snd_i2 = &sndColIdx2;
		snd_j2 = &ctr;
		rcv_i = &rcvColIdx;
		rcv_j = &ctr;
		rcv_i2 = &rcvColIdx2;
		rcv_j2 = &ctr;
		break;
	case 3:
		vel1 = &V;
		vel2 = &U;
		snd_i = &ctr;
		snd_j = &sndColIdx;
		snd_i2 = &ctr;
		snd_j2 = &sndColIdx2;
		rcv_i = &ctr;
		rcv_j = &rcvColIdx;
		rcv_i2 = &ctr;
		rcv_j2 = &rcvColIdx2;
		break;
	case 4:
		vel1 = &V;
		vel2 = &U;
		snd_i = &ctr;
		snd_j = &sndColIdx;
		snd_i2 = &ctr;
		snd_j2 = &sndColIdx2;
		rcv_i = &ctr;
		rcv_j = &rcvColIdx;
		rcv_i2 = &ctr;
		rcv_j2 = &rcvColIdx2;

		break;
	}

	// Store 1VelComp values
	for (ctr = lowBound; ctr < upBound; ctr++) {
		bufSend[ctr2] = (*vel1)[(*snd_i)][(*snd_j)];
		ctr2++;
	}
	// Store 2VelComp values
	for (ctr = lowBound2; ctr < upBound2; ctr++) {
		bufSend[ctr2] = (*vel2)[(*snd_i2)][(*snd_j2)];
		ctr2++;
	}
	ctr2 = 0;
	MPI_Sendrecv(bufSend, elCnt, MPI_DOUBLE, sndID, 1, bufRecv, elCnt,
	MPI_DOUBLE, rcvID, MPI_ANY_TAG, MPI_COMM_WORLD, status);
	if (status->MPI_SOURCE != MPI_PROC_NULL) {
		// Receive 1VelComp values
		for (ctr = lowBound; ctr < upBound; ctr++) {
			(*vel1)[(*rcv_i)][(*rcv_j)] = bufRecv[ctr2];
			ctr2++;
		}
		// Receive 2VelComp values
		for (ctr = lowBound2; ctr < upBound2; ctr++) {
			(*vel2)[(*rcv_i2)][(*rcv_j2)] = bufRecv[ctr2];
			ctr2++;
		}
	}
}

void uv_comm(double**U, double**V, int il, int ir, int jb, int jt, int rank_l,
		int rank_r, int rank_b, int rank_t, double *bufSend, double *bufRecv,
		MPI_Status *status, int chunk) {
	int iMaxU = (ir + 1) - (il - 2); //Last i index in U matrix
	int jMaxU = (jt + 1) - (jb - 1); //Last j index in U matrix
	int iMaxV = (ir + 1) - (il - 1); //Last i index in v matrix
	int jMaxV = (jt + 1) - (jb - 2); //Last j index in v matrix

	uv_MPI_SndRcv(U, V, 1, 1, jMaxU, 1, jMaxV, rank_l, 2, 1, bufSend, rank_r,
			iMaxU, iMaxV, bufRecv, chunk, status); //Send to Left Neighbour, Receive from Right Neighbour
	uv_MPI_SndRcv(U, V, 2, 1, jMaxU, 1, jMaxV, rank_r, (iMaxU - 1), (iMaxV - 1),
			bufSend, rank_l, 0, 0, bufRecv, chunk, status); //Send to Right Neighbour, Receive from Left Neighbour
	uv_MPI_SndRcv(U, V, 3, 1, iMaxV, 1, iMaxU, rank_b, 2, 1, bufSend, rank_t,
			jMaxV, iMaxU, bufRecv, chunk, status); //Send to Bottom Neighbour, Receive from Top Neighbour
	uv_MPI_SndRcv(U, V, 4, 1, iMaxV, 1, iMaxU, rank_t, (jMaxV - 2), (jMaxU - 1),
			bufSend, rank_b, 0, 0, bufRecv, chunk, status); //Send to Top Neighbour, Receive from Bottom Neighbour

	/* //Old Approach

	 int iBuf2VelComp, jBuf2VelComp, iBuf1VelComp, jBuf1VelComp;
	 iBuf2VelComp = 2*((jt-jb)+1); //No. of V velocity elements (in the i-passing direction, i.e 2 velocity components passed per cell)
	 jBuf2VelComp = 2*((ir-il)+1); //No. of U velocity elements (in the j-passing direction, i.e 2 velocity components passed per cell)
	 iBuf1VelComp = iBufTotal-iBuf2VelComp; //No. of U velocity elements (in the i-passing direction, i.e one velocity component passed per cell)
	 jBuf1VelComp = jBufTotal-jBuf2VelComp; //No. of V velocity elements (in the j-passing direction, i.e one velocity component passed per cell)

	 // Pass to Left
	 //Store V (2VelComp) values
	 for(i=jb-1;i < jt + 1;i++) {
	 bufSend[ctr] = V[il][i];
	 ctr++;
	 }
	 //Store U (1VelComp) values
	 for (i = jb; i < jt + 1; i++) {
	 bufSend[ctr] = U[il][i];
	 ctr++;
	 }
	 ctr = 0;
	 MPI_Sendrecv(bufSend, iBufTotal, MPI_DOUBLE, rank_l, 1, bufRecv,
	 iBufTotal,
	 MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD, status);
	 //Receive V (2VelComp) values
	 for (i = jb - 1; i < jt + 1; i++) {
	 V[ir + 1][i] = bufRecv[ctr];
	 ctr++;
	 }
	 //Receive U (1VelComp) values
	 for (i = jb; i < jt + 1; i++) {
	 U[ir + 1][i] = bufRecv[ctr];
	 ctr++;
	 }
	 ctr = 0;

	 // Pass to Right
	 //Store V (2VelComp) values
	 for (i = jb - 1; i < jt + 1; i++) {
	 bufSend[ctr] = V[ir][i];
	 ctr++;
	 }
	 //Store U (1VelComp) values
	 for (i = jb; i < jt + 1; i++) {
	 bufSend[ctr] = U[ir - 1][i];
	 ctr++;
	 }
	 ctr = 0;
	 MPI_Sendrecv(bufSend, iBufTotal, MPI_DOUBLE, rank_r, 1, bufRecv,
	 iBufTotal,
	 MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD, status);
	 //Receive V (2VelComp) values
	 for (i = jb - 1; i < jt + 1; i++) {
	 V[il - 1][i] = bufRecv[ctr];
	 ctr++;
	 }
	 //Receive U (1VelComp) values
	 for (i = jb; i < jt + 1; i++) {
	 U[il - 2][i] = bufRecv[ctr];
	 ctr++;
	 }
	 ctr = 0;

	 // Pass to bottom
	 //Store U (2VelComp) values
	 for (i = il - 1; i < ir + 1; i++) {
	 bufSend[ctr] = U[i][jb];
	 ctr++;
	 }
	 //Store V (1VelComp) values
	 for (i = il; i < ir + 1; i++) {
	 bufSend[ctr] = V[i][jb];
	 ctr++;
	 }
	 ctr = 0;
	 MPI_Sendrecv(bufSend, jBufTotal, MPI_DOUBLE, rank_t, 1, bufRecv,
	 jBufTotal,
	 MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD, status);
	 //Receive U (2VelComp) values
	 for (i = il - 1; i < ir + 1; i++) {
	 U[i][jt + 1] = bufRecv[ctr];
	 ctr++;
	 }
	 //Receive V (1VelComp) values
	 for (i = il; i < ir + 1; i++) {
	 V[i][jt + 1] = bufRecv[ctr];
	 ctr++;
	 }
	 ctr = 0;

	 // Pass to Top
	 //Store U (2VelComp) values
	 for (i = il - 1; i < ir + 1; i++) {
	 bufSend[ctr] = U[i][jt];
	 ctr++;
	 }
	 //Store V (1VelComp) values
	 for (i = il; i < ir + 1; i++) {
	 bufSend[ctr] = V[i][jt - 1];
	 ctr++;
	 }
	 ctr = 0;
	 MPI_Sendrecv(bufSend, jBufTotal, MPI_DOUBLE, rank_t, 1, bufRecv,
	 jBufTotal,
	 MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD, status);
	 //Receive U (2VelComp) values
	 for (i = il - 1; i < ir + 1; i++) {
	 U[i][jb - 1] = bufRecv[ctr];
	 ctr++;
	 }
	 //Receive V (1VelComp) values
	 for (i = il; i < ir + 1; i++) {
	 V[i][jb - 2] = bufRecv[ctr];
	 ctr++;
	 }
	 ctr = 0;

	 */

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
