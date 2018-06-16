#include "precice_adapter.h"
#include "boundary_val.c"


int *precice_set_interface_vertices(int imax, int jmax, double dx, double dy, double x_origin, double y_origin, int num_coupling_cells, int meshID, int **flag)
{
    int dimension =3;
    int coupledcell=0;
    int* vertexIDs  = (int*)malloc(num_coupling_cells*sizeof(int));
    double* vertices = (double*) malloc(sizeof(double) * num_coupling_cells * dimension);
    
    for (int j=1; j<=jmax; j++){
        if(flag[imax+1][j]&1<<9 && flag[imax+1][j]&1<<8){      //East boundary
        vertices[dimension*coupledcell]     = x_origin + (imax - 1)*dx;
        vertices[dimension*coupledcell + 1] = y_origin + (j - 0.5)*dy;
        vertices[dimension*coupledcell + 2] = 0;
        coupledcell++;
    }
    }
    for (int j=1; j<=jmax; j++){
        if(flag[0][j]&1<<9 && flag[0][j]&1<<7){     // West boundary
            vertices[dimension*coupledcell]     = 0;
            vertices[dimension*coupledcell + 1] = y_origin + (j - 0.5)*dy;
        vertices[dimension*coupledcell + 2] = 0;
        coupledcell++;
    }
    }
    for (int i=1; i<=imax; i++){
        if(flag[i][0]&1<<9 && flag[i][0]&1<<6){    // Bottom boundary
        vertices[dimension*coupledcell] = x_origin + (i - 0.5)*dx;
            vertices[dimension*coupledcell + 1]     = 0;
        vertices[dimension*coupledcell + 2] = 0;
        coupledcell++;
    }
    }
    for (int i=1; i<=imax; i++){
        if(flag[i][jmax+1]&1<<9 && flag[i][jmax+1]&1<<5){   // Top boundary
            vertices[dimension*coupledcell] = x_origin + (i - 0.5)*dx;
        vertices[dimension*coupledcell + 1]     = y_origin + (jmax-1)*dy;
        vertices[dimension*coupledcell + 2] = 0;

        coupledcell++;
    }
}
            precicec_setMeshVertices(meshID, num_coupling_cells, vertices, vertexIDs);
    return vertexIDs;
}


void precice_write_temperature(int imax, int jmax, int num_coupling_cells, double *temperature, int *vertexIDs, int temperatureID, double **T, int **flag)
{
    int count=0;
    for (int i=1;i<=imax;i++)
    {
            if(flag[i][jmax+1]&1<<9 && flag[i][jmax+1]&1<<5){
             temperature[count]=T[i][jmax];
                count++;
            }
            if(flag[i][0]&1<<9 && flag[i][0]&1<<6){
                temperature[count]=T[i][1];
                count++;
        }
}
        for (int j=1; j<=jmax; j++) {
        if(flag[imax+1][j]&1<<9 && flag[imax+1][j]&1<<8){
                temperature[count]=T[imax][j];
            count++;
            }
            if(flag[0][j]&1<<9 && flag[0][j]&1<<7){
                temperature[count]=T[1][j];
                count++;
            }
        }
        precicec_writeBlockScalarData(temperatureID, num_coupling_cells, vertexIDs, temperature);
}


void write_checkpoint(double time, double **U, double **V, double **T, double *time_cp, double **U_cp, double **V_cp,
                      double **T_cp, int imax, int jmax)
{
    time_cp = time;
    
    for (int i=1; i<=imax; i++) {
        for (int j=1; j<=jmax; j++) {
            T_cp[i][j] = T[i][j];
            U_cp[i][j] = U[i][j];
            V_cp[i][j] = V[i][j];
        }
    }
}


void restore_checkpoint(double *time, double **U, double **V, double **T, double time_cp, double **U_cp,
                        double **V_cp, double **T_cp, int imax, int jmax)
{
    
    time = time_cp;
    
    for (int i=1; i<=imax; i++) {
        for (int j=1; j<=jmax; j++) {
            T[i][j] = T_cp[i][j];
            U[i][j] = U_cp[i][j];
            V[i][j] = V_cp[i][j];
        }
    }
}

