#include "precice_adapter.h"
#include "boundary_val.c"


int *precice_set_interface_vertices(int imax, int jmax, double dx, double dy, double x_origin, double y_origin, int num_coupling_cells, int meshID, int **FLAG)
{
    int dimension =3;
    int* vertexIDs  = (int*)malloc(num_coupling_cells*sizeof(int));
    double* vertices = (double*) malloc(sizeof(double) * num_coupling_cells * dimension);
    
    for (int j=1; j<=jmax; j++){
        if(flag[i][j]&1<<9 && flag[i][j]&1<<8){      //East boundary
        vertices[dimension*coupledcell]     = 0;
        vertices[dimension*coupledcell + 1] = y_origin + (j - 0.5)*dy;
        vertices[dimension*coupledcell + 2] = 0;
        coupledcell++;
    }
    }
    for (int j=1; j<=jmax; j++){
        if(flag[i][j]&1<<9 && flag[i][j]&1<<7){     // West boundary
            vertices[dimension*coupledcell]     = x_origin + (imax - 1)*dx;
            vertices[dimension*coupledcell + 1] = y_origin + (j - 0.5)*dy;
        vertices[dimension*coupledcell + 2] = 0;
        coupledcell++;
    }
    }
    for (int i=1; i<=imax; i++){
        if(flag[i][j]&1<<9 && flag[i][j]&1<<6){    // Bottom boundary
        vertices[dimension*coupledcell] = x_origin + (i - 0.5)*dx;
            vertices[dimension*coupledcell + 1]     = 0;
        vertices[dimension*coupledcell + 2] = 0;
        coupledcell++;
    }
    }
    for (int i=1; i<=imax; i++){
        if(flag[i][j]&1<<9 && flag[i][j]&1<<5){   // Top boundary
            vertices[dimension*coupledcell] = x_origin + (i - 0.5)*dx;
        vertices[dimension*coupledcell + 1]     = y_origin + (jmax-1)*dy;
        vertices[dimension*coupledcell + 2] = 0;

        coupledcell++;
    }
}
            precicec_setMeshVertices(meshID, num_coupling_cells, vertices, vertexIDs);
}


void precice_write_temperature(int imax, int jmax, int num_coupling_cells, double *temperature, int *vertexIDs, int temperatureID, double **T, int **flag)
{
    for (int i=1;i<=imax;i++)
    {
            if(flag[i][j]&1<<9 && flag[i][j]&1<<5){
             temperature[count]=T[i][jmax];
                count++;
            }
            if(flag[i][j]&1<<9 && flag[i][j]&1<<6){
                temperature[count]=T[i][0];
                count++;
        }
}
        for (int j=1; j<=jmax; j++) {
        if(flag[i][j]&1<<9 && flag[i][j]&1<<8){
                temperature[count]=T[imax][j];
            count++;
            }
            if(flag[i][j]&1<<9 && flag[i][j]&1<<7){
                temperature[count]=T[0][j];
                count++;
            }
        }
        precicec_writeBlockScalarData(temperatureID, num_coupling_cells, vertexIDs, temperature);
}

