#include "helper.h"
#include "init.h"

int read_parameters(    const char *szFileName,
               //     char *problem,
                 //   char *geometry,
                    int  *imax,
                    int  *jmax,
                    double *xlength,
                    double *ylength,
                    double *dt,
                    double *t_end,
                    double *tau,
                    double *dt_value,
                    double *eps,
                    double *omg,
                    double *alpha,
                    int  *itermax,
                    double *GX,
                    double *GY,
                    double *Re,
                    double *Pr,
                    double *UI,
                    double *VI,
                    double *PI,
                    double *TI,
                    double *TH,
                    double *TC,
                    double *beta,
                    double *dx,
                    double *dy
                    )
{
   READ_DOUBLE( szFileName, *xlength );
   READ_DOUBLE( szFileName, *ylength );

   READ_DOUBLE( szFileName, *Re    );
   READ_DOUBLE( szFileName, *t_end );
   READ_DOUBLE( szFileName, *dt    );

   READ_INT   ( szFileName, *imax );
   READ_INT   ( szFileName, *jmax );

   READ_DOUBLE( szFileName, *omg   );
   READ_DOUBLE( szFileName, *eps   );
   READ_DOUBLE( szFileName, *tau   );
   READ_DOUBLE( szFileName, *alpha );

   READ_INT   ( szFileName, *itermax );
   READ_DOUBLE( szFileName, *dt_value );

   READ_DOUBLE( szFileName, *UI );
   READ_DOUBLE( szFileName, *VI );
   READ_DOUBLE( szFileName, *GX );
   READ_DOUBLE( szFileName, *GY );
   READ_DOUBLE( szFileName, *PI );

   *dx = *xlength / (double)(*imax);
   *dy = *ylength / (double)(*jmax);

   return 1;
}

void init_uvp (double UI, double VI, double PI, int imax, int jmax, double **U, double **V, double **P, int **flag)
{
    for(int i=0; i<imax; i++){
        for(int j=0; j<jmax; j++){
            if(flag[i][j]&(1<<0)){
                
                U[i][j] = UI;
                V[i][j] = VI;
                P[i][j] = PI;
            }
        }
    }
}

void init_uvpt (double UI, double VI, double PI, double TI, int imax, int jmax, double **U, double **V, double **P, double **T, int **flag)
{
    for(int i=0; i<imax; i++){
        for(int j=0; j<jmax; j++){
            if(flag[i][j]&(1<<0)){
                
                U[i][j] = UI;
                V[i][j] = VI;
                P[i][j] = PI;
                T[i][j] = TI;
            }
        }
    }
}

int  isfluid(int pic){
    if((pic == 2)||(pic == 3)||(pic == 4)) {return 1;}
    else {return 0;}
}

void assert_error()
{
    char szBuff[80];
    sprintf( szBuff, "Forbidden geometry. \n");
    ERROR( szBuff );
}
int forbidden_WO(int **pic, int i, int j)
{
    if( (pic[i-1][j]==4)&&(pic[i+1][j]==4) ) {return 1;}
    else {return 0;}
}

// fluid on opposite sides top and bottom
int forbidden_NS(int **pic, int i, int j)
{
    if( (pic[i][j+1]==4)&&(pic[i][j-1]==4) ) {return 1;}
    else {return 0;}
}

//Avoids any forbidden configuration
void assert_forbid(int imax, int jmax, int **pic)
{
    //inner obstacles
    for(int i=1; i<imax-1; i++)
    {
        for(int j=1; j<jmax-1; j++)
        {
            if(pic[i][j]!=4)
            {
                if(forbidden_WO(pic,i,j)||forbidden_NS(pic,i,j))
                {
                    assert_error();
                }
            }
        }
    }
    //left boundary
    for(int j=1; j<jmax-1; j++)
    {
        if(pic[0][j]!=4)
        {
            if(forbidden_NS(pic,0,j))
            {
                assert_error();
            }
        }
    }
    //right boundary
    for(int j=1; j<jmax-1; j++)
    {
        if(pic[imax-1][j]!=4)
        {
            if(forbidden_NS(pic,imax-1,j))
            {
                assert_error();
            }
        }
    }
    //top boundary
    for(int i=1; i<imax-1; i++)
    {
        if(pic[i][jmax-1]!=4)
        {
            if(forbidden_WO(pic,i,jmax-1))
            {
                assert_error();
            }
        }
    }
    //bottom boundary
    for(int i=1; i<imax-1; i++)
    {
        if(pic[i][0]!=4)
        {
            if(forbidden_WO(pic,i,0))
            {
                assert_error();
            }
        }
    }
    
}

void init_flag(char* problem, char* geometry, int imax, int jmax, int **flag)
{
    int **pic = read_pgm(geometry);
    
    for (int i=0; i<imax; i++)
    {
        for (int j=0; j<jmax; j++)
        {
            
            flag[i][j] = 0;
            
            assert_forbid(imax, jmax, pic);

            switch(pic[i][j])
            {
                case 0: //fluid
                    flag[i][j] = 1<<1;
                    break;
                    
                case 1: //no-slip
                    flag[i][j] = 1<<2;
                    break;
                    
                case 2: //free-slip
                    flag[i][j] = 1<<3;
                    break;
                    
                case 3: //outflow
                    flag[i][j] = 1<<4;
                    break;
                    
                case 4: //inflow
                    flag[i][j] = 1<<0;
                    break;
            }
            
            if(!isfluid(pic[i][j])) //set boundaries if not
            {
                
                if(i<imax-1 && pic[i+1][j]==4)
                {
                    flag[i][j] |= 1<<8;
                }
                if( i>0 && pic[i-1][j]==4)
                {
                    flag[i][j] |= 1<<7;
                }
                if(j<jmax-1 && pic[i][j+1]==4)
                {
                    flag[i][j] |= 1<<5;
                }
                if(j>0 && pic[i][j-1]==4)
                {
                    flag[i][j] |= 1<<6;
                }
            }
            
            
        }
        
    }
    
    printf("PROGRESS: flags are set \n");
    
}
