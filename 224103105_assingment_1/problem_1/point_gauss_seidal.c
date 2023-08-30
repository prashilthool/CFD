#include<stdio.h>
#include<stdlib.h>
#include<math.h>
# define m 31
# define n 21

int main()
{
     double dx = 6.0/(m-1),dy = 4.0/(n-1),psi_old[n][m],psi_new[n][m],as,aw,ap,ae,an, Error,tem[n][m];
    int iteration = 0,i,j;
    Error=1;
    ap=2.0*((1.0/pow(dx,2.0))+(1.0/pow(dy,2.0)));
    as=(1.0/pow(dy,2.0));
    aw=(1.0/pow(dx,2.0));
    ae=(1.0/pow(dx,2.0));
    an=(1.0/pow(dy,2.0)); 
    for(int i=0; i<n ; i++)
        {
        for(int j=0 ; j<m ; j++)
        {
            if(i==0)
            {
                if(j<=4){
                    psi_new[i][j]=0;
                }
            else if(j >= 5){
                psi_new[i][j]=100.0;
                }
            }
            else if(i == n-1 ){
                psi_new[i][j]=0.0;
                }
            else if(j == 0){
                psi_new[i][j]=0.0;
                }
            else{
                psi_new[i][j]=0.0 ;
                  }
                }
    }
    //point gauss-seidal method
    FILE *file1;
    file1= fopen("error.txt","w");
    fprintf(file1,"Iter\t Error\n");
     while(Error>1e-6)
    {
        Error=0;
        for( i=0; i<n ; i++)
        {
            for( j=0 ; j<m ; j++)
            {
                tem[i][j]=psi_new[i][j];
            }
            }

        for( i=1; i<n-1 ; i++)
        {
            for(j=1 ; j<m-1 ; j++)
            {
                psi_new[i][j]=(1/ap)*((ae*psi_new[i][j+1])+(aw*psi_new[i][j-1])+(an*psi_new[i+1][j])+(as*psi_new[i-1][j]));
                Error = Error + pow((psi_new[i][j]-tem[i][j]),2);

            }
        }
        for(i=0; i<n ; i++)
        {
            psi_new[i][m-1]=psi_new[i][m-2];
        }
       
                Error =sqrt(Error/(m*n));
                printf("iteration %d\t",iteration);
                printf("error %.10lf\n",Error);
                fprintf(file1, "%d\t%.10f\n",iteration,Error);
                iteration++;
    }
    FILE *file2;
   file2=fopen("Stream1b.plt","w");
    fprintf(file2,"VARIABLES= \"x\",\"y\",\"PHI\"\n");
    fprintf(file2, "ZONE t=\"BLOCK1\", J=31,I=21,F= POINT \n\n");

    for(int i=0; i<n ; i++)
    {
            for(int j=0 ; j<m ; j++)
            {
               fprintf(file2, "%lf\t%lf\t%lf\n",j*dy,i*dx,psi_new[i][j]);
                printf("%f\t",psi_new[i][j]);
            }
    }

    return 0;

}



