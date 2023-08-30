#include<stdio.h>
#include<stdlib.h>
#include<math.h>
# define m 31
# define n 21
int main(){
   double psi_new[n][m],B,tem[n][m],ai,bi,ci,d[m-2],P[m-2],Q[m-2],dy=4.0/(n-1), dx=6.0/(m-1), error;
   int i,j,iteration=0;
    B = (dx/dy);
    ai=-2*(1+pow(B,2));
    bi=1;
    ci=1;
     for(i=0; i<n ; i++)
    {
        for(j=0 ; j<m; j++)
        {
            if(i==0)
            {
                //bottom boundary
                if(j<=5){
                    psi_new[i][j]=0;
                }
                if(j >= 6){
                psi_new[i][j]=100.0;
                }
            }
            //top boundary
            else if(i == n-1 )
            {
                psi_new[i][j]=0.0;
            }
            //left boundary
            else if(j == 0){
                psi_new[i][j]=0.0;
                }
            else{
                psi_new[i][j]=0.0;
                }
        }
                
    }
    FILE *file1;
    file1=fopen("error.txt","w");
    error=1;
   while(error > 1e-6)
   {
     error = 0.0;
        for(j=1; j<m-1;j ++)
        {
        d[0]=-(pow(B,2))*(psi_new[0][j+1]+psi_new[0][j-1]);
        P[0]=-bi/ai;
        Q[0]=d[0]/ai;

        for(i=1; i<n-1; i++){
            d[i]=-(pow(B,2))*(psi_new[i][j+1]+psi_new[i][j-1]);;

                   P[i]=-(bi/(ai+ci*P[i-1]));
                   Q[i]=(d[i]-ci*Q[i-1])/(ai+ci*P[i-1]);
        }
         for(i=n-2; i>0; i--){
            tem[i][j]=psi_new[i][j];
                psi_new[i][j]=P[i]*psi_new[i+1][j]+Q[i];
              error+=pow((psi_new[i][j]- tem[i][j]),2.0);
        }

    }

      for (i=0; i<n; i++){
        psi_new[i][m-1]=psi_new[i][m-2];
      }

      error=sqrt(error/((m-2)*(n-2)));
      iteration = iteration+1;
      printf("iteration %d\t",iteration);
      printf("error %.9f\n",error);
      fprintf(file1,"%d\t%.9lf\n",iteration,error);
   }

     FILE *file2;
   file2=fopen("tdmam1.plt","w");
    fprintf(file2,"VARIABLES= \"x\",\"y\",\"PHI\"\n");
    fprintf(file2, "ZONE t=\"BLOCK1\", J=31,I=21,F= POINT \n\n");

    for(int i=0; i<n; i++)
    {
            for(int j=0 ; j<m; j++)
            {
               fprintf(file2, "%lf\t%lf\t%lf\n",j*dy,i*dx,psi_new[i][j]);
                printf("%f\t",psi_new[i][j]);
            }
    }

    return 0;

}