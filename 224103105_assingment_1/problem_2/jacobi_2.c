#include<stdio.h>
#include<stdlib.h>
#include<math.h>
# define m 41
int main()
{
    double dx = 1.0/(m-1),dy = 1.0/(m-1),psi_old[m][m],psi_new[m][m],as,aw,ap,ae,an, Error;
     int iteration = 0,i,j;
    Error=1;
    ap=2.0*((1.0/pow(dx,2.0))+(1.0/pow(dy,2.0)));
    as=(1.0/pow(dy,2.0));
    aw=(1.0/pow(dx,2.0));
    ae=(1.0/pow(dx,2.0));
    an=(1.0/pow(dy,2.0)); 
    //boudary conditions
       
    for ( i = 0; i < m; i++)
    {
        for ( j = 0; j <m; j++)
        {
           
           psi_new[0][j]=1.0;
           psi_new[m-1][j]=1.0;
           psi_new[i][0]=1.0;
          if (j==m-1)
          {
             psi_new[i][j]=0.0;
          }
          else
          {
            psi_new[i][j]=0.0;
          } 
           
        }
        
    }
    FILE *file1;
    file1= fopen("error.txt","w");
    fprintf(file1,"Iter\t Error\n");
    while (Error>1e-6)
    {
        for ( i = 0; i <m; i++)
        {
            for ( j = 0; j <m; j++)
            {
                psi_old[i][j]=psi_new[i][j];
            }
            
        }
        for ( i = 1; i <m-1; i++)
        {
            for ( j = 0; j <m-1; j++)
            {
                psi_new[i][j]=(1/ap)*((ae*psi_old[i][j+1])+(aw*psi_old[i][j-1])+(an*psi_old[i+1][j])+(as*psi_old[i-1][j]));

            }
            
        }
        Error=0.0;
        for ( i = 0; i < m; i++)
        {
            for ( j =0; j< m; j++)
            {
                Error = Error + pow((psi_new[i][j]-psi_old[i][j]),2);
            }
            
        }
        Error=sqrt(Error/((m-2)*(m-2)));
        iteration=iteration+1;
        printf("Iteration= %d\t", iteration);
        printf("Error= %.9f\n", Error);
        fprintf(file1, "%d \t %.9f \n", iteration, Error);
        
        
        
    }
    
    FILE *file2;
    file2=fopen("Stream1Ta.plt","w");
    fprintf(file2,"VARIABLES= \"X\",\"Y\",\"T\"\n");
    fprintf(file2, "ZONE t=\"BLOCK1\", J=41,I=41,F= POINT \n\n");

    for( i=0; i<m ; i++)
    {
            for(j=0 ; j<m ; j++)
            {
                fprintf(file2, "%lf \t %lf \t %lf \n",j*dy,i*dx,psi_new[i][j]);
                printf("%f\t",psi_new[i][j]);
            }
    }
    
    
    
    return 0;
}
