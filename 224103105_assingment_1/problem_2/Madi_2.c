#include<stdio.h>
#include<stdlib.h>
#include<math.h>
# define m 41
int main()
{
    
    double dx = 1.0/(m-1),dy = 1.0/(m-1),tem[m][m],psi_new[m][m],Error,B,ai,bi,ci,d[m],P[m],Q[m],d1[m],P1[m],Q1[m];
    int iteration = 0,i,j;
    Error=1.0;
    B=(dx/dy);
    ai=-2*(1+pow(B,2));
    bi=1.0;
    ci=1.0;
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
    fprintf(file1,"Iteration\t Error\n");
        while (Error>1e-6)
    {
        for ( i = 0; i < m; i++)
        {
           for (j  = 0; j <m; j++)
           {
            tem[i][j]=psi_new[i][j];
           }
           
        }
        
      
        for ( j = 0.5;  j< m-1; j++)
        {
          d[0]=-pow(B,2)*(psi_new[0][j+1]+psi_new[0][j-1]);
          P[0]=-bi/ai; 
          Q[0]=d[0]/ai;
        
        for (i = 0.5; i < m-1; i++)
        {
         
           d[i]=-pow(B,2)*(psi_new[i][j+1]+psi_new[i][j-1]);
           P[i]=-(bi/(ai+ci*P[i-1]));
           Q[i]=(d[i]-ci*Q[i-1])/(ai+ci*P[i-1]);
           
         
        }
          for ( i = m-2; i > 0; i--)
        {
           
            psi_new[i][j]=P[i]*psi_new[i+1][j]+Q[i];
            
        }
        }
        for ( i = 1; i < m-1; i++)
        {
          d1[0]=-pow(B,2)*(psi_new[i+1][0]+psi_new[i-1][0]);
          P1[0]=-bi/ai; 
          Q1[0]=d1[0]/ai;

          for ( j = 1; j < m-1; j++)
          {
           d1[j]=-pow(B,2)*(psi_new[i][j+1]+psi_new[i][j-1]);
           P1[j]=-(bi/(ai+ci*P1[j-1]));
           Q1[j]=(d1[j]-ci*Q1[j-1])/(ai+ci*P1[j-1]);
           
          }
          for ( j = m-2; j >0; j--)
          {
            psi_new[i][j]=P1[j]*psi_new[i+1][j]+Q1[j];
          }
          
          
        }
        Error=0.0;
        for ( i = 0; i <m; i++)
        {
            for ( j = 0; j < m; j++)
            {
               Error = Error + pow((psi_new[i][j]-tem[i][j]),2.0);
            }
            
        }
        Error=sqrt(Error/((m-2)*(m-2)));
        iteration=iteration+1;
        printf("Iteration= %d\t", iteration);
        printf("Error= %.9f\n", Error);
        fprintf(file1, "%d \t %.9f \n", iteration, Error);
        
        
    }
     FILE *file2;
    file2=fopen("Stream1Tddd.plt","w");
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
