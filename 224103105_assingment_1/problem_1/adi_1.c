#include<stdio.h>
#include<stdlib.h>
#include<math.h>
# define m 31
# define n 21
int main()
{
    double tem[m][n], psi_new[m][n],P[m-2],Q[m-2],d[m-2],d1[n-2],P1[n-2],Q1[n-2],dx=6.0/(m-1),dy=4.0/(n-1),as,aw,ap,ae,an,aj,bj,cj,ai,bi,ci ,error;
    int i,j, iteration = 0;
    as=(1.0/pow(dy,2.0));
    aw=(1.0/pow(dx,2.0));
    ap=2.0*((1.0/pow(dx,2.0))+(1.0/pow(dy,2.0)));
    ae=(1.0/pow(dx,2.0));
    an=(1.0/pow(dy,2.0));
    ai=ap;
    bi=-ae;
    ci=-aw;
    aj=ap;
    bj=-an;
    cj=-as;

    for(i=0; i<m ; i++)
    {
        for(j=0 ; j<n; j++)
        {
            if(j==0)
            {
                //bottom boundary
                if(i<=5){
                    psi_new[i][j]=0;
                }
                if(i>= 6){
                psi_new[i][j]=100.0;
                }
            }
            //top boundary
            else if(j == n-1 )
            {
                psi_new[i][j]=0.0;
            }
            //left boundary
            else if(i == 0){
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



          for(i=0; i<m; i++)
          {
         for(j=0; j<n; j++)
         {
          tem[i][j]= psi_new[i][j];
        }
}
   for(j=1; j<n-1; j++)
   {
    d[0]=an*(psi_new[0][j+1]+psi_new[0][j-1]);
    P[0]=-bi/ai;
    Q[0]=d[0]/ai;

    for(i=1; i<m-1; i++)
    {
        d[i]=an*(psi_new[i][j+1]+psi_new[i][j-1]);
        P[i]=-(bi/(ai+ci*P[i-1]));
        Q[i]=(d[i]-ci*Q[i-1])/(ai+ci*P[i-1]);
    }
    for(i=m-2; i>0; i--)
    {

         psi_new[i][j]=P[i]*psi_new[i+1][j]+Q[i];

    }
    }
    for(i=1; i<m-1; i++)
    {
        d1[0]=ae*(psi_new[i+1][0]+psi_new[i-1][0]);
        P1[0]=-(bj/aj);
        Q1[0]=(d1[0])/aj;

    for(j=1; j<n-1; j++)
    {
        d1[j]=ae*(psi_new[i+1][j]+psi_new[i-1][j]);     
        P1[j]=-(bj/(aj+cj*P1[j-1]));
        Q1[j]=(d1[j]-cj*Q1[j-1])/(aj+cj*P1[j-1]);
    }
    for(j=n-2; j>0; j--){

    psi_new[i][j]=P1[j]*psi_new[i][j+1]+Q1[j];
    }
    }

    for (j=0; j<n; j++)
    {
    psi_new[m-1][j]=psi_new[m-2][j] ;
    }
    error=0.0;
       for(i=0; i<m; i++)
       {
            for(j=0; j<n; j++){
            error+=pow((psi_new[i][j]-tem[i][j]),2);
        }
    }
      iteration=iteration+1 ;
      error=sqrt(error/((m-2)*(n-2)));
      printf("iteration %d\t",iteration);
      printf("error %.9f\n",error);
      fprintf(file1,"%d\t%lf\n",iteration,error);
      }
 FILE *file2;
file2=fopen("adi_1.plt","w");
fprintf(file2,"VARIABLES= \"x\",\"y\",\"PHI\"\n");
fprintf(file2, "ZONE t=\"BLOCK1\", J=31,I=21,F= POINT \n\n");
for(i = 0; i < m; i++)  {
        for(j = 0; j < n; j++)
{
    fprintf(file2,"%lf\t%lf\t%lf\n",dx*i,dy*j,psi_new[i][j]);
     printf("%f\t",psi_new[i][j]);

 }
} 
return 0;

}

