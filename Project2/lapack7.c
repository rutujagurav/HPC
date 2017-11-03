#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
//#include"lapacke.h"
//#include"blas.h"

/*void kinverse(double *a, int n){
    int i,j;
    for(i=0;i<n;i++)
        {
        for(j=0;j<n;j++)
        {
            if(i==j)
            {
                a[i*n+j] = 1;
            }
            if(i>j)
            {
                a[i*n+j] *= (-1);
            }
            if(i<j)
                {
                    a[i*n+j] = 0;
                }

        }
    }
}*/

void mydgetrfBLOCK(double* A,double* tempv, int* pvt, int n, int B)
{
        int i,j,k,t;
        int ib,l,m,p,q,end;
        int maxind, temps;
        double max, msum;
        double *ll;

for(ib=0; ib<n; ib+=B)
{
    end = ib + B-1;

        for(i=ib; i<=end; i++)
        {
            //pivoting
            maxind = i; max = abs(A[i*n+i]);

            for(t=i+1; t<n; t++)
            {
                if(abs(A[t*n+i]) > max)
                {
                    maxind = t;
                    max = abs(A[t*n+i]);
                }
            }

            //at this point I have the maximum number in the column and its position in the column

            if(max==0.0)
            {
                printf("LUfactoration failed: coefficient matrix is singular.");
                return;
            }

            else
            {
                if(maxind != i)
                {
                    //means we need to pivot
                    //printf("Hello");
                    temps=pvt[i];
                    pvt[i]=pvt[maxind];
                    pvt[maxind]=temps;

                    for(j=0; j<n; j++)
                    {
                        tempv[j] = A[i*n+j];
                        A[i*n+j] = A[maxind*n+j];
                        A[maxind*n+j] = tempv[j];
                    }
                }
            }

            //blas 2 factorization
            for(j=i+1; j<n; j++)
            {
                A[j*n+i] = A[j*n+i] / A[i*n+i];
                for(k=i+1; k<=end; k++)
                    A[j*n+k] = A[j*n+k] - A[j*n+i] * A[i*n+k];
            }

            //
            ll = (double*)calloc(sizeof(double), B*B);
            p=0;q=0;
            for(l=ib; l<=end; l++){
                for(m=ib; m<=end; m++){
                    if(l>m)
                    {
                        ll[p*B+q] = A[l*B+m] * (-1);
                    }
                    else if(l==m)
                    {
                        ll[p*B+q] = 1;
                    }
                    else
                    {
                        ll[p*B+q] = 0;
                    }

                    }
                    q++;
                }
                p++;
                q=0;
            }
            //kinverse(ll,B);
            //matInverse(ll,B);
             //printf("************");
            //llinverse(ll,B);
            //printf("************");
            //pink code line
            p=0;q=0;
            for(j=ib;j<=end;j++)
            {
                for(k=end+1;k<n;k++)
                {
                    msum = 0.0;
                    for(m=ib;m<=end;m++)
                    {
                        msum += ll[p*B+q] * A[m*n+k];
                        q++;
                    }
                A[j*n+k] = msum;
                q=0;
                }
                p++;
                q=0;
            }

            for(j=end+1;j<n;j++)
            {
                for(k=end+1;k<n;k++)
                {
                    double gmat = 0.0;
                    for(l=ib;l<=end;l++)
                    {
                        gmat += A[j*n+l] * A[l*n+k];
                    }
                    A[j*n+k] -= gmat;
                }
            }

        }
        //return A;
    }

double* genRand(double *z, int rows, int columns, double lower, double upper)
{

    int i;

    for(i=0; i<rows*columns; i++)
    {
        z[i] = ((double)rand()/RAND_MAX)*(upper-lower) + lower; //double precision floating point random numbers between [x,y]
        //printf("%f ", z[i]);

    }
    //printf("\nDone printing...\n");
    return z;
}

void FillMatrix(double* z, int rows, int columns)
{
    //printf("\n");
    double lower = 0.0;
    double upper = 10.0;
    z = genRand(z, rows, columns, lower, upper);

}

void mydgetrf(double* A,double* tempv, int* pvt, int n)
{
        int i,j,k,t;
        int maxind, temps;
        double max;
        for(i=0; i<n; i++)
        {
            //pivoting
            maxind = i; max = abs(A[i*n+i]);

            for(t=i+1; t<n; t++)
            {
                if(abs(A[t*n+i]) > max)
                {
                    maxind = t;
                    max = abs(A[t*n+i]);
                }
            }

            //at this point I have the maximum number in the column and its position in the column

            if(max==0.0)
            {
                printf("LUfactoration failed: coefficient matrix is singular.");
                return;
            }

            else
            {
                if(maxind != i)
                {
                    //means we need to pivot
                    //printf("Hello");
                    temps=pvt[i];
                    pvt[i]=pvt[maxind];
                    pvt[maxind]=temps;

                    for(j=0; j<n; j++)
                    {
                        tempv[j] = A[i*n+j];
                        A[i*n+j] = A[maxind*n+j];
                        A[maxind*n+j] = tempv[j];
                    }
                }
            }

            //for(i=0;i<n;i++)
                //printf("%d ",pvt[i]);

            //factorization
            for(j=i+1; j<n; j++)
            {
                A[j*n+i] = A[j*n+i] / A[i*n+i];
                for(k=i+1; k<n; k++)
                    A[j*n+k] = A[j*n+k] - A[j*n+i] * A[i*n+k];
            }

        }
        //return A;
}

void  mydtrsmF(double* x, double* y, double* A, double* B, int* pvt, int N)
{
    //printf("Inside mydtrsmF...");
    //printf("\n ##### %f %d",y[0],pvt[0]);
    y[0] = B[pvt[0]];
    //printf("\n ##### %d ",y[0]);
    int i,j,k;
    //int t = 1;
    //printf("Inside loop...");
    //double* sum = (double *) calloc(sizeof(double), t);
    //int NN = *N; // Problem here, runtime error

    //printf("%d ",NN);
    double sum;
    //printf("%f ",sum);
    //printf("Outside loop...");
    for(i=1; i<N; i++)
    {
                sum = 0.0;
        //printf("Inside loop...");
        for(j=0;j<i;j++)
              sum  = sum + y[j] * A[i*N+j];
        y[i] = B[pvt[i]] - sum;
    }
    //return y;
}

void  mydtrsmB(double* x, double* y, double* A, double* B, int* pvt, int N)
{
    //int NN  = *N;

    x[N-1] = y[N-1] / A[(N-1)*N+(N-1)];
    int i,j,k;
    double sum;
    for(i=N-2; i>=0; i--)
    {
                sum = 0.0;
        for(j=i+1;j<N;j++)
            sum  = sum + x[j] * A[i*N+j];
        x[i] = (y[i] - sum) / A[i*N+i];
    }
    //return x;
}

double* transpose(double* Mat, int n)
{
    int i,j;
    double temp;
    double *trans = (double *) calloc(sizeof(double), n*n);

    /*printf("\nMat---\n");
    for(i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
            printf("%f ", Mat[i*n+j]);
        printf("\n");
    }*/
    for(i=0; i<n; i++)
        for(j=0; j<n; j++)
        {
            trans[j*n+i] = Mat[i*n+j];
        }
    /*printf("\nTranspose---\n");
    for(i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
            printf("%f ", trans[i*n+j]);
        printf("\n");
    }*/

    return trans;
}

void checkCorrectness(double* B, double* x, int n)
{
        int i;
        double diff;
        double error = 0.0;
        for(i=0; i<n; i++)
        {
            diff = abs(B[i]-x[i]);
            //printf("diff = %f",diff);
            if(diff>error)
                error = diff;
        }

        printf("Error = %f \n", error);
        if(error < 0.0001)
            printf("Error < 1e-3 \n");
}



int main()
{

    int n[] = {4,1000,2000,3000,4000,5000};
    int i,j,k;
        srand((unsigned)time(NULL)); //setting seed to system time so that it gives different random numbers each time.


        double *A, *A2, *B, *B2, *tempv, *x, *y, *x1, *y1;
    int *pvt;

    double *z = (double *) calloc(sizeof(double), 1);
    z = genRand(z, 1, 1, 1.0, 10.0); //wasted float value

for(k=0; k<1; k++)
{
        printf("---------------------------------------------------------\n");
        printf("Number of Equations: %d \n",n[k]);
        //clock_t start, end;
    double cpu_time_used, gflops, error = 0.0, total_time = 0.0;
    double diff;

        A = (double *) calloc(sizeof(double), n[k]*n[k]);
        B = (double *) calloc(sizeof(double), n[k]*1);
        //Abk = (double *) calloc(sizeof(double), n[k]*n[k]);
        tempv = (double *) calloc(sizeof(double), n[k]*1);
        pvt = (int *) calloc(sizeof(int), n[k]*1);
        y = (double *) calloc(sizeof(double), n[k]*1);
        x = (double *) calloc(sizeof(double), n[k]*1);
        y1 = (double *) calloc(sizeof(double), n[k]*1);
        x1 = (double *) calloc(sizeof(double), n[k]*1);

        A2 = (double *) calloc(sizeof(double), n[k]*n[k]);
        B2 = (double *) calloc(sizeof(double), n[k]*1);

        //A3 = (double *) calloc(sizeof(double), n[k]*n[k]);
        //B3 = (double *) calloc(sizeof(double), n[k]*1);

        for(i=0;i<n[k];i++)
        pvt[i] = i;

    FillMatrix(A,n[k],n[k]);
        for(i=0; i<n[k]*n[k]; i++)
        {
                A2[i] = A[i];
                //A3[i] = A[i];
        }
        //FillMatrix(A2,n[k],n[k]);
    FillMatrix(B,n[k],1);
        //FillMatrix(B2,n[k],1);
        for(i=0; i<n[k]; i++)
        {
                B2[i] = B[i];
                //B3[i] = B[i];
        }
        A = transpose(A,n[k]);

        // note, to understand this part take a look in the MAN pages, at section of parameters.
    char    TRANS = 'N';
    int     INFO=n[k];
    int     LDA = n[k];
    int     LDB = n[k];
    int     N = n[k];
    int     NRHS = 1;
    int     IPIV[n[k]] ;


        //struct timespec tstart={0,0}, tend={0,0};

//Library Lapack -------------------------------------------------------
/*        clock_gettime(CLOCK_MONOTONIC, &tstart);
        //start = clock();
                // LU factorization
                LAPACK_dgetrf(&N,&N,A,&LDA,IPIV,&INFO);
        clock_gettime(CLOCK_MONOTONIC, &tend);
        cpu_time_used = ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec);
        //end = clock();
    //cpu_time_used = ((double) (end - start));
        gflops = (2*pow(n[k],3)) /(3*cpu_time_used * pow(10,9));

        total_time = total_time + cpu_time_used;
        printf("Part1.1: Library Functions\n");
        printf("\n Solved LU Factorization --- Exec. time = %.16f \n", cpu_time_used / CLOCKS_PER_SEC);
    printf("\n Solved LU Factorization --- Gflops = %.16f \n", gflops);

    char     SIDE = 'L';
    char     UPLO = 'L';
    char     DIAG = 'U';
    int      M    = 1;
    double   a    = 1.0;


           // change the order of B according to IPIV[] from LU factorization
                double tmp;

                for(i = 0; i < N; i++)
                {
                        tmp = B[IPIV[i]-1];
                        B[IPIV[i]-1] = B[i];
                        B[i] = tmp;
                }

                // forward  L(Ux) = B => y = Ux

        //start = clock();
        clock_gettime(CLOCK_MONOTONIC, &tstart);
                dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&N,&M,&a,A, &N, B, &N);
        clock_gettime(CLOCK_MONOTONIC, &tend);
        cpu_time_used = ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec);
        //end = clock();
    //cpu_time_used = ((double) (end - start));
        gflops = ((double)pow(n[k],2)) /(double)(cpu_time_used * pow(10,9));
        total_time = total_time + cpu_time_used;

        //printf("\n Solved Forward Substitution --- Exec. time = %.16f \n", cpu_time_used / CLOCKS_PER_SEC);
    //printf("\n Solved Forward Substitution --- Gflops = %.16f \n", gflops);


                UPLO = 'U';
                DIAG = 'N';

                // backward Ux = y

        //start = clock();
        clock_gettime(CLOCK_MONOTONIC, &tstart);
                dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&N,&M,&a,A, &N, B, &N);
        clock_gettime(CLOCK_MONOTONIC, &tend);
        cpu_time_used = ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec);
        //end = clock();
    //cpu_time_used = ((double) (end - start));
        gflops = ((double)pow(n[k],2)) /(double)(cpu_time_used * pow(10,9));
        total_time = total_time + cpu_time_used;

        //printf("\n Solved Backward Substitution --- Exec. time = %.16f \n", cpu_time_used / CLOCKS_PER_SEC);
    //printf("\n Solved Backward Substitution --- Gflops = %.16f \n", gflops);

        //printf("Total time of execution = %.16f \n",total_time);

        //printf("print the result from library : {");
    //for (i=0;i<N;i++)
    //{
      //          printf("%f ",B[i]);
    //}
    //printf("}");

printf("\n");

*/
//My lapack ----------------------------------------------------------

printf("\n");
printf("Part1.2: My GEPP Functions\n");

//clock_gettime(CLOCK_MONOTONIC, &tstart);
        mydgetrf(A,tempv,pvt,n[k]);

        for(i=0;i<n[k];i++)
        {
            for(j=0;j<n[k];j++)
                printf("%f ",A[i*n[k]+j]);
            printf("\n");
        }

//clock_gettime(CLOCK_MONOTONIC, &tend);
//cpu_time_used = ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec);
//gflops = (2*pow(n[k],3)) /(3*cpu_time_used * pow(10,9));
//printf("\n Solved LU Factorization --- Exec. time = %.16f \n", cpu_time_used / CLOCKS_PER_SEC);
//printf("\n Solved LU Factorization --- Gflops = %.16f \n", gflops);

        mydtrsmF(x,y,A,B,pvt,n[k]);
        printf("----- \n");
        for(i=0;i<n[k];i++)
            printf("%f ",y[i]);

        mydtrsmB(x,y,A,B,pvt,n[k]);


        printf("\n My Result = { \n");
        for(i=0; i<n[k]; i++)
            printf("%f ",x[i]);
    printf("\n}");

//printf("\n");

        //checkCorrectness(B,x,n[k]);

//free(A); free(A2);
//free(B); free(B2);
//free(tempv);
//free(pvt);
//free(x);
//free(y);
//free(IPIV);

//Blocked lapack ----------------------------------------------------------

printf("\n");
printf("Part2: Blocked GEPP Functions\n");

//clock_gettime(CLOCK_MONOTONIC, &tstart);
        mydgetrfBLOCK(A2,tempv,pvt,n[k],2);
        for(i=0;i<n[k];i++)
        {
            for(j=0;j<n[k];j++)
                printf("%f ",A[i*n[k]+j]);
            printf("\n");
        }
//clock_gettime(CLOCK_MONOTONIC, &tend);
//cpu_time_used = ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec);
//gflops = (2*pow(n[k],3)) /(3*cpu_time_used * pow(10,9));
//printf("\n Solved LU Factorization --- Exec. time = %.16f \n", cpu_time_used / CLOCKS_PER_SEC);
//printf("\n Solved LU Factorization --- Gflops = %.16f \n", gflops);


        mydtrsmF(x1,y1,A2,B2,pvt,n[k]);
        printf("----- \n");
        for(i=0;i<n[k];i++)
            printf("%f ",y[i]);

        mydtrsmB(x1,y1,A2,B2,pvt,n[k]);


        printf("\n My Result = { \n");
        for(i=0; i<n[k]; i++)
            printf("%f ",x[i]);
    printf("\n}");
	checkCorrectness(x1,x,n[k]);

}
    return 0;
}
