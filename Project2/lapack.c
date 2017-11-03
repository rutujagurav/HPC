#include<stdio.h>
#include<stdlib.h>
#include <time.h>
#include"lapacke.h"
#include"blas.h"

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
    printf("\n");
    double lower = 0.0;
    double upper = 10.0;
    z = genRand(z, rows, columns, lower, upper);

}


int main()
{
    int n[] = {3,1000,2000,3000,4000,5000};
    int i,j,k;
	srand((unsigned)time(NULL)); //setting seed to system time so that it gives different random numbers each time.
    
	
	double *A, *B;

    double *z = (double *) calloc(sizeof(double), 1);
    z = genRand(z, 1, 1, 1.0, 10.0); //wasted float value
	
for(k=0; k<1; k++)
{
	clock_t start, end;
    double cpu_time_used, gflops, error = 0.0;
    double diff;
	
    A = (double *) calloc(sizeof(double), n[k]*n[k]);
    B = (double *) calloc(sizeof(double), n[k]*1);

    FillMatrix(A,n[k],n[k]);
    FillMatrix(B,n[k],1);

	// note, to understand this part take a look in the MAN pages, at section of parameters.
    char    TRANS = 'N';
    int     INFO  = n[k];
    int     LDA   = n[k];
    int     LDB   = n[k];
    int     N     = n[k];
    int     NRHS  = 1;
    int     IPIV[n[k]] ;

    
	
	
	
	start = clock();
		// LU factorization
		LAPACK_dgetrf(&N,&N,A,&LDA,IPIV,&INFO);
		

		// change the order of B according to IPIV[] from LU factorization
		double tmp;
		
		for(i = 0; i < N; i++)
		{
			tmp = B[IPIV[i]-1];
			B[IPIV[i]-1] = B[i];
			B[i] = tmp;
		}

		
	char     SIDE = 'L';
    char     UPLO = 'L';
    char     DIAG = 'U';
    int      M    = 1;
    double   a    = 1.0;
		
		// forward  L(Ux) = B => y = Ux
		dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&N,&M,&a,A, &N, B, &N);
		UPLO = 'U';
		DIAG = 'N';
		// backward Ux = y
		dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&N,&M,&a,A, &N, B, &N);
	
	end = clock();
    cpu_time_used = ((double) (end - start));
    gflops = ((double)2*pow(n[k],3) /(double)(cpu_time_used * pow(10,9)));
	
	printf("\n Solved Linear Equations System --- Exec. time = %.16f \n", cpu_time_used / CLOCKS_PER_SEC);
    printf("\n Solved Linear Equation System --- Gflops = %.16f \n", gflops);

	
	printf("print the result : {");
    for (i=0;i<N;i++)
    {
		printf("%f ",B[i]);
    }
    printf("}");



	/*printf("compute the LU factorization...\n");
    //void LAPACK_dgetrf( lapack_int* m, lapack_int* n, double* A, lapack_int* lda, lapack_int* ipiv, lapack_int *info );
    LAPACK_dgetrf(&N,&N,A,&LDA,IPIV,&INFO);

    // checks INFO, if INFO != 0 something goes wrong, for more information see the MAN page of dgetrf.
    if(INFO)
    {
        printf("an error occured : \n");
    }else{
        printf("solving the system...\n");
        //void LAPACK_dgetrs( char* trans, lapack_int* n, lapack_int* nrhs, const double* a, lapack_int* lda, const lapack_int* ipiv,double* b, lapack_int* ldb, lapack_int *info );
        dgetrs_(&TRANS,&N,&NRHS,A,&LDA,IPIV,B,&LDB,&INFO);
        if(INFO)
        {
            // checks INFO, if INFO != 0 something goes wrong, for more information see the MAN page of dgetrs.
            printf("an error occured : \n");
        }else{
            printf("print the result : {\n");
            int i;
            for (i=0;i<N;i++)
            {
                printf("%f",B[i]);
            }
            printf("}\n");
        }
    }

    printf("program terminated.");*/
}
    return 0;
}
