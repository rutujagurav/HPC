#include<stdio.h>
#include"lapacke.h"

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
    int n[] = {5,1000,2000,3000,4000,5000};
    int i,j,k;

    // note, to understand this part take a look in the MAN pages, at section of parameters.
    char    TRANS = 'N';
    int     INFO=n[0];
    int     LDA = n[0];
    int     LDB = n[0];
    int     N = n[0];
    int     NRHS = 1;
    int     IPIV[1000] ;

    double *A, *B;

    double *z = (double *) calloc(sizeof(double), 1);
    z = genRand(z, 1, 1, 1.0, 10.0); //wasted float value


    A = (double *) calloc(sizeof(double), n[0]*n[0]);
    B = (double *) calloc(sizeof(double), n[0]*1);

    FillMatrix(A,n[0],n[0]);
    FillMatrix(B,n[0],1);

//    for(i=0; i<n[0]*n[0]; i++)
//        printf("%f \n ",A[i]);
//
//    printf("\n");
//    for(i=0; i<n[0]; i++)
//        printf("%f \n",B[i]);


    /*double  A[9] =
    {
    1, 2, 3,
    2, 3, 4,
    3, 4, 1
    };

    double B[3] =
    {
    -4,
    -1,
    -2
    };*/
// end of declarations

   printf("compute the LU factorization...\n");
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

    printf("program terminated.");
    return 0;
}
