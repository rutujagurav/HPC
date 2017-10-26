#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double* genRand(double *z, int n, double lower, double upper)
{

    int i;

    for(i=0; i<n*n; i++)
        z[i] = ((double)rand()/RAND_MAX)*(upper-lower) + lower; //double precision floating point random numbers between [x,y]

    return z;
}


double* FillMatrix(double* z, int n)
{
    printf("\n");
    double lower = 0.0;
    double upper = 10.0;
    z = genRand(z, n, lower, upper);
    return z;
}


//------------------------------------------------------------------------------------------------------------------------------------------

double* ver01(double *a, double *b, double *c, int n)
{
    /* ijk – simple triple loop algorithm with simple single register reuse*/
    int i, j, k;
    for (i=0; i<n; i++)
        for (j=0; j<n; j++)
        {
            register double r=c[i*n+j];
            //register double r=0.0;
            for (k=0; k<n; k++)
                r += a[i*n+k] * b[k*n+j];
            c[i*n+j]=r;
        }
    return c;
}

double* ver02(double *a, double *b, double *c, int n)
{
    /* jik – simple triple loop algorithm with simple single register reuse*/
    int i, j, k;
    for (j=0; j<n; j++)
        for (i=0; i<n; i++)
        {
            register double r=c[i*n+j];
            //register double r=0.0;
            for (k=0; k<n; k++)
                r += a[i*n+k] * b[k*n+j];
            c[i*n+j]=r;
        }
    return c;
}


double* ver03(double *a, double *b, double *c, int n)
{
    /* kij – simple triple loop algorithm with simple single register reuse*/
    int i, j, k;
    for (k=0; k<n; k++)
        for (i=0; i<n; i++)
        {
            //register double r=c[i*n+j];
            register double r=a[i*n+k];
            for (j=0; j<n; j++)
                c[i*n+j] += r * b[k*n+j];
        }
    return c;
}

double* ver04(double *a, double *b, double *c, int n)
{
    /* ikj – simple triple loop algorithm with simple single register reuse*/
    int i, j, k;
    for (i=0; i<n; i++)
        for (k=0; k<n; k++)
        {
            //register double r=c[i*n+j];
            register double r=a[i*n+k];
            for (j=0; j<n; j++)
                c[i*n+j] += r * b[k*n+j];
        }
    return c;
}

double* ver05(double *a, double *b, double *c, int n)
{
    /* jki – simple triple loop algorithm with simple single register reuse*/
    int i, j, k;
    for (j=0; j<n; j++)
        for (k=0; k<n; k++)
        {
            //register double r=c[i*n+j];
            register double r=b[k*n+j];
            for (i=0; i<n; i++)
                c[i*n+j] += a[i*n+k] * r;
        }
    return c;
}

double* ver06(double *a, double *b, double *c, int n)
{
    /* kji – simple triple loop algorithm with simple single register reuse*/
    int i, j, k;
    for (k=0; k<n; k++)
        for (j=0; j<n; j++)
        {
            //register double r=c[i*n+j];
            register double r=b[k*n+j];
            for (i=0; i<n; i++)
                c[i*n+j] += a[i*n+k] * r;
        }
    return c;
}

//----------------------------------------------------------------------------------------------------------------------------------------

int main()
{
    srand((unsigned)time(NULL)); //setting seed to system time so that it gives different random numbers each time.
    //int n[] = {64,128,256,512,1024,2048};
    int n = 2048;
    int i,j,k;
    double *z = (double *) calloc(sizeof(double), 1);

    z = genRand(z, 1, 1.0, 10.0); //wasted float value


        double *a, *b, *c1, *c2;

        a = (double *) calloc(sizeof(double), n*n);
        b = (double *) calloc(sizeof(double), n*n);
        c1 = (double *) calloc(sizeof(double), n*n);
        c2 = (double *) calloc(sizeof(double), n*n);


        double lower = 0.0;
        double upper = 10.0;

        clock_t start, end;
        double cpu_time_used, gflops, error = 0.0;
        double diff;

        //Generating random numbers and filling arrays a and b.

        a = FillMatrix(a,n);
        b = FillMatrix(b,n);
        c1 = FillMatrix(c1,n);

        /*for(i=0; i<n*n; i++)
            printf("%f \n",a[i]);*/

        printf("**********CACHE REUSE WITHOUT BLOCKING*********** \n");

        start = clock();
        c1 = ver01(a, b, c1, n);
        end = clock();
        cpu_time_used = ((double) (end - start));
        gflops = ((double)2*pow(n,3) /(double)(cpu_time_used * pow(10,9)));
        printf("MATRIX SIZE = %d \n", n);

        printf("\n Returned from ver01 --- Exec. time = %f \n", cpu_time_used / CLOCKS_PER_SEC);
        printf("\n Returned from ver01 --- Gflops = %.16f \n", gflops);

        start = clock();
        c2 = ver02(a, b, c1, n);
        end = clock();
        cpu_time_used = ((double) (end - start));
        gflops = ((double)2*pow(n,3) /(double)(cpu_time_used * pow(10,9)));
        printf("MATRIX SIZE = %d \n", n);
        printf("\n Returned from ver02 --- Exec. time = %f \n", cpu_time_used / CLOCKS_PER_SEC );
        printf("\n Returned from ver02 --- Gflops = %.16f \n", gflops);

        for(i=0; i<n*n; i++)
        {
            diff = abs(c1[i]-c2[i]);
            //printf("Diff = %f",diff);
            if(diff>error)
                error = diff;
        }

        //for(i=0; i<n*n; i++)
           // error = maximum(abs(c1[i]-c2[i]), error);


        printf("Error for ver01 and ver02 = %f \n", error);
        if(error < 0.0001)
            printf("Error < 1e-3 \n");

//----------------------------------------------------------------------------------------
        start = clock();
        c1 = ver03(a, b, c1, n);
        end = clock();
        cpu_time_used = ((double) (end - start));
        gflops = ((double)2*pow(n,3) /(double)(cpu_time_used * pow(10,9)));
        printf("MATRIX SIZE = %d \n", n);
        printf("\n Returned from ver03 --- Exec. time = %f \n", cpu_time_used / CLOCKS_PER_SEC);
        printf("\n Returned from ver03 --- Gflops = %.16f \n", gflops);


        start = clock();
        c2 = ver04(a, b, c1, n);
        end = clock();
        cpu_time_used = ((double) (end - start));
        gflops = ((double)2*pow(n,3) /(double)(cpu_time_used * pow(10,9)));
        printf("MATRIX SIZE = %d \n", n);
        printf("\n Returned from ver04 --- Exec. time = %f \n", cpu_time_used / CLOCKS_PER_SEC);
        printf("\n Returned from ver04 --- Gflops = %.16f \n", gflops);

        for(i=0; i<n*n; i++)
        {
            diff = abs(c1[i]-c2[i]);
            if(diff>error)
                error = diff;
        }

        //for(i=0; i<n*n; i++)
          //  error = maximum(abs(c1[i]-c2[i]), error);


        printf("Error for ver03 and ver04 = %f \n", error);
        if(error < 0.0001)
            printf("Error < 1e-3 \n");

//--------------------------------------------------------------------------------------------

        start = clock();
        c1 = ver05(a, b, c1, n);
        end = clock();
        cpu_time_used = ((double) (end - start));
        gflops = ((double)2*pow(n,3) /(double)(cpu_time_used * pow(10,9)));
        printf("MATRIX SIZE = %d \n", n);
        printf("\n Returned from ver05 --- Exec. time = %f \n", cpu_time_used / CLOCKS_PER_SEC);
        printf("\n Returned from ver05 --- Gflops = %.16f \n", gflops);


        start = clock();
        c2 = ver06(a, b, c1, n);
        end = clock();
        cpu_time_used = ((double) (end - start));
        gflops = ((double)2*pow(n,3) /(double)(cpu_time_used * pow(10,9)));
        printf("MATRIX SIZE = %d \n", n);
        printf("\n Returned from ver06 --- Exec. time = %f \n", cpu_time_used / CLOCKS_PER_SEC);
        printf("\n Returned from ver06 --- Gflops = %.16f \n", gflops);

        for(i=0; i<n*n; i++)
        {
            diff = abs(c1[i]-c2[i]);
            if(diff>error)
                error = diff;
        }
        //for(i=0; i<n*n; i++)
            //error = maximum(abs(c1[i]-c2[i]), error);


        printf("Error = %f \n", error);
        if(error < 0.0001)
            printf("Error < 1e-3 \n");



    return 0;
}
