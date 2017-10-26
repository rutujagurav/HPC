#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double* genRand(double *z, int n, double lower, double upper)
{

    int i;

    for(i=0; i<n*n; i++)
    {
        z[i] = ((double)rand()/RAND_MAX)*(upper-lower) + lower; //double precision floating point random numbers between [x,y]
        //printf("%f ", z[i]);

    }
    //printf("\nDone printing...\n");
    return z;
}

double* FillMatrix(double* z, int n)
{
    printf("\n");
    double lower = 0.0;
    double upper = 10.0;
    z = genRand(z, n, lower, upper);

    //printf("...Printing matrix generated...");
    //int i;
    /*for(i=0; i<n*n; i++)
    {
        printf("%f ",z[i]);
    }

    return z;*/
}

void dgemm0(double *a, double *b, double *c, int n)
    {

        int i, j, k;
        for (i=0; i<n; i++)
            for (j=0; j<n; j++)
                for (k=0; k<n; k++)
                {
                    c[i*n+j] += a[i*n+k] * b[k*n+j];
                }

    }

void EightReg(double* a, double* b, double* c2, int n)
{
    int i,j,k;
    for(i = 0; i < n; i += 2)
       for(j = 0; j < n; j += 2)  {
            register int t = i*n+j; register int tt = t+n;
            register double c00 = c2[t]; register double c01 = c2[t+1];  register double c10 = c2[tt]; register double c11 = c2[tt+1];

            for(k = 0; k < n; k += 2) {
                /* 2 by 2 mini matrix multiplication using registers*/
                register int ta = i*n+k; register int tta = ta+n; register int tb = k*n+j; register int ttb = tb+n;
                register double a00 = a[ta]; register double a10 = a[tta]; register double b00 = b[tb]; register double b01 = b[tb+1];

                c00 += a00*b00 ; c01 += a00*b01 ; c10 += a10*b00 ; c11 += a10*b01 ;

                a00 = a[ta+1]; a10 = a[tta+1]; b00 = b[ttb]; b01 = b[ttb+1];

                c00 += a00*b00 ; c01 += a00*b01 ; c10 += a10*b00 ; c11 += a10*b01 ;
             }
             c2[t] = c00;
             c2[t+1] = c01;
             c2[tt] = c10;
             c2[tt+1] = c11;
        }


}

void dgemm2(double* a, double* b, double* c2, int n)
    {
        int i, j, k;
        for(i = 0; i < n; i += 2)
        {
            for(j = 0; j < n; j += 2)
                {
                    register int t = i*n+j; register int tt = t+n;
                    register double c00 = c2[t]; register double c01 = c2[t+1];  register double c10 = c2[tt]; register double c11 = c2[tt+1];
                    //register double c00 = 0.0; register double c01 = 0.0;  register double c10 = 0.0; register double c11 = 0.0;

                    for(k = 0; k < n; k += 2)
                    {
                        //2 by 2 mini matrix multiplication using registers
                        register int ta = i*n+k; register int tta = ta+n; register int tb = k*n+j; register int ttb = tb+n;
                        register double a00 = a[ta]; register double a01 = a[ta+1]; register double a10 = a[tta]; register double a11 = a[tta+1];
                        register double b00 = b[tb]; register double b01 = b[tb+1]; register double b10 = b[ttb]; register double b11 = b[ttb+1];
                        c00 += a00*b00 + a01*b10;
                        c01 += a00*b01 + a01*b11;
                        c10 += a10*b00 + a11*b10;
                        c11 += a10*b01 + a11*b11;
                    }

                     c2[t] = c00;
                     c2[t+1] = c01;
                     c2[tt] = c10;
                     c2[tt+1] = c11;
                }
        }
        //return c2;
    }

void dgemm3(double* a, double* b, double* c2, int n)
    {

        {
        int i, j, k;
        for(i = 0; i < n; i += 2)
        {
            for(j = 0; j < n; j += 2)
                {
                    register int t = i*n+j; register int tt = t+n;
                    register double c00 = c2[t]; register double c01 = c2[t+1];  register double c10 = c2[tt]; register double c11 = c2[tt+1];
                    //register double c00 = 0.0; register double c01 = 0.0;  register double c10 = 0.0; register double c11 = 0.0;

                    for(k = 0; k < n; k += 2)
                    {
                        //2 by 2 mini matrix multiplication using registers
                        register int ta = i*n+k; register int tta = ta+n; register int tb = k*n+j; register int ttb = tb+n;
                        register double a00 = a[ta]; register double a01 = a[ta+1]; register double a10 = a[tta]; register double a11 = a[tta+1];
                        register double b00 = b[tb]; register double b01 = b[tb+1]; register double b10 = b[ttb]; register double b11 = b[ttb+1];

                        register double temp1;
                        register double temp2;
                        register double temp3;
                        register double temp4;

                        /*c00 += a00*b00 + a01*b10;
                        c01 += a00*b01 + a01*b11;
                        c10 += a10*b00 + a11*b10;
                        c11 += a10*b01 + a11*b11;*/

                        temp1 = a00*b00;
                        temp2 = a00*b01;
                        temp3 = a10*b00;
                        temp4 = a10*b01;

                        temp1 += a01*b10;
                        temp2 += a01*b11;
                        temp3 += a11*b10;
                        temp4 += a11*b11;

                        c00 += temp1;
                        c01 += temp2;
                        c10 += temp3;
                        c11 += temp4;


                        /*temp1 = a00*b00;
                        temp2 = a01*b10;


                        c00 += temp1 + temp2;

                        temp1 = a00*b01;
                        temp2 = a01*b11;

                        c01 += temp1 + temp2;

                        temp1 = a10*b00;
                        temp2 = a11*b10;

                        c10 += temp1 + temp2;

                        temp1 = a10*b01;
                        temp2 = a11*b11;

                        c11 += temp1 + temp2;*/
                    }

                     c2[t] = c00;
                     c2[t+1] = c01;
                     c2[tt] = c10;
                     c2[tt+1] = c11;
                }
        }
        //return c2;


    }

    }


int main()
{
    srand((unsigned)time(NULL)); //setting seed to system time so that it gives different random numbers each time.
    int n[] = {64,128,256};
    //int n = 64;
    int i,j,k;
    double *z = (double *) calloc(sizeof(double), 1);

    z = genRand(z, 1, 1.0, 10.0); //wasted float value

    printf("********** PART #3 16 Register Reuse ********* \n");


    for(k=0; k<3; k++)
    {
        printf("MATRIX SIZE = %d \n", n[k]);
        double *a, *b, *c1, *c2, *c3;

        a = (double *) calloc(sizeof(double), n[k]*n[k]);
        b = (double *) calloc(sizeof(double), n[k]*n[k]);
        c1 = (double *) calloc(sizeof(double), n[k]*n[k]);
        c2 = (double *) calloc(sizeof(double), n[k]*n[k]);
        c3 = (double *) calloc(sizeof(double), n[k]*n[k]);

        double lower = 0.0;
        double upper = 10.0;

        clock_t start, end;
        double cpu_time_used, gflops, error = 0.0;
        double diff;

        //Generating random numbers and filling arrays a and b.


        a = FillMatrix(a,n[k]);
        b = FillMatrix(b,n[k]);
        c1 = FillMatrix(c1,n[k]);
        //c2 = FillMatrix(c2,n[k]);


        for(i=0; i<n[k]*n[k]; i++)
        {
            c2[i] = c1[i];
            c3[i] = c1[i];
        }


        //Matrix Multiplication algorithms

        //start = clock();
        struct timespec tstart={0,0}, tend={0,0};
		clock_gettime(CLOCK_MONOTONIC, &tstart);
        dgemm0(a, b, c1, n[k]);
        clock_gettime(CLOCK_MONOTONIC, &tend);
		cpu_time_used = ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec);
        //end = clock();
        //cpu_time_used = ((double) (end - start));
        gflops = ((double)2*pow(n[k],3) /(double)(cpu_time_used * pow(10,9)));
        //printf("MATRIX SIZE = %d \n", n[k]);
        printf("\n Returned from dgemm0() --- Exec. time = %.16f \n", cpu_time_used / CLOCKS_PER_SEC);
        printf("\n Returned from dgemm0() --- Gflops = %.16f \n", gflops);


        //start = clock();
        clock_gettime(CLOCK_MONOTONIC, &tstart);
        dgemm2(a, b, c2, n[k]);
        clock_gettime(CLOCK_MONOTONIC, &tend);
		cpu_time_used = ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec);
        //end = clock();
        //cpu_time_used = ((double) (end - start));
        gflops = ((double)2*pow(n[k],3) /(double)(cpu_time_used * pow(10,9)));
        //printf("MATRIX SIZE = %d \n", n[k]);
        printf("\n Returned from dgemm2() --- Exec. time = %.16f \n", cpu_time_used / CLOCKS_PER_SEC);
        printf("\n Returned from dgemm2() --- Gflops = %.16f \n", gflops);

        //printf("dgemm0 PRINTING ANSWER MATRIX \n");

//        for(i=0; i<n[k]*n[k]; i++)
//            printf("%f \n",c1[i]);


        //start = clock();
        clock_gettime(CLOCK_MONOTONIC, &tstart);
        dgemm3(a, b, c3, n[k]);
        clock_gettime(CLOCK_MONOTONIC, &tend);
		cpu_time_used = ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec);
        //end = clock();
        //cpu_time_used = ((double) (end - start));
        gflops = ((double)2*pow(n[k],3) /(double)(cpu_time_used * pow(10,9)));
        //printf("MATRIX SIZE = %d \n", n[k]);
        printf("\n Returned from dgemm3() --- Exec. time = %.16f \n", cpu_time_used / CLOCKS_PER_SEC);
        printf("\n Returned from dgemm3() --- Gflops = %.16f \n", gflops);

        //for(i=0; i<n[k]*n[k]; i++)
          //  printf("%f %f \n",c1[i],c2[i]);



        for(i=0; i<n[k]*n[k]; i++)
        {
            diff = abs(c1[i]-c3[i]);
            if(diff>error)
                error = diff;
        }

        printf("Error for dgemm0 and dgemm3 = %f \n", error);
        if(error < 0.0001)
            printf("Error < 1e-3 \n");

        for(i=0; i<n[k]*n[k]; i++)
        {
            diff = abs(c2[i]-c3[i]);
            if(diff>error)
                error = diff;
        }

        printf("Error for dgemm2 and dgemm3 = %f \n", error);
        if(error < 0.0001)
            printf("Error < 1e-3 \n");


        printf("\n");

    }
    return 0;
}
