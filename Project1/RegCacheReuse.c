#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

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

void ver01(double *a, double *b, double *c, int n)
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
}


void RegCacheMM(double* a, double* b, double* c2, int n, int B)
{

    int i, j, k;
    int i1,j1,k1;
    //int B = 16;
    for (i=0;i<n;i+=B)
    for (j=0;j<n;j+=B)
    for (k=0;k<n;k+=B)
    for (i1=i; i1<i+B; i1+=2)
        for (j1=j; j1<j+B; j1+=2)
        {
            register int t = i1*n+j1; register int tt = t+n;
            register double c00 = c2[t]; register double c01 = c2[t+1];  register double c10 = c2[tt]; register double c11 = c2[tt+1];
            //register double c00 = 0.0; register double c01 = 0.0;  register double c10 = 0.0; register double c11 = 0.0;
            for (k1=k; k1<k+B; k1+=2)
            {
                //2 by 2 mini matrix multiplication using registers
                        register int ta = i1*n+k1; register int tta = ta+n; register int tb = k1*n+j1; register int ttb = tb+n;
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
    //return c;
}




int main()
{
    srand((unsigned)time(NULL)); //setting seed to system time so that it gives different random numbers each time.
    //int n[] = {64,128,256,512,1024,2048};
    int n = 2048;
    int i,j,k,m;
    double *z = (double *) calloc(sizeof(double), 1);
    z = genRand(z, 1, 1.0, 10.0); //wasted float value

    //printf("%f ",*z);

    //for(k=0; k<6; k++)
    //{

    printf("********** PART #4: CACHE BLOCKING & REGISTER BLOCKING ********** \n");
        double *a, *b, *c, *c1;

        a = (double *) calloc(sizeof(double), n*n);
        b = (double *) calloc(sizeof(double), n*n);
        c = (double *) calloc(sizeof(double), n*n);
        c1 = (double *) calloc(sizeof(double), n*n);


        double lower = 0.0;
        double upper = 10.0;

        clock_t start, end;
        double cpu_time_used, gflops, error = 0.0;
        double diff;
        int B[] = {4,8,16,32,64,128,256,512};
        //Generating random numbers and filling arrays a and b.


        a = FillMatrix(a,n);
        b = FillMatrix(b,n);
        c = FillMatrix(c,n);

        for(i=0; i<n*n; i++)
            {
                c1[i] = c[i];
                //printf("A: %f & B: %f & C1: %f & C2: %f\n",a[i], b[i], c[i], c1[i]);
            }


         for(m=0; m<3; m++)
        {
            printf("MATRIX SIZE = %d \n", n);
            printf("BLOCK SIZE = %d \n", B[m]);
        //Matrix Multiplication algorithm calls

            ver01(a, b, c1, n); //just computing for correctness
    //        for(i=0; i<n*n; i++)
    //            printf("c1: %f \n", c1[i]);

                printf("\n");
            start = clock();
            RegCacheMM(a, b, c, n, B[m]);
            end = clock();
    //        for(i=0; i<n*n; i++)
    //            printf("c: %f \n", c[i]);

            cpu_time_used = ((double) (end - start));
            gflops = ((double)2*pow(n,3) /(double)(cpu_time_used * pow(10,9)));
            printf("Returned from Register and Cache blocking Algorithm RegCacheMM() \n");
            printf("--- Exec. time = %.16f \n", cpu_time_used / CLOCKS_PER_SEC);
            printf("--- Gflops = %.16f \n", gflops);

            for(i=0; i<n*n; i++)
            {

                diff = abs(c1[i]-c[i]);
                //printf("%f\n",diff);
                if(diff>error)
                    error = diff;
            }

            //for(i=0; i<n; i++)
              //      error = maximum(abs(c1[i]-c[i]), error);

            printf("Error for ver01() and RegCacheMM() = %f \n", error);
            if(error < 0.0001)
                printf("Error < 1e-3 \n");
        }


    return 0;
}
