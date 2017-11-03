#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
// #include "lapacke.h"
// #include "blas.h"

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

double randnum(){
    double z;
    int lower = 0;
    int upper = 10;
    z = ((double)rand()/(RAND_MAX))*(upper-lower)+lower;
    return z;
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


void FillMatrix(double *a, int n){
    int i;
    for(i=0;i<n;i++){
        a[i] = randnum();
    }
}

//void CopyMatrix(double *X, double *Y, int n){
//    int i,j;
//    for(i=0;i<n;i++){
//        Y[i] = X[i];
//    }
//}

//void printArray(double *a, int n, int d){
//    int i,j;
//    if(d==2){
//        for(i=0;i<n;i++){
//            for(j=0;j<n;j++){
//                printf("%f ",a[i*n+j]);
//            }
//            printf("\n");
//        }
//    }
//    else{
//        for(i=0;i<n;i++){
//            printf("%f ",a[i]);
//        }
//        printf("\n");
//    }
//}

double checkCorrectness(double *x1, double *x2, int n){
    int i,j;
    double diff,error = 0.0;
//    for(i=0;i<n;i++){
//        if(error < abs(a[i]-b[i]))
//            error = abs(a[i]-b[i]);
//    }
//    printf("Error = %f\n",error);

for(i=0; i<n; i++)
        {
            diff = abs(x1[i]-x2[i]);
            if(diff>error)
                error = diff;
        }
        printf("Error = %f \n", error);
        if(error < 0.0001)
            printf("Error < 1e-3 \n");

    printf("\n");
}




void mydgetrfBLOCKED(double *A,int *pvt, double *tempv, int n, int b){
    int i,t,l,m,p,q,j,k,maxind,temps,end,ib;
    double max, matsum;
    double *ll;
    for(ib=0;ib<n;ib+=b){
        end = ib + b-1;
        for(i=ib;i<=end;i++){
            maxind = i;
            max=abs(A[i*n+i]);
            for(t=i+1;t<n;t++){
                if(abs(A[t*n+i])>max){
                    maxind = t;
                    max = abs(A[t*n+i]);
                }
            }
            if(max==0.0){
                printf("LU factorization failed: coefficient matrix is singular\n");
                return;
            }
            else{
                if(maxind != i){
                    //Save pivoting information
                    temps = pvt [i];
                    pvt[i] = pvt[maxind];
                    pvt[maxind] = temps;
                    //Swap rows
                    for(k=0;k<n;k++){
                        tempv[k] = A[i*n+k];
                        A[i*n+k] = A[maxind*n+k];
                        A[maxind*n+k] = tempv[k];
                    }
                }
            }

            // printf("PVTed\n");
            // printArray(arrA,n,2);
            // printf("\n");
            for(j=i+1;j<n;j++){
                A[j*n+i] = A[j*n+i]/A[i*n+i];
                for(k=i+1;k<=end;k++){
                    A[j*n+k] = A[j*n+k] - A[j*n+i] * A[i*n+k];
                }
            }
            // printf("FacA\n");
            // printArray(arrA,n,2);
            // printf("\n");
        }

        ll = (double*)calloc(sizeof(double), b*b);
        //double *lln = (double*)calloc(sizeof(double), b*b);
        p=0;q=0;
        for(l=ib;l<=end;l++){
            for(m=ib;m<=end;m++){
                if(l>m){
                    ll[p*b+q] = A[l*n+m] * (-1);
                }
                else if(l==m){
                    ll[p*b+q] = 1;
                }
                else{
                    ll[p*b+q] = 0;
                }
                q++;
            }
            p++;
            q=0;
        }

        // printf("LL\n");
        // printArray(ll,b,2);
        // printf("\n");
        //matInverse(ll,lln,b);
        p=0;q=0;
        for(j=ib;j<=end;j++){
            //arrA[j*n+i] = arrA[j*n+i]/arrA[i*n+i];
            for(k=end+1;k<n;k++){
                matsum = 0.0;
                for(m=ib;m<=end;m++){
                    matsum += ll[p*b+q] * A[m*n+k];
                    q++;
                }
                A[j*n+k] = matsum;
                q=0;
            }
            p++;
            q=0;
        }

        // printf("LL mult\n");
        // printArray(arrA,n,2);
        // printf("\n");
        for(j=end+1;j<n;j++){
            for(k=end+1;k<n;k++){
                double gmat = 0.0;
                for(l=ib;l<=end;l++){
                    gmat += A[j*n+l] * A[l*n+k];
                }
                A[j*n+k] -= gmat;
            }
        }
        //
        // printf("Green\n");
        // printArray(arrA,n,2);
        // printf("\n");
        //free(ll);
    }

            //free(lln);
}


void mydgetrf(double *A,int *pvt, double *tempv, int n){
    int i,t,j,k,maxind,temps;
    double max;
    for(i=0;i<n-1;i++){
        maxind = i;
        max=abs(A[i*n+i]);
        for(t=i+1;t<n;t++){
            if(abs(A[t*n+i])>max){
                maxind = t;
                max = abs(A[t*n+i]);
            }
        }
        if(max==0.0){
            printf("LU factorization failed: coefficient matrix is singular\n");
            return;
        }
        else{
            if(maxind != i){
                //Save pivoting information
                temps = pvt[i];
                pvt[i] = pvt[maxind];
                pvt[maxind] = temps;
                //Swap rows
                for(k=0;k<n;k++){
                    tempv[k] = A[i*n+k];
                    A[i*n+k] = A[maxind*n+k];
                    A[maxind*n+k] = tempv[k];
                }
            }
        }
        for(j=i+1;j<n;j++){
            A[j*n+i] = A[j*n+i]/A[i*n+i];
            for(k=i+1;k<n;k++){
                A[j*n+k] = A[j*n+k] - A[j*n+i] * A[i*n+k];
            }
        }
    }
}

void mydtrsmF(int n, double *A, double *B, int *pvt, double *x, double *y){
    double sum = 0.0, temp;
    int i,k;
    y[0] = B[pvt[0]];
    for(i=1;i<n;i++){
        sum = 0.0;
        for(k=0;k<i;k++){
            sum += y[k]*A[i*n+k];
        }
        y[i] = B[pvt[i]]-sum;
    }
}

void mydtrsmB(int n, double *A, double *B, int *pvt, double *x, double *y){
    double sum = 0.0, temp;
    int i,k;
    x[n-1] = y[n-1]/A[(n-1)*n+(n-1)];
    for(i=n-2;i>=0;i--){
        sum=0.0;
        for(k=i+1;k<n;k++){
            sum+= x[k]*A[i*n+k];
        }
        x[i] = (y[i]-sum)/A[i*n+i];
    }
}

int main(){
    srand((double)time(NULL));
    int *pvt,n,i,j,k,l,m;
    //int ubound = 100, lbound = 0;
    int N[] = {3000};//00,2000,3000,4000,5000};
    int blockSize[] = {10,50,100,200,500};//0,100,200,300,400,500};
    //double random = randomNumber(ubound,lbound);
    double z = randnum();
    double cpu_time_used,gflops;
    int num = sizeof(N)/sizeof(N[0]);
    int bnum = sizeof(blockSize)/sizeof(blockSize[0]);
    for(i=0;i<num;i++)
     {
        n = N[i];
        struct timespec tstart={0,0},tend={0,0};

        double  *A, *A1, *A2, *B, *B1, *B2, *x, *y, *tempv, *x1;
        int *pvt;

        char TRANS = 'N';
        int INFO = n;
        int LDA = n;
        int LDB = n;
        int N = n;
        int NRHS = 1;
        int *IPIV = (int *)calloc(sizeof(int),n);

        char     SIDE = 'L';
        char     UPLO = 'L';
        char     DIAG = 'U';
        int      M    = 1;
        double   a    = 1.0;
        A = (double *)calloc(sizeof(double),n*n);
        A1 = (double *)calloc(sizeof(double),n*n);
        B = (double *)calloc(sizeof(double),n);
        B1 = (double *)calloc(sizeof(double), n);
        tempv = (double *)calloc(sizeof(double), n);
        x = (double *)calloc(sizeof(double), n);
        y = (double *)calloc(sizeof(double), n);
        pvt = (int *)calloc(sizeof(int), n);
        for(m=0;m<n;m++){
            pvt[m]=m;
        }
        FillMatrix(A,n*n);
        //CopyMatrix(A,A1,n*n);

        for(l=0;l<n*n;l++)
            A1[l] = A[l];

        FillMatrix(B,n);
        // printArray(arrB,n,1);
        for(l=0;l<n;l++)
            B1[l] = B[l];

 printf("\n");
printf("Part1.2: My GEPP Functions\n");
printf("Matrix size = %d",n);

        clock_gettime(CLOCK_MONOTONIC, &tstart);
        mydgetrf(A,pvt,tempv,n);
        clock_gettime(CLOCK_MONOTONIC, &tend);

        //cpu_time_used = ((double) (end - start));
        cpu_time_used = ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec);
        gflops = (2*pow(n,3)) /(3*cpu_time_used * pow(10,9));


printf("\n Solved LU Factorization --- Exec. time = %.16f \n", cpu_time_used / CLOCKS_PER_SEC);
printf("\n Solved LU Factorization --- Gflops = %.16f \n", gflops);
        mydtrsmF(n,A,B,pvt,x,y);
        mydtrsmB(n,A,B,pvt,x,y);




                // printf("\n");
//        printf("\n");
//        //printArray(x,n,1);
//        printf("\n");
//        printf("\nBLOCKED GEPP \n");
//        printf("\nSize N = %d\n",n);


        free(y);
        free(tempv);
        free(pvt);

printf("\n");
printf("Part2: Blocked GEPP Functions\n");
printf("Matrix size = %d",n);

        for(k=0;k<bnum;k++){

            printf("\n Block size = %d",blockSize[k]);
            A2 = (double *)calloc(sizeof(double),n*n);
            B2 = (double *)calloc(sizeof(double),n);
            tempv = (double *)calloc(sizeof(double),n);
            x1 = (double *)calloc(sizeof(double), n);
            y = (double *)calloc(sizeof(double), n);
            pvt = (int *)calloc(sizeof(int), n);
            for(m=0;m<n;m++){
                pvt[m]=m;
            }
            //CopyMatrix(A1,A2,n*n);
            //CopyMatrix(B1,B2,n);

            for(l=0;l<n*n;l++)
                A2[l] = A1[l];

            for(l=0;l<n;l++)
                B2[l] = B1[l];


            clock_gettime(CLOCK_MONOTONIC, &tstart);
            mydgetrfBLOCKED(A2,pvt,tempv,n,blockSize[k]);
            clock_gettime(CLOCK_MONOTONIC, &tend);

        cpu_time_used = ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec);
        gflops = (2*pow(n,3)) /(3*cpu_time_used * pow(10,9));


printf("\n Solved LU Factorization --- Exec. time = %.16f \n", cpu_time_used / CLOCKS_PER_SEC);
printf("\n Solved LU Factorization --- Gflops = %.16f \n", gflops);


            mydtrsmF(n,A2,B2,pvt,x1,y);
            mydtrsmB(n,A2,B2,pvt,x1,y);
            printf("\n");




            checkCorrectness(x1,x,n);



        }

        printf("\n");
    }
    return 0;
}
[rgura001@head Project_2]$
