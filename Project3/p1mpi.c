#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#define MIN(a,b)  ((a)<(b)?(a):(b))
//newest

#define BLOCK_LOW(id, p, n) ((id)*(n)/(p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW((id)+1, p,n)-1)
#define BLOCK_SIZE(id, p, n) (BLOCK_LOW((id)+1,p,n)-BLOCK_LOW((id),p,n))
#define BLOCK_OWNER(index,p,n) ((((p) * index)+1)-1)/(n)

int main(int argc, char *argv[])
{
    double elapsed_time;
    int id, index,p,count,nodes;
    unsigned long long int y,z,n,low_value, global_count, high_value, size, proc0_size,i,prime,first;
    char *marked;
    
    nodes = atoi(argv[2]);
    MPI_Init(&argc, &argv);
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    if (argc != 3) {
          if (!id) printf ("Command line: %s <m>\n", argv[0]);
          MPI_Finalize(); exit(1);
    }
    n = atoll(argv[1]);
    //printf("\n n is = %llu \n",n);
    //printf("\n argv[1] = %s \n",argv[1]);

    low_value =  3 + BLOCK_LOW(id,p,n-2) + BLOCK_LOW(id,p,n-2) % 2; //(2 * BLOCK_LOW(id,p,n-1)) + 3;
    high_value = 3 + BLOCK_HIGH(id,p,n-2) - BLOCK_HIGH(id,p,n-2) % 2; //(2 * BLOCK_HIGH(id,p,n-1)) + 3;
    size = (high_value - low_value) / 2 + 1;    //BLOCK_SIZE(id,p,n-1);
    proc0_size = ((n-2)/(2*p)); //(n-1)/p;
    if ((3 + proc0_size) < (int) sqrt((double) n)) {
        if (!id) printf ("Too many processes\n");
        MPI_Finalize();
        exit(1);
    }

    marked = (char *) malloc (size);
    if (marked == NULL) {
        printf("Cannot allocate enough memory\n");
        MPI_Finalize();
        exit(1);
    }
    for (i = 0; i < size; i++){
        marked[i] = 0;
    }
    if (!id) index = 0;
    prime = 3;
    do{
        if (prime * prime > low_value){
            first = (prime * prime - low_value) / 2;
        }
        else{
            if (!(low_value % prime)){
                first = 0;
            }

            else{

                if((low_value % prime) % 2 == 0)
                {

                        first = prime - (low_value % prime) / 2;
                }

            else{
                        first = (prime - (low_value % prime)) / 2;
                }

            }
        }
        for (i = first; i < size; i += prime){
            marked[i] = 1;
        }
        if (!id) {
            while (marked[++index]);
            prime = 2 * index + 3;
        }
        if(p>1)
                MPI_Bcast(&prime,  1, MPI_INT, 0, MPI_COMM_WORLD);
    }while (prime * prime <= n);
    count = 0;
    for (i = 0; i < size; i++){
        if (!marked[i]){
            count++;
        }
    }
    MPI_Reduce (&count, &global_count, 1, MPI_INT, MPI_SUM,0, MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();
    if (!id) {
                global_count++;
        printf("Total number of primes: %llu, Total time: %10.6f sec, Total nodes: %d\n",global_count,elapsed_time,nodes);
        // printf ("Total elapsed time: %10.6f\n", elapsed_time);
    }
    MPI_Finalize();
    return 0;
}
