#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#define MIN(a,b)  ((a)<(b)?(a):(b))
//new

#define BLOCK_LOW(id, p, n) ((id)*(n)/(p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW((id)+1, p,n)-1)
#define BLOCK_SIZE(id, p, n) (BLOCK_LOW((id)+1,p,n)-BLOCK_LOW((id),p,n))
#define BLOCK_OWNER(index,p,n) ((((p) * index)+1)-1)/(n)

int main(int argc, char *argv[])
{
    
    int id, index,p,count, nodes;
    unsigned long long int n,i,j,low_value, high_value, size, proc0_size,prime,first;
    char *marked, *sievemarked;
    unsigned long long int global_count;
	unsigned long long int sievelow,sievehigh,sievesize,sievefirst;
    double elapsed_time;

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
    low_value = 3 + BLOCK_LOW(id,p,n-2) + BLOCK_LOW(id,p,n-2) % 2;
    high_value = 3 + BLOCK_HIGH(id,p,n-2) - BLOCK_HIGH(id,p,n-2) % 2;
    size = (high_value - low_value) / 2 + 1;
    proc0_size = ((n-2)/(2*p));

    sievelow = 3;
   // sievehigh = 3 + BLOCK_HIGH(0,p,n-2) - BLOCK_HIGH(0,p,n-2) % 2;
    sievehigh = (unsigned long long int)sqrt(n);
    sievesize = (sievehigh - sievelow) / 2 + 1;
    sievemarked = (char*)malloc(sievesize);

    if (sievemarked == NULL) {
        printf("Cannot allocate memory to local array for seiving primes\n");
        MPI_Finalize();
        exit(1);
    }

    for(j=0;j<sievesize;j++){
        sievemarked[j] = 0;
    }

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

    index = 0;
    prime = 3;
    do{
        if (prime * prime > low_value){
            first = (prime * prime - low_value)/2;
        }
        else{
            if (!(low_value % prime)){
                first = 0;
            }
            else{
                if ((low_value%prime) % 2 == 0)
                                {
                                        first = prime - (low_value%prime) / 2;
                                }
                                else
                                {
                                        first = (prime - (low_value%prime)) / 2;
                                }
            }
        }
        for (i = first; i < size; i += prime){
            marked[i] = 1;
        }
        if (!id) {
            while (marked[++index]);
            prime = 2*index + 3;
        }
        //if(p>1)
            //MPI_Bcast(&prime,  1, MPI_INT, 0, MPI_COMM_WORLD);
        if(id){
            sievefirst = (prime * prime - sievelow)/2;
            for(j=sievefirst;j<sievesize;j+=prime){
                sievemarked[j] = 1;
            }
            while(sievemarked[++index]);
            prime = 2 * index + 3;
        }
    }while (prime * prime <= n);
    count = 0;
    for (i = 0; i < size; i++){
        if (!marked[i]){
            count++;
        }
    }
    if(p>1){
        MPI_Reduce (&count, &global_count, 1, MPI_INT, MPI_SUM,0, MPI_COMM_WORLD);
    }
    elapsed_time += MPI_Wtime();
    if (!id) {
        global_count++;
        printf("Total number of primes: %llu, Total time: %10.6f sec, Total nodes: %d\n",global_count,elapsed_time,nodes);
    }
    MPI_Finalize();
    return 0;
}
