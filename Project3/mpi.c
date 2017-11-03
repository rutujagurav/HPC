#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include "MyMPI.h“ slide 19
#define MIN(a,b)  ((a)<(b)?(a):(b))

int main (int argc, char *argv[])
{
   ...
   MPI_Init (&argc, &argv);
   MPI_Barrier(MPI_COMM_WORLD);
   elapsed_time = -MPI_Wtime();
   MPI_Comm_rank (MPI_COMM_WORLD, &id);
   MPI_Comm_size (MPI_COMM_WORLD, &p);
    if (argc != 2) {
          if (!id) printf ("Command line: %s <m>\n", argv[0]);
          MPI_Finalize(); exit (1);
    }
    n = atoi(argv[1]);
   low_value = 2 + BLOCK_LOW(id,p,n-1);
   high_value = 2 + BLOCK_HIGH(id,p,n-1);
   size = BLOCK_SIZE(id,p,n-1);
   proc0_size = (n-1)/p;
   if ((2 + proc0_size) < (int) sqrt((double) n)) {
      if (!id) printf ("Too many processes\n");
      MPI_Finalize();
      exit (1);
   }

   marked = (char *) malloc (size); //local array of each processor
   if (marked == NULL) {
      printf ("Cannot allocate enough memory\n");
      MPI_Finalize();
      exit (1);
   }
    for (i = 0; i < size; i++) marked[i] = 0; //init mark array to 0
       if (!id) index = 0;//”if(!id)” means you are process 0.
       prime = 2; //initially we try 2 as the testing number (also called sieving prime)
       do { //try to find first number who is multiplier of 2 (or any of the sieving prime). ‘first’ is the index of that first number
          if (prime * prime > low_value) //in case low_value is not the multiple of the current sieving prime, you’ll have to do linear search from low_value till you find it. But instead do something smart like this.)
             first = prime * prime - low_value;
          else {
             if (!(low_value % prime)) first = 0;
             else first = prime - (low_value % prime);
          }
          for (i = first; i < size; i += prime) marked[i] = 1; //marking the multipliers of the sieving prime
          if (!id) { //finding the next sieving prime. In Project 3 Part 2, we make all processors do this at the same time as poc0 and not keep waiting for broadcast. That way you save on communication time as well.
             while (marked[++index]);
             prime = index + 2;
          }
          MPI_Bcast (&prime,  1, MPI_INT, 0, MPI_COMM_WORLD);
       } while (prime * prime <= n);
    count = 0; //now need to count the total no. of primes
       for (i = 0; i < size; i++)
          if (!marked[i]) count++;
       MPI_Reduce (&count, &global_count, 1, MPI_INT, MPI_SUM,
          0, MPI_COMM_WORLD);
       elapsed_time += MPI_Wtime();
       if (!id) {
          printf ("%d primes are less than or equal to %d\n",
             global_count, n);
          printf ("Total elapsed time: %10.6f\n", elapsed_time);
       }
       MPI_Finalize ();
       return 0;
}

