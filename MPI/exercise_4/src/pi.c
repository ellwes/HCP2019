#include "pi.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define SEED     921


void init_pi(int set_seed, char *outfile)
{
	if (filename != NULL) {
		free(filename);
		filename = NULL;
	}

	if (outfile != NULL) {
		filename = (char*)calloc(sizeof(char), strlen(outfile)+1);
		memcpy(filename, outfile, strlen(outfile));
		filename[strlen(outfile)] = 0;
	}
	seed = set_seed;
}

void cleanup_pi()
{
	if (filename != NULL)
		free(filename);
}

void compute_pi(int flip, int *local_count, double *answer)
{
	//printf("hello");
 	
        double x, y, z;
        ////int local_nr = NUM_ITER/MPI_COMM_WORLD;
	int world_rank;
	int nr_ranks;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nr_ranks);       

	srand(SEED*world_rank);
        //räkna ihop för en look
        int iter;
	for (iter = 0; iter < flip; iter++) {
                x = (double)random() / (double)RAND_MAX;
                y = (double)random() / (double)RAND_MAX;
                z = sqrt((x*x) + (y*y));

         if (z <= 1.0) {
                        (*local_count)++;
                }
        }
	printf("SENT COUNT: %d\n", *local_count);

	int sum = 0;
       	MPI_Reduce(local_count, &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
       	*answer =  (double)sum / ((double)(flip*nr_ranks)) * 4.0;

}
