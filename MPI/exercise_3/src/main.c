#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "pi.h"

#define SEED            921
#define NUM_ITER        1000000

void print_usage();

int main(int argc, char *argv[])
{
        int opt;
        int world_rank;
        int count = 0, flip = 10000, seed = 1;
        double pi = 0.0;
        char *filename = NULL;
        int nr_ranks;

        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nr_ranks);       
	

        while ((opt = getopt(argc, argv, "s:f:o:")) != -1) {
                switch (opt) {
                        case 's':
                                seed = atoi(optarg);
                                break;
                        case 'f':
                                flip = atoi(optarg);
                                break;
                        case 'o':
                                filename = optarg;
                                break;
                        default:
                                if (world_rank == 0) print_usage(argv[0]);
                }
        }

        init_pi(seed, filename);

        flip = NUM_ITER/nr_ranks;

        //compute_pi(flip, &count, &pi);
        
        // Code starts
        double x, y, z;

        srand(SEED*world_rank);
        //räkna ihop för en look
        int iter;
	for (iter = 0; iter < flip; iter++) {
                x = (double)random() / (double)RAND_MAX;
                y = (double)random() / (double)RAND_MAX;
                z = sqrt((x*x) + (y*y));

         if (z <= 1.0) {
                        count++;
                }
        }
        printf("rank %d: %d / %d = %f\n", world_rank, count, flip, (double)count / (double)flip);
        int recv_count[nr_ranks];
        
	MPI_Gather(&count, 1, MPI_INT, recv_count, 1, MPI_INT, 0, MPI_COMM_WORLD);	
        double sum = 0;
        int i;
	for(i = 0; i < nr_ranks; i++) {
                        sum += (double)recv_count[i];                                                   
        }
        pi =  (sum / (double)NUM_ITER) * 4.0;
        printf("pi: %f\n", pi);
        

        //Code ends

        cleanup_pi();
        MPI_Finalize();

        return 0;
}

void print_usage(char *program)
{
        fprintf(stderr, "Usage: %s [-s seed] [-f trials]\n", program);
        exit(1);
}









