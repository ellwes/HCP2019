#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "pi.h"

#define SEED            921
#define NUM_ITER        10000000

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
        
	MPI_File file; 
	MPI_File_open(MPI_COMM_WORLD, "result.txt", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file);
	int size = 4+sizeof(float) + 3;
	char str[size]; 
	int offset = world_rank*size; 
	
	sprintf(str, "%d %f\n", world_rank, (float)count/flip); 
	printf(str);

	//char input[] = "%d\n";  

	MPI_File_write_at(file, offset, str, size, MPI_CHAR, MPI_STATUS_IGNORE);	
	printf("rank %d: %d / %d = %f\n", world_rank, count, flip, (double)count / (double)flip);
        
	int sum = 0; 
	MPI_Reduce(&count, &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        
	if (world_rank == 0) {
		printf("HEEEEEEEJ %d", sum); 
		pi = (float)sum/(float)NUM_ITER*4; 
		offset = size*nr_ranks;
		size = 5+7;
		char res_str[size];
		sprintf(res_str, "pi = %f", (float)pi);
		MPI_File_write_at(file, offset, res_str, size, MPI_CHAR, MPI_STATUS_IGNORE);


	}
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









