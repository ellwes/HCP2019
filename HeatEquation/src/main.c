#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "calc_heat.h"

//#define K 0.02

int main(int argc, char *argv[]) {
	/* Variables: */
	const double K = 0.02;
	int size = 1296;
	const double heat_start = 350;
	const double cold_start = 250;
	const double time = 0.05;
	char * outfile = "test.ans"; 
	int repeat = 10;

	
	MPI_Init(&argc, &argv);	
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	double avg_time = 0.0, prev_avg_time = 0.0, stddev_time = 0.0;
 	
	for(int i = 0; i < repeat; i++) {
		/* init solver */
		init_solver(K, size, heat_start, cold_start, time, outfile);	
		
		/* wait for all process to be done with init */
		MPI_Barrier(MPI_COMM_WORLD);
		/* start timer */
		double start_time = MPI_Wtime();
		/* calc heat equation */
		calc_heat();
		/* end timer */
		MPI_Barrier(MPI_COMM_WORLD);
		double end_time = MPI_Wtime();

		/* print given time*/
		if (world_rank == 0) {
			double run_time = end_time - start_time;
			printf("run %d, time: %fs\n", i, run_time);
		}
		prev_avg_time = avg_time; 
		avg_time = avg_time + ( (end_time - start_time) - avg_time) / (i + 1);
		stddev_time = stddev_time + ( (end_time - start_time) - avg_time) * ( (end_time - start_time) - prev_avg_time);	
	}	
	
	if (world_rank == 0) {
		stddev_time = sqrt(stddev_time / (repeat - 1));	
		printf("Duration:  %f±%f\n", avg_time, stddev_time);
	}

	/*cleanup*/
	cleanup();
	
	MPI_Finalize();

	return 0;
}

