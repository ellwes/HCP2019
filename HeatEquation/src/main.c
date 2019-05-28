#include <stdio.h>
#include "mpi.h"
#include "calc_heat.h"

#define K 0.02

int main(int argc, char *argv[]) {
	/* Variables: */
	//const double K = 0.02;
	int size = 100;
	const double heat_start = 350;
	const double cold_start = 250;
	const double time = 10;
	char * outfile = "outfile.ans"; 
	
	MPI_Init(&argc, &argv);
	
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
	
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	
	/* print given time*/
	if (world_rank == 0) {
		double run_time = end_time - start_time;
		printf("duration: %f", run_time);
	}
	
	MPI_Finalize();

	return 0;
}

