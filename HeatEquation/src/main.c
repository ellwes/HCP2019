include "calc_heat.h"


int main() {
	/* Variables: */
	const double K = 0.02;
	int size = 100;
	const double heat_start = 350;
	const double cold_start = 250;
	const double time = 10;
	char outfile = "outfile.ans" 

	/* init solver */
	init_solver(K, size, heat_start, cold_start, time, outfile);	
		
	/* wait for all process to be done with init */
	MPI_Barrier(MPI_COMM_WORLD);

	/* start timer */
	start_time = MPI_Wtime();

	/* calc heat equation */
	calc_heat();
	
	/* end timer */
	MPI_Barrier(MPI_COMM_WORLD);
	end_time = MPI_Wtime();


	/* print given time*/
	if (world_rank == 0) {
		run_time = end_time - start_time;
		printf("duration: %f", run_time);
	}

	return 0;
}
