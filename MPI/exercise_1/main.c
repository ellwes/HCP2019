#include <mpi.h>
#include <stdio.h>

int main(int argc, char **argv)
{
	int rank, size, prev_rank, next_rank, recv_rank;
	MPI_Init(&argc, &argv);
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	next_rank = (rank+1) % size;
	prev_rank = (rank-1 + size) % size;
	
	if (rank == 0) {
		MPI_Send(&rank, 1, MPI_INT, next_rank, 0, MPI_COMM_WORLD);
		MPI_Recv(&recv_rank, 1, MPI_INT, prev_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);	
	} else {
		MPI_Recv(&recv_rank, 1, MPI_INT, prev_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);	
		MPI_Send(&rank, 1, MPI_INT, next_rank, 0, MPI_COMM_WORLD);
	}
	printf("%d received %d from rank=%d\n", rank, recv_rank, prev_rank);
	MPI_Finalize();
}
