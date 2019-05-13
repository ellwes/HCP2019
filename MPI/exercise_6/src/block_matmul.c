#include "block_matmul.h"

struct Config {
	/* MPI Files */
	MPI_File A_file, B_file, C_file;
	char *outfile;

	/* MPI Datatypes for matrix blocks */
	MPI_Datatype block;

	/* Matrix data */
	double *A, *A_tmp, *B, *C;

	/* Cart communicators */
	MPI_Comm grid_comm;
	MPI_Comm row_comm;
	MPI_Comm col_comm;

	/* Cart communicator dim and ranks */
	int dim[2], coords[2];
	int world_rank, world_size, grid_rank;
	int row_rank, row_size, col_rank, col_size;

	/* Full matrix dim */
	int A_dims[2];
	int B_dims[2];
	int C_dims[2];
	int matrix_size;

        
	int sum; 
	/* Process local matrix dim */
	int local_dims[2];
	int local_size;
};

struct Config config;

void init_matmul(char *A_file, char *B_file, char *outfile)
{
	/* Copy output file name to configuration */
	config.outfile = outfile; 	
	/* Get matrix size header */
	MPI_File_open(MPI_COMM_WORLD, A_file, (MPI_MODE_RDWR | MPI_MODE_CREATE), MPI_INFO_NULL, &config.A_file);		  
	MPI_File_read(config.A_file, config.A_dims, 2, MPI_INT, MPI_STATUS_IGNORE);	
	
	MPI_File_open(MPI_COMM_WORLD, B_file, (MPI_MODE_RDWR | MPI_MODE_CREATE), MPI_INFO_NULL, &config.B_file);
	MPI_File_read(config.B_file, config.B_dims, 2, MPI_INT, MPI_STATUS_IGNORE);	
	config.matrix_size = config.A_dims[0] * config.A_dims[0];

	if(config.world_rank == 0) {
		double test[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
		MPI_File_read(config.A_file, test, 16, MPI_DOUBLE, MPI_STATUS_IGNORE);
		int j; 
		for(j = 0; j < 16; j++) {
			printf("rank; %d, j: %d, value: %f..... \n", config.world_rank, j, test[j]);
		}
	}
	
	
	/* Broadcast global matrix sizes */
	MPI_Bcast(&config.matrix_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

	/* Set dim of tiles relative to the number of processes as NxN where N=sqrt(world_size) */
	MPI_Comm_size(MPI_COMM_WORLD, &(config.world_size));
	MPI_Dims_create(config.world_size, 2, config.dim);
	
	/* Verify dim of A and B matches for matul and both are square*/
	if(!(config.A_dims[0] == config.B_dims[0] && config.A_dims[0] == config.A_dims[1] && config.B_dims[0] == config.B_dims[1])) {
		printf("WRONG DIMENSIONS\n");
	
	}
	/* Create Cart communicator for NxN processes */
	int period[2] = {1, 1};
	MPI_Cart_create(MPI_COMM_WORLD, sqrt(config.world_size), config.dim, period, 1, &(config.grid_comm));

	/* Sub div cart communicator to N row communicator */
	MPI_Comm_rank(MPI_COMM_WORLD, &(config.world_rank));
	
	int remain_dims[2] = {0, 1};
	MPI_Cart_sub(config.grid_comm, remain_dims, &(config.row_comm));
	//MPI_Comm_split(MPI_COMM_WORLD, config.world_rank / config.A_dims[0] , config.world_rank, &(config.row_comm));	
	MPI_Comm_rank(config.row_comm, &(config.row_rank));
	MPI_Comm_size(config.row_comm, &(config.row_size));
	
	/* Sub div cart communicator to N col communicator */
	int remain_dims2[2] = {1, 0};
	MPI_Cart_sub(config.grid_comm, remain_dims2, &(config.col_comm));	
	//MPI_Comm_split(MPI_COMM_WORLD, config.world_rank % config.world_size, config.world_rank, &(config.col_comm));	
	MPI_Comm_rank(config.col_comm, &(config.col_rank));
	MPI_Comm_size(config.col_comm, &(config.col_size));
	
	/* Setup sizes of full matrices */
	//GJORT??!!
	
	/* Setup sizes of local matrix tiles */
	config.local_dims[0] = config.A_dims[0] / config.dim[0];
	config.local_dims[1] = config.A_dims[1] / config.dim[1];
	config.local_size = config.local_dims[0] * config.local_dims[1];

	/* Create subarray datatype for local matrix tile */
	int start[2] = {config.col_rank, config.row_rank};
	MPI_Type_create_subarray(2, config.A_dims, config.local_dims, start, MPI_ORDER_C, MPI_DOUBLE, &(config.block));
	MPI_Type_commit(&(config.block));

	/* Create data array to load actual block matrix data */
	double * A_matrix = (double *) malloc(config.local_size * sizeof(double));	
	double * B_matrix = (double *) malloc(config.local_size * sizeof(double));
	
	config.A = A_matrix;	
	config.B = B_matrix;

	/* Set fileview of process to respective matrix block */
	MPI_Offset viewOffset = config.world_rank * config.local_size * sizeof(double) + sizeof(int) * 2;
	
	MPI_File_set_view(config.A_file, viewOffset, MPI_DOUBLE, config.block, "native", MPI_INFO_NULL);	
	MPI_File_set_view(config.B_file, viewOffset, MPI_DOUBLE, config.block, "native", MPI_INFO_NULL);

	/* Collective read blocks from files */
	MPI_File_read_all(config.A_file, config.A, config.local_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_read_all(config.B_file, config.B, config.local_size, MPI_DOUBLE, MPI_STATUS_IGNORE);

 	//JUST FOR TESTING:	
	int i;
	for (i = 0; i < config.local_size; i++) {
		printf("Rank: %d, col_rank: %d, row_rank: %d, value: %f \n", config.world_rank, config.col_rank, config.row_rank, config.A[i]); 
		
	}
	

	/* Close data source files */
}

void cleanup_matmul()
{
	/* Rank zero writes header specifying dim of result matrix C */

	/* Set fileview of process to respective matrix block with header offset */

	/* Collective write and close file */

	/* Cleanup */
}

void compute_fox()
{

	/* Compute source and target for verticle shift of B blocks */
	int i; 
	for (i = 0; i < config.dim[0]; i++) {
		/* Diag + i broadcast block A horizontally and use A_tmp to preserve own local A */

		/* dgemm with blocks */
		
		/* Shfting block B upwards and receive from process below */

	}
}
