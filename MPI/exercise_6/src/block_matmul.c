#include "block_matmul.h"
//#include "mkl.h"
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

	/* Process local matrix dim */
	int local_dims[2];
	int local_size;
};

struct Config config;

void init_matmul(char *A_file, char *B_file, char *outfile)
{
/* Copy output file name to configuration */
	config.outfile = outfile;
	MPI_Comm_rank(MPI_COMM_WORLD, &(config.world_rank)); 	
	/* Get matrix size header */

		MPI_File_open(MPI_COMM_WORLD, A_file, (MPI_MODE_RDWR | MPI_MODE_CREATE), MPI_INFO_NULL, &config.A_file);
		MPI_File_read(config.A_file, config.A_dims, 2, MPI_INT, MPI_STATUS_IGNORE);	

		MPI_File_open(MPI_COMM_WORLD, B_file, (MPI_MODE_RDWR | MPI_MODE_CREATE), MPI_INFO_NULL, &config.B_file);
		MPI_File_read(config.B_file, config.B_dims, 2, MPI_INT, MPI_STATUS_IGNORE);	
		config.matrix_size = config.A_dims[0] * config.A_dims[0];
	
//	printf("rank: %d ", config.world_rank);

//	if(config.world_rank == 0) {
//		double test[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
//		MPI_File_read(config.A_file, test, 16, MPI_DOUBLE, MPI_STATUS_IGNORE);
//		int j; 
//		for(j = 0; j < 16; j++) {
//			printf("rank; %d, j: %d, value: %f..... \n", config.world_rank, j, test[j]);
//		}
//	}
	
	
	/* Broadcast global matrix sizes */
	MPI_Bcast(&config.matrix_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&(config.A_dims), 2, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&(config.B_dims), 2, MPI_INT, 0, MPI_COMM_WORLD);
	/* Set dim of tiles relative to the number of processes as NxN where N=sqrt(world_size) */
	MPI_Comm_size(MPI_COMM_WORLD, &(config.world_size));
	MPI_Dims_create(config.world_size, 2, config.dim);
	//printf("config.dim = (%d, %d)\n", config.dim[0], config.dim[1]);
	//printf("World_size: %d", config.world_size);	
	/* Verify dim of A and B matches for matul and both are square*/
	if(!(config.A_dims[0] == config.B_dims[0] && config.A_dims[0] == config.A_dims[1] && config.B_dims[0] == config.B_dims[1])) {
		printf("WRONG DIMENSIONS\n");
	
	}
	/* Create Cart communicator for NxN processes */
	int period[2] = {1, 1};
	MPI_Cart_create(MPI_COMM_WORLD, 2, config.dim, period, 1, &(config.grid_comm));
	MPI_Comm_rank(config.grid_comm, &(config.grid_rank));
	/* Sub div cart communicator to N row communicator */
//	MPI_Comm_rank(MPI_COMM_WORLD, &(config.world_rank));
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
	int start[2] = {config.row_rank, config.col_rank};
	MPI_Type_create_subarray(2, config.A_dims, config.local_dims, start, MPI_ORDER_C, MPI_DOUBLE, &(config.block));
	MPI_Type_commit(&(config.block));

	/* Create data array to load actual block matrix data */
	double * A_matrix = (double *) malloc(config.local_size * sizeof(double));	
	double * B_matrix = (double *) malloc(config.local_size * sizeof(double));
	double * C_matrix = (double *) malloc(config.local_size * sizeof(double));
	double * A_tmp_matrix = (double *) malloc(config.local_size * sizeof(double));
	config.A_tmp = A_tmp_matrix;
	config.C = C_matrix;
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
	printf("Rank: %d", config.grid_rank);
	for (i = 0; i < 4; i++) {
		printf(" %f", config.B[i]); 
		
	}
	printf("\n");

	/* Close data source files */
        MPI_File_close(&(config.A_file));
        MPI_File_close(&(config.B_file));
}

void cleanup_matmul()
{
        /* Rank zero writes header specifying dim of result matrix C */
        if (config.world_rank == 0) {
                MPI_File_open(MPI_COMM_WORLD, config.outfile, (MPI_MODE_RDWR | MPI_MODE_CREATE), MPI_INFO_NULL, &(config.C_file));
                MPI_File_write(config.C_file, config.A_dims, 2, MPI_INT, MPI_STATUS_IGNORE);
        }

        /* Set fileview of process to respective matrix block with header offset */
        MPI_Offset viewOffset = (config.local_size*config.dim[0] + config.local_dims[0] * config.row_rank) * config.local_size * sizeof(double) + sizeof(int) * 2;
        MPI_File_set_view(config.C_file, viewOffset, MPI_DOUBLE, config.block, "native", MPI_INFO_NULL);

        /* Collective write and close file */
        MPI_File_write_all(config.C_file, config.C, config.local_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
        MPI_File_close(&config.C_file);
        /* Cleanup */

}


void compute_fox()
{

	/* Compute source and target for verticle shift of B blocks */
	int src;
	int dest;
	int diag;
	double alpha = 1.0;
	double beta = 0.0;
	MPI_Request req;
	//int root;
	int i;	
	//printf("rank: %d, local_size: %d", config.grid_rank, config.local_size);
	for (i = 0; i < config.dim[0]; i++) {
		
		/* Diag + i broadcast block A horizontally and use A_tmp to preserve own local A */
		//if((config.row_rank+i)%config.dim[0] == config.col_rank)
		//{
		//	printf("(%d, %d)\n", config.row_rank, config.col_rank);
		//	config.A_tmp = config.A;	
		//}
		for(diag = 0; diag < config.dim[0]; diag++){
			if((diag+i)%config.dim[0] == config.col_rank){
				config.A_tmp = config.A;
			}		
		
			if(config.col_rank == diag){
	//			printf("diag + i=%d\n",diag+i);	
				MPI_Bcast(config.A_tmp, config.local_size, MPI_DOUBLE, (diag+i)%config.dim[0], config.row_comm);
			}
		}
		/* dgemm with blocks */
		//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, &config.local_dims[0], &config.local_dims[0], &config.local_dims[0]);//, &alpha, config.A_tmp, &config.local_dims[0], config.B, &config.local_dims[1], &beta, config.C, &config.local_dims[0]);	
		
		int i, j, k;
 
		for (i = 0 ; i < config.local_size ; i+=config.local_dims[0]) {
			for (j = 0 ; j < config.local_dims[1] ; j++) {
				for (k = 0 ; k < config.local_dims[1] ; k++) {
					config.C[i+j] += config.A[i+k] * config.B[k*config.local_dims[0]+j];
					//matrix_r[i][j] += matrix_a[i][config.local_dims*k] * matrix_b[k][j];
				}
			}
		}
		
		/* Shfting block B upwards and receive from process below */
		MPI_Cart_shift(config.col_comm, 0, 1, &src, &dest);
		MPI_Send(config.B, config.local_size, MPI_DOUBLE, src, 0, config.col_comm);
		MPI_Recv(config.B, config.local_size, MPI_DOUBLE, dest, 0, config.col_comm, MPI_STATUS_IGNORE);
		if(config.grid_rank == 1 || config.grid_rank == 4 || config.grid_rank == 7){
			//printf("Current shift: %d", i+1);
		//	printf("src: %d, dest %d, rank: %d, col_rank: %d\n", src, dest, config.grid_rank, config.col_rank);
			printf("printing C for grid %d ", config.grid_rank);
			int j;
			for(j = 0; j < 4; j++){
			printf("%f ", config.C[j]);
			}
			printf("\n");
			
		}
		//MPI_Barrier(MPI_COMM_WORLD);
	}
}
