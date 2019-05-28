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

void print_matrix(double* mat)
{
	int i;
	for(i=0; i<config.local_size; i++)
	{
		printf(" %f", mat[i]);
	}
	printf("\n");
}
void init_matmul(char *A_file, char *B_file, char *outfile)
{
	/* Copy output file name to configuration */
	config.outfile = outfile;
	MPI_Comm_rank(MPI_COMM_WORLD, &(config.world_rank)); 	
	/* Get matrix size header */
	
	MPI_File_open(MPI_COMM_WORLD, A_file, (MPI_MODE_RDONLY ), MPI_INFO_NULL, &config.A_file);
	MPI_File_open(MPI_COMM_WORLD, B_file, (MPI_MODE_RDONLY ), MPI_INFO_NULL, &config.B_file);
	if(config.world_rank==0){
		MPI_File_read(config.A_file, config.A_dims, 2, MPI_INT, MPI_STATUS_IGNORE);
		MPI_File_read(config.B_file, config.B_dims, 2, MPI_INT, MPI_STATUS_IGNORE);
		config.matrix_size = config.A_dims[0] * config.A_dims[0];
	}

	/* Broadcast global matrix sizes */
	MPI_Bcast(&config.matrix_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(config.A_dims, 2, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(config.B_dims, 2, MPI_INT, 0, MPI_COMM_WORLD);
	

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
	printf("Size of big A: (%d, %d), B: (%d, %d), %d", config.A_dims[0], config.A_dims[1], config.B_dims[0], config.B_dims[1], config.grid_rank);	
	/* Setup sizes of local matrix tiles */
	config.local_dims[0] = config.A_dims[0] / config.dim[0];
	config.local_dims[1] = config.A_dims[1] / config.dim[1];
	config.local_size = config.local_dims[0] * config.local_dims[1];

	/* Create subarray datatype for local matrix tile */
	int start[2];
	MPI_Cart_coords(config.grid_comm, config.grid_rank, 2, start);
	start[0] *= config.local_dims[0];
	start[1] *= config.local_dims[1];
	MPI_Type_create_subarray(2, config.A_dims, config.local_dims, start,  MPI_ORDER_C, MPI_DOUBLE, &(config.block));
	MPI_Type_commit(&(config.block));

	/* Create data array to load actual block matrix data */
	double * A_matrix = (double *) calloc(config.local_size, sizeof(double));	
	double * B_matrix = (double *) calloc(config.local_size, sizeof(double));
	double * C_matrix = (double *) calloc(config.local_size, sizeof(double));
	double * A_tmp_matrix = (double *) calloc(config.local_size,  sizeof(double));
	config.A_tmp = A_tmp_matrix;
	config.C = C_matrix;
	config.A = A_matrix;	
	config.B = B_matrix;

	/* Set fileview of process to respective matrix block */
	MPI_Offset viewOffset = 2*sizeof(int);	

	MPI_File_set_view(config.A_file, viewOffset, MPI_DOUBLE, config.block, "native", MPI_INFO_NULL);
	MPI_File_set_view(config.B_file, viewOffset, MPI_DOUBLE, config.block, "native", MPI_INFO_NULL);
	
	/* Collective read blocks from files */
	int errA = MPI_File_read_all(config.A_file, config.A, config.local_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
	int errB = MPI_File_read_all(config.B_file, config.B, config.local_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
	printf("Finished reading %d, return_val A:%d, return_val B:%d\n", config.grid_rank, errA, errB);
	printf("Printing matrix A:\n");
	print_matrix(config.A);
	//MPI_Barrier(config.grid_comm);
	/* Close data source files */
        MPI_File_close(&(config.A_file));
      	MPI_File_close(&(config.B_file));
	printf("closed the file %d\n", config.grid_rank);
	//MPI_Barrier(config.grid_comm);
}

void cleanup_matmul()
{      
	/* Rank zero writes header specifying dim of result matrix C */
       
        MPI_File_open(MPI_COMM_WORLD, config.outfile, (MPI_MODE_RDWR | MPI_MODE_CREATE), MPI_INFO_NULL, &(config.C_file));
	if (config.world_rank == 0) {

                MPI_File_write(config.C_file, config.A_dims, 2, MPI_INT, MPI_STATUS_IGNORE);
        }

        /* Set fileview of process to respective matrix block with header offset */
        MPI_Offset viewOffset =  sizeof(int)*2;

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

	int i;
	printf("Inside compute_fox! %d, config.dim is:%d\n", config.grid_rank, config.dim[0]);
	for (i = 0; i < config.dim[0]; i++) {
		printf("updating A_tmp i:%d, %d\n", i, config.grid_rank);
		/* Diag + i broadcast block A horizontally and use A_tmp to preserve own local A */
		if((config.col_rank+i)%config.dim[0] == config.row_rank){
			int q;
			for(q=0; q<config.local_size;q++){
				config.A_tmp[q] = config.A[q];	
			}
		}

		printf("initiating broadcast i:%d, %d\n", i, config.grid_rank);
		for(diag = 0; diag < config.dim[0]; diag++){
			if((int)config.grid_rank/config.dim[0] == diag){
				//broadcast along the row
				MPI_Bcast(config.A_tmp, config.local_size, MPI_DOUBLE, (diag+i)%config.dim[0], config.row_comm);
			}
		}		
		printf("finished broadcast i:%d, %d\n", i, config.grid_rank);		
		//Matrix multiplication	
		int ii, j, k;
		for (ii = 0 ; ii < config.local_dims[0]; ii++) {
			for (j = 0 ; j < config.local_dims[0]; j++) {
				for (k = 0 ; k < config.local_dims[1] ; k++) {
					config.C[ii*config.local_dims[0]+j] += config.A_tmp[k+ii*config.local_dims[0]] * config.B[k*config.local_dims[0]+j];
				}
			}
		}
		
		/* Shfting block B upwards and receive from process below */
		MPI_Cart_shift(config.col_comm, 0, 1, &src, &dest);
		MPI_Send(config.B, config.local_size, MPI_DOUBLE, src, config.col_rank, config.col_comm);
		MPI_Recv(config.B, config.local_size, MPI_DOUBLE, dest, dest, config.col_comm, MPI_STATUS_IGNORE);	
	}
	
}
