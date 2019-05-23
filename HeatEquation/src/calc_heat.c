

struct Config {
	char *outfile;
	
	/* Full matrix */
	double ** matrix;
	double ** matrix_tmp;
	int matrix_dim[2];
	int matrix_size;

	/* Communicators */
	MPI_Comm grid_comm;
	MPI_Comm part_comm; //Not sure if this can be used
	int grid_dim[2]; //Dimensions of processes	

	/* Communicators dims and ranks */
	int dim[2];
	int world_rank, world_size, grid_rank;
	int part_rank;

	/* Local dims and local type*/
	int local_dims[2];
	int block;

	/* Time */	
	double time;
	double d_time; //Time step
} Config;


void init_solver(const double K, int size, const double heat_start, const double cold_start, const double time, char outfile) {
 	config.outfile = outfile;
	
	/* Init time */	
	config.time = 0; //How should this be handeled between processes?? 
	
 	/* Create Cart communicator for NxN processes */
	MPI_Comm_size(MPI_COMM_WORLD, &config.world_size);
	MPI_Dims_create(config.world_size, 2, config.grid_dim);
 
	int period[2] = {0, 0};
	MPI_Cart_create(MPI_COMM_WORLD, 2, config.grid_dim, period, 1, &config.grid_comm);
	MPI_Comm_rank(config.grid_comm, &condig.grid_rank);

	/* Setup sizes of local matrix tiles */
        config.local_dim[0] = config.dim[0] / config.grid_dim[0];
	config.local_dim[1] = config.dim[1] / config.grid_dim[1];

	/* Create subarray datatype for local matrix tile */	
	int start[2];
	MPI_Cart_coords(config.grid_comm, config.grid_rank, 2, start);
	start[0] *= config.local_dims[0];
	start[1] *= config.local_dims[0];
	MPI_Type_create_subarray(2, config.A_dims, config.local_dims, start,  MPI_ORDER_C, MPI_DOUBLE, &config.block);

 	/* Create data array */
	double m[size][size];
	double m_tmp[size][size];

	config.matrix = &m; 
	config.matrix_tmp = &m_tmp;

	/* Init data array with heat_start and cold_start */
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if(j == 0 || j == size || i == 0 || i == size) {
				config.matrix[i] = heat_start; 
				config.matrix_tmp[i] = heat_start;
			} else {
				config.matrix[i] = cold_start;
				config.matrix_tmp[i] = cold_start;
			}
			
		}
	}


}




void cleanup() {
 /* Set file view and write to file - använd grid_comm?  */
 



 /* Close file */
}





void calc_heat() {
 
 for(double t = config.time; t > 0; t--) {
 	/* Broadcast values to close cells */
 	/* Calculate current theta */
 }

}
