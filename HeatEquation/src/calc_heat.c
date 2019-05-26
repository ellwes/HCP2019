

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

	/*spatial*/
	double hx;
	double hy

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

	/*spatial time*/
	config.hx = 1/size;
	config.hy = 1/size;

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
 
 for(double t = config.time; t > 0; t--) { //Maybe change this to the other way around? (t = 0; t < config.time; t++)
 	/* Broadcast values to close cells */
 	/* Calculate current theta */
	step();
	MPI_Barrier();
 }

}

void step() {
	
	//	up, down, left, right
	int dir[4] = {-1, -1, -1, -1};
	int offset_row[4] = {0, config.local_dims[0], 1, 1}; //row offset start
	int offset_col[4] = {1, 1, 0, config.local_dims[1]}; //col offset start
	int rec_row, rec_col;	//receiving row & column pos.	

	//Observe that this might be wrong and more Cart_shifts may be needed.
	MPI_Cart_shift(config.grid_comm, 0, 1, &dir[0], &dir[1]);	//shift up/down
	MPI_Cart_shift(config.grid_comm, 1, 1, &dir[2], &dir[3]);	//shift left/right
	printf("I am %d, up/down is %d/%d. And my left/right is %d/%d", config.grid_rank, dir[0], dir[1], dir[2], dir[3]);

	int i;
	double send_buffer, recv_buffer;
	for(i = 0; i < 4; i++){
		//Get the element to send and put it in buffer:
		int j;
		for(j = 0; j < config.local_size-2; j++){
			//send buffer depending if it's a row or col ghost cell.
			//send_buffer = i < 2 ? config.matrix[offset_row[i]][offset_col[i]+j] : config.matrix[offset_row[i]+j][offset_col[i]]
			//set send_buffer and receiving row/col to correct values depending on current row that is being sent.
			if(i < 2){
				send_buffer = config.matrix[offset_row[i]][offset_col[i] + j];
				rec_row = abs(offset_row[i] - 1);
				rec_col = offset_col[i] + j;
			}
			else{
				send_buffer = config.matrix[offset_row[i] + j][offset_col[i]];
				rec_row = offset_row[i] + j;
				rec_col = abs(offset_col[i] - 1);
			}
			
			if(dir[i] != MPI_PROC_NULL){ //send & recv.
				MPI_Send(&send_buffer, 1, MPI_DOUBLE, dir[i], config.grid_rank, config.grid_comm);
				MPI_Recv(&recv_buffer, 1, MPI_DOUBLE, dir[i], config.grid_comm, MPI_STATUS_IGNORE);
			}
			else{
				//Process tries to send to something outside of the grid space.
				// set recv_buffer to send buffer..
				recv_buffer = send_buffer;
			}
			//Update the matrix value
			config.matrix[rec_row][rec_col] = recv_buffer;					
		}
	}


	int j; 
	//Loop through all the inner points of the matrix. The outer points will either be updated by other processes or remain constant.
	for(i = 1; i < config.local_dims[0] - 1; i++){
		for(j = 1; j < config.local_dims[1] - 1; j++){
			//Compute the stepping. Where to define K.
			config.matrix_tmp[i][j] = config.matrix[i][j] + K * config.d_time * ((config.matrix[i + 1][j] - 2 * config.matrix[i][j] + config.matrix[i - 1][j]) / config.hx / config.hx + (config.matrix[i][j + 1] - 2 * config.matrix[i][j] + config.matrix[i][j -1]) / config.hy / config.hy);
		}
	} 

	/*update the config.matrix*/
	for(i = 1; i < config.local_dims[0] - 1; i++){
		for(j = 1; j < config.local_dims[0] -1; j++){
			config.matrix[i][j] = config.matrix_tmp[i][j];
		}
	}
}





