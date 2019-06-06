#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "calc_heat.h"

struct Config {
	MPI_File out;

	char *outfile;
	
	/* Full matrix */
	int matrix_dim[2];
	int matrix_size;
	/* Communicators */
	MPI_Comm grid_comm;
	MPI_Comm part_comm; //Not sure if this can be used
	int grid_dim[2]; //Dimensions of processes	

	/* Communicators dims and ranks */
	int world_rank, world_size, grid_rank;
	int part_rank;

	/* Matrix and local dims and local type*/
	double * local_matrix;
	double * local_matrix_tmp;
	int local_dims[2]; // Dimension including ghost cells
	int write_local_dims[2]; // When writing to file, we do not want to use ghost cells. Therefor this is local dims without counting the ghost cells
	int block;

	/*spatial*/
	double hx;
	double hy;
	
	double K;
	/* Time */	
	double time;
	double d_time; //Time step
} config;

double get_element(int row, int col, double * matrix, int dim0) {
	double res = matrix[row * dim0 + col];
	return res;
}

void set_element(int row, int col, double element, double * matrix, int dim0) {
	matrix[row * dim0 + col] = element;
}


int is_boundary_cell(int rank, int row, int col){ //function that checks if the cel at (row, cell) in the tile with rank is a boundary cell:
	/*Check if rank belongs to an outer tile*/
	if(rank < config.grid_dim[0]){ //Upper tile
		if(row == 1){
			return 1;
		}
	}
	if(rank > (config.grid_dim[0] - 1) * config.grid_dim[1] - 1){ //Bottom tile
		if(row == config.write_local_dims[0]){
			return 1;
		}
	}
	if(rank % config.grid_dim[0] == 0){//Left tile
		if(col == 1){
			return 1;
		}
	}
	if((rank + 1) % config.grid_dim[1] == 0){//Right tile
		if(col == config.write_local_dims[1]){
			return 1;
		}
	}
	return 0;	

}


void init_solver(const double K, int size, const double heat_start, const double cold_start, const double time, char * outfile) {
 	config.outfile = outfile;
	
	/* Init time */	
	config.time = time; //How should this be handeled between processes?? 	
	config.K = K;	
	config.matrix_size = size;

	/* intit world config */	
	MPI_Comm_size(MPI_COMM_WORLD, &config.world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &config.world_rank);

 	/* Create Cart communicator for NxN processes */
	MPI_Dims_create(config.world_size, 2, config.grid_dim);
 
	int period[2] = {0, 0};
	MPI_Cart_create(MPI_COMM_WORLD, 2, config.grid_dim, period, 1, &config.grid_comm);
	MPI_Comm_rank(config.grid_comm, &config.grid_rank);

	/* Setup sizes of local matrix tiles */
	config.matrix_dim[0] = sqrt(size);
	config.matrix_dim[1] = sqrt(size);	

        config.local_dims[0] = config.matrix_dim[0] / config.grid_dim[0] + 2;
	config.local_dims[1] = config.matrix_dim[1] / config.grid_dim[1] + 2;
	
	config.write_local_dims[0] = config.local_dims[0] - 2;
	config.write_local_dims[1] = config.local_dims[1] - 2;

 	/* Create data array */
	double * m = (double *) calloc(config.local_dims[0] * config.local_dims[1], sizeof(double));
	double * m_tmp = (double *) calloc(config.local_dims[0] * config.local_dims[1], sizeof(double));

	/*spatial & time steps*/
	config.hx = 1/(double)size;
	config.hy = 1/(double)size;
	config.d_time = config.hx * config.hx / 4 / config.K;

	int dim0 = config.local_dims[0];
	/* Init data array with heat_start and cold_start */
	for (int i = 1; i < config.local_dims[0]-1; i++) {
		for (int j = 1; j < config.local_dims[1]-1; j++) {
			if (config.grid_rank < config.grid_dim[0]) {
				// Set upper to heat:
				if(i == 1) {
					set_element(i, j, heat_start, m, dim0);	
					set_element(i, j, heat_start, m_tmp, dim0);
				} else if (get_element(i, j, m, dim0) != heat_start) { // in case we should set both a col and row
					set_element(i, j, cold_start, m, dim0);
					set_element(i, j, cold_start, m_tmp, dim0);
				}
			}
			
			if ((config.grid_rank + 1) % config.grid_dim[1] == 0 ) {	
				// Set right to heat:				
				if(j == config.local_dims[1]-2) {
					set_element(i, j, heat_start, m, dim0);	
					set_element(i, j, heat_start, m_tmp, dim0);
				} else if (get_element(i, j, m, dim0) != heat_start) { // in case we should set both a col and row
					set_element(i, j, cold_start, m, dim0);
					set_element(i, j, cold_start, m_tmp, dim0);
				}
			}
			
			if (config.grid_rank % config.grid_dim[1] == 0) {
				//Set left to heat
				if(j == 1) {
					set_element(i, j, heat_start, m, dim0);
					set_element(i, j, heat_start, m_tmp, dim0);				
				} else if (get_element(i, j, m, dim0) != heat_start) {
					set_element(i, j, cold_start, m, dim0);
					set_element(i, j, cold_start, m_tmp, dim0);
				}
			}
			
			if (config.grid_rank > (config.grid_dim[0] - 1) * (config.grid_dim[1])-1) {
				//Set lower row
				if (i == config.local_dims[0] - 2) {
					set_element(i, j, heat_start, m, dim0);
					set_element(i, j, heat_start, m_tmp, dim0);				
				} else if (get_element(i, j, m, dim0) != heat_start) {
					set_element(i, j, cold_start, m, dim0);
					set_element(i, j, cold_start, m, dim0);
				}
			}
			
			if (get_element(i, j, m, dim0) != heat_start) {
				set_element(i, j, cold_start, m, dim0);
				set_element(i, j, cold_start, m_tmp, dim0);
			}
		}
	}

	config.local_matrix = m; 
	config.local_matrix_tmp = m_tmp; 
}

void print_matrix( double * matrix, int rows, int cols) {
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			printf("%f ", get_element(i, j, matrix, rows));	
			if( j == config.local_dims[1] -1 ) {
				printf("\n");
			} 
		}
	}

}

void cleanup() {
	int write_matrix_size = config.write_local_dims[0]*config.write_local_dims[1];
	
	double * write_matrix = (double *) malloc( write_matrix_size * sizeof(double));
	for(int i = 1; i < config.local_dims[0] - 1; i++) {
		for(int j = 1; j < config.local_dims[1]- 1; j++) {
			double element = get_element(i, j, config.local_matrix, config.local_dims[0]);
			set_element(i-1, j-1, element, write_matrix, config.write_local_dims[0]);
		}
	}

	/* Create subarray datatype for local matrix tile */	
	int * start = (int *) malloc(2*sizeof(int));	
	MPI_Cart_coords(config.grid_comm, config.grid_rank, 2, start);
	start[0] *= config.write_local_dims[0];  
	start[1] *= config.write_local_dims[1];  

	MPI_Type_create_subarray(2, config.matrix_dim, config.write_local_dims, start,  MPI_ORDER_C, MPI_DOUBLE, &(config.block));
	MPI_Type_commit(&(config.block));	

	MPI_File_open(MPI_COMM_WORLD, config.outfile, (MPI_MODE_RDWR | MPI_MODE_CREATE), MPI_INFO_NULL, &config.out);
	//Write time to file
	if (config.world_rank == 0) {
		MPI_File_write(config.out, &config.time, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
	}	
	MPI_Offset view_offset = 1*sizeof(double);
	MPI_File_set_view(config.out, view_offset, MPI_DOUBLE, (config.block), "native",  MPI_INFO_NULL);
	MPI_File_write_all(config.out, write_matrix, write_matrix_size, MPI_DOUBLE, MPI_STATUS_IGNORE);	
	
	/* Close file */
	MPI_File_close(&config.out);
}

void print_file() {	
	MPI_File_open(MPI_COMM_WORLD, config.outfile, (MPI_MODE_RDONLY), MPI_INFO_NULL, &config.out);
	if (config.world_rank == 0) {
		printf("\nMATRIX_SIZE: %d\n", config.matrix_size);		
		double buf[config.matrix_size];
		MPI_File_read(config.out, &buf, config.matrix_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
		for (int i = 0; i < config.matrix_size; i++) {
			printf("%f ", buf[i]);
		}
	}
}


void step() {
	//	up, down, left, right
	int dir[4] = {-1, -1, -1, -1};
	
	//These arrays are the same, no matter what. Maybe put them in config and initialize in init to save time?
	int offset_row[4] = {1, config.write_local_dims[0], 1, 1};
	int offset_col[4] = {1, 1, 1, config.write_local_dims[1]};

	//top row, bottom row, left col, right col
	int recv_cell[4] = {0, config.local_dims[0] - 1, 0, config.local_dims[1] - 1};
	int rec_row, rec_col;	//receiving row & column pos.	

	//Observe that this might be wrong and more Cart_shifts may be needed.
	MPI_Cart_shift(config.grid_comm, 0, 1, &dir[0], &dir[1]);	//shift up/down
	MPI_Cart_shift(config.grid_comm, 1, 1, &dir[2], &dir[3]);	//shift left/right

	int i;
	double send_buffer;//, recv_buffer;
	
	//MPI_Request * req_send = (int*) calloc((local.write_local_dims[0])*4, sizeof(int));  
	//MPI_Request * req_recv =  (int*) calloc((local.write_local_dims[0])*4, sizeof(int));  
	MPI_Request req_send[config.write_local_dims[0]*4];
	double recv_buffer[config.write_local_dims[0]*4];

	for(i = 0; i < 4; i++){
		//Get the element to send and put it in buffer:
		int j;
		for(j = 0; j < config.local_dims[(int) i/2] - 2; j++){ //i = 0/1 are rows, 2/4 are cols
			//send buffer depending if it's a row or col ghost cell.

			int req_index = i*config.write_local_dims[0] + j;
			if(i < 2){
				send_buffer = get_element(offset_row[i], offset_col[i] + j, config.local_matrix, config.local_dims[0]);
				rec_row = recv_cell[i];
				rec_col = offset_col[i] + j;

			}
			else{
				send_buffer = get_element(offset_row[i] + j, offset_col[i], config.local_matrix, config.local_dims[0]);
				rec_row = offset_row[i] + j;
				rec_col = recv_cell[i];
			}
			
			if(dir[i] != MPI_PROC_NULL){ //send & recv.
				MPI_Isend(&send_buffer, 1, MPI_DOUBLE, dir[i], config.grid_rank, config.grid_comm, &(req_send[req_index]));
				MPI_Irecv(&(recv_buffer[req_index]), 1, MPI_DOUBLE, dir[i], dir[i], config.grid_comm, &(req_send[req_index]));
				//set_element(rec_row, rec_col, recv_buffer, config.local_matrix, config.local_dims[0]);

			} else {
				MPI_Isend(&send_buffer, 1, MPI_DOUBLE, config.grid_rank, config.grid_rank, config.grid_comm, &(req_send[req_index]));
				MPI_Irecv(&recv_buffer[req_index], 1, MPI_DOUBLE, config.grid_rank, config.grid_rank, config.grid_comm, &(req_send[req_index]));
			}					
		}
	}

	if(config.grid_rank == 0) {
		print_matrix(config.local_matrix, config.local_dims[0], config.local_dims[1]);
	}
	

	MPI_Waitall(config.write_local_dims[0] * 4, req_send, MPI_STATUS_IGNORE);
	//MPI_Waitall(config.write_local_dims[0] * 4, req_send, MPI_STATUS_IGNORE);
	
	for(int i = 0; i < config.write_local_dims[0] * 4; i++) {
		int j = i % config.write_local_dims[0];
		
		//top
		if(i < config.write_local_dims[0]){
			set_element(0, j + 1, recv_buffer[i], config.local_matrix, config.local_dims[0]);
		}
		//bottom
		else if (i-config.write_local_dims[0] < config.write_local_dims[0]){
			set_element(config.local_dims[0]-1, j+1, recv_buffer[i], config.local_matrix, config.local_dims[0]);

		}
		//left
		else if(i - 2*config.write_local_dims[0] < config.write_local_dims[0]){
			set_element(j+1, 0, recv_buffer[i], config.local_matrix, config.local_dims[0]);


		}
		//right
		else if(i - 3*config.write_local_dims[0] < config.write_local_dims[0]){
			set_element(j+1, config.local_dims[1]-1, recv_buffer[i], config.local_matrix, config.local_dims[0]);
		
		}	
	}

	int j; 
	//Loop through all the inner points of the matrix. The outer points will either be updated by other processes or remain constant.
	for(i = 1; i < config.local_dims[0] - 1; i++){
		for(j = 1; j < config.local_dims[1] - 1; j++){
			//check that it is not a boundary value:
			double tmp = 0;
			if(!is_boundary_cell(config.grid_rank, i, j)){
				//Compute the stepping. Where to define K.
				double u_11 = get_element(i, j, config.local_matrix, config.local_dims[0]);
				double u_01 = get_element(i - 1, j, config.local_matrix, config.local_dims[0]);
				double u_21 = get_element(i + 1, j, config.local_matrix, config.local_dims[0]);
				double u_10 = get_element(i, j - 1, config.local_matrix, config.local_dims[0]);
				double u_12 = get_element(i, j + 1, config.local_matrix, config.local_dims[0]); 			
				tmp = u_11 + config.K * config.d_time * ((u_21 - 2 * u_11 + u_01) / config.hx / config.hx  +  (u_12 -2 * u_11 + u_10) / config.hy / config.hy);
				
 			}
			else {
				tmp = get_element(i, j, config.local_matrix, config.local_dims[0]);
			}
			set_element(i, j, tmp, config.local_matrix, config.local_dims[0]);

		}	

	}	
	
}


void calc_heat() {
 
 for(double t = 0; t < config.time; t+=config.d_time) { //Maybe change this to the other way around? (t = 0; t < config.time; t++)
	step();
	MPI_Barrier(config.grid_comm);
 }
}


