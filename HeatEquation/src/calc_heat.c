#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "calc_heat.h"

struct Config {
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

double get_element(int row, int col, double * matrix) {
	double res = matrix[row * config.local_dims[0] + col];
	return res;
}

void set_element(int row, int col, double element, double * matrix) {
	matrix[row * config.local_dims[0] + col] = element;
}

void init_solver(const double K, int size, const double heat_start, const double cold_start, const double time, char * outfile) {
 	config.outfile = outfile;
	
	/* Init time */	
	config.time = 0; //How should this be handeled between processes?? 
	
	config.K = K;	

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
	
	/* Create subarray datatype for local matrix tile */	
	int start[2];
	MPI_Cart_coords(config.grid_comm, config.grid_rank, 2, start);
	start[0] *= config.write_local_dims[0];
	start[1] *= config.write_local_dims[0];


	MPI_Type_create_subarray(2, config.matrix_dim, config.write_local_dims, start,  MPI_ORDER_C, MPI_DOUBLE, &config.block);

 	/* Create data array */
	double * m = (double *) calloc(config.local_dims[0] * config.local_dims[1], sizeof(double));
	double * m_tmp = (double *) calloc(config.local_dims[0] * config.local_dims[1], sizeof(double));

	/*spatial time*/
	config.hx = 1/size;
	config.hy = 1/size;

	/* Init data array with heat_start and cold_start */
	for (int i = 1; i < config.local_dims[0]-1; i++) {
		for (int j = 1; j < config.local_dims[1]-1; j++) {
			if (config.grid_rank < config.grid_dim[0]) {
				// Set upper to heat:
				if(i == 1) {
					set_element(i, j, heat_start, m);	
					set_element(i, j, heat_start, m_tmp);
				} else if (get_element(i, j, m) != heat_start) { // in case we should set both a col and row
					set_element(i, j, cold_start, m);
					set_element(i, j, cold_start, m_tmp);
				}
			}
			if (config.grid_rank % config.grid_dim[1] == 0 ) {	
				// Set right to heat:				
				if(j == config.local_dims[1]) {
					set_element(i, j, heat_start, m);	
					set_element(i, j, heat_start, m_tmp);
				} else if (get_element(i, j, m) != heat_start) { // in case we should set both a col and row
					set_element(i, j, cold_start, m);
					set_element(i, j, cold_start, m_tmp);
				}
			}
			if (config.grid_rank % config.grid_dim[1] == 1) {
				//Set left to heat
				if(j == 1) {
					set_element(i, j, heat_start, m);
					set_element(i, j, heat_start, m_tmp);				
				} else if (get_element(i, j, m) != heat_start) {
					set_element(i, j, cold_start, m);
					set_element(i, j, cold_start, m_tmp);
				}
			}
			
			if (config.grid_rank > (config.grid_dim[0] - 1) * (config.grid_dim[1])) {
				//Set lower row
				if (i == config.local_dims[0] - 1) {
					set_element(i, j, heat_start, m);
					set_element(i, j, heat_start, m_tmp);				
				} else if (get_element(i, j, m) != heat_start) {
					set_element(i, j, cold_start, m);
					set_element(i, j, cold_start, m);
				}
			}
			
			if (get_element(i, j, m) != heat_start) {
				set_element(i, j, cold_start, m);
				set_element(i, j, cold_start, m_tmp);
			}
		}
	}

	printf("NU JÄVLAR\n");	
	config.local_matrix = m; 
	config.local_matrix_tmp = m_tmp; 
	printf("DEHÄR DÅÅÅ????\n");

	if (config.grid_rank == 4) {
		print_matrix(config.local_matrix);
	}

}

void print_matrix( double * matrix) {
	for (int i = 0; i < config.local_dims[0]; i++) {
		for (int j = 0; j < config.local_dims[1]; j++) {
			printf("%f ", get_element(i, j, matrix));	
			if( j == config.local_dims[1] -1 ) {
				printf("\n");
			} 
		}
	}

}

void cleanup() {
 /* Set file view and write to file - använd grid_comm?  */
 



 /* Close file */
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
		for(j = 0; j < config.local_dims[(int) i/2] - 2; j++){ //i = 0/1 are rows, 2/4 are cols
			//send buffer depending if it's a row or col ghost cell.
			//send_buffer = i < 2 ? config.matrix[offset_row[i]][offset_col[i]+j] : config.matrix[offset_row[i]+j][offset_col[i]]
			//set send_buffer and receiving row/col to correct values depending on current row that is being sent.
			if(i < 2){
				send_buffer = get_element(offset_row[i], offset_col[i] + j, config.local_matrix);
				//send_buffer = config.local_matrix[offset_row[i]][offset_col[i] + j];
				rec_row = abs(offset_row[i] - 1);
				rec_col = offset_col[i] + j;
			}
			else{
				send_buffer = get_element(offset_row[i] + j, offset_col[i], config.local_matrix);
				//send_buffer = config.local_matrix[offset_row[i] + j][offset_col[i]];
				rec_row = offset_row[i] + j;
				rec_col = abs(offset_col[i] - 1);
			}
			
			if(dir[i] != MPI_PROC_NULL){ //send & recv.
				MPI_Send(&send_buffer, 1, MPI_DOUBLE, dir[i], config.grid_rank, config.grid_comm);
				MPI_Recv(&recv_buffer, 1, MPI_DOUBLE, dir[i], dir[i],config.grid_comm, MPI_STATUS_IGNORE);
			}
			else{
				//Process tries to send to something outside of the grid space.
				// set recv_buffer to send buffer..
				recv_buffer = send_buffer;
			}
			//Update the matrix value
			set_element(rec_row, rec_col, recv_buffer, config.local_matrix);
			//config.local_matrix[rec_row][rec_col] = recv_buffer;					
		}
	}


	int j; 
	//Loop through all the inner points of the matrix. The outer points will either be updated by other processes or remain constant.
	for(i = 1; i < config.local_dims[0] - 1; i++){
		for(j = 1; j < config.local_dims[1] - 1; j++){
			//Compute the stepping. Where to define K.
			double u_11 = get_element(i, j, config.local_matrix);
			double u_01 = get_element(i - 1, j, config.local_matrix);
			double u_21 = get_element(i + 1, j, config.local_matrix);
			double u_10 = get_element(i, j - 1, config.local_matrix);
			double u_12 = get_element(i, j + 1, config.local_matrix); 			

			double tmp = u_11 + config.K * config.d_time * ((u_21 - 2 * u_11 + u_01) / config.hx / config.hx  +  (u_12 -2 * u_11 + u_10) / config.hy / config.hy);
			set_element(i, j, tmp, config.local_matrix_tmp);
 
			//config.local_matrix_tmp[i][j] = config.local_matrix[i][j] + config.K * config.d_time * ((config.local_matrix[i + 1][j] - 2 * config.local_matrix[i][j] + config.local_matrix[i - 1][j]) / config.hx / config.hx + (config.local_matrix[i][j + 1] - 2 * config.local_matrix[i][j] + config.local_matrix[i][j -1]) / config.hy / config.hy);
		}
	} 

	/*update the config.matrix*/
	for(i = 1; i < config.local_dims[0] - 1; i++){
		for(j = 1; j < config.local_dims[0] -1; j++){
			double tmp = get_element(i, j, config.local_matrix_tmp);
			set_element(i, j, tmp, config.local_matrix);
		}

	}
}


void calc_heat() {
 
 for(double t = config.time; t > 0; t--) { //Maybe change this to the other way around? (t = 0; t < config.time; t++)
 	/* Broadcast values to close cells */
 	/* Calculate current theta */
	step();
	MPI_Barrier(config.grid_rank);
 }

}



