#include <stdio.h>
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
	double ** local_matrix;
	double ** local_matrix_tmp;
	int local_dims[2]; // Dimension including ghost cells
	int write_local_dims[2]; // When writing to file, we do not want to use ghost cells. Therefor this is local dims without counting the ghost cells
	int block;

	/* Time */	
	double time;
	double d_time; //Time step
} config;


void init_solver(const double K, int size, const double heat_start, const double cold_start, const double time, char * outfile) {
 	config.outfile = outfile;
	
	/* Init time */	
	config.time = 0; //How should this be handeled between processes?? 
	
 	/* Create Cart communicator for NxN processes */
	MPI_Comm_size(MPI_COMM_WORLD, &config.world_size);
	MPI_Dims_create(config.world_size, 2, config.grid_dim);
 
	int period[2] = {0, 0};
	MPI_Cart_create(MPI_COMM_WORLD, 2, config.grid_dim, period, 1, &config.grid_comm);
	MPI_Comm_rank(config.grid_comm, &config.grid_rank);


	/* Setup sizes of local matrix tiles */
        config.local_dims[0] = config.matrix_dim[0] / config.grid_dim[0] + 2;
	config.local_dims[1] = config.matrix_dim[1] / config.grid_dim[1] + 2;


	config.write_local_dims[0] = config.local_dims[0] - 2; 
	config.write_local_dims[1] = config.local_dims[1] - 2; 

	/* Create subarray datatype for local matrix tile */	
	int start[2];
	MPI_Cart_coords(config.grid_comm, config.grid_rank, 2, start);
	start[0] *= config.write_local_dims[0];
	start[1] *= config.write_local_dims[0];
	MPI_Type_create_subarray(2, config.grid_dim, config.write_local_dims, start,  MPI_ORDER_C, MPI_DOUBLE, &config.block);

 	/* Create data array */
	double m[config.local_dims[0]][config.local_dims[1]];
	double m_tmp[config.local_dims[0]][config.local_dims[1]];




	/* Init data array with heat_start and cold_start */
	for (int i = 1; i < config.local_dims[0]-1; i++) {
		for (int j = 1; j < config.local_dims[1]-1; j++) {
			if (config.grid_rank < config.grid_dim[0]) {
				// Set upper to heat:				
				if(i == 1) {
					m[i][j] = heat_start;	
					m_tmp[i][j] = heat_start;	
				} else if (m[i][j] != heat_start) { // in case we should set both a col and row
					m[i][j] = cold_start;
					m_tmp[i][j] = cold_start;
				}
			}
			if (config.grid_rank % config.grid_dim[1] == 0 ) {	
				// Set right to heat:				
				if(j == config.local_dims[1]) {
					m[i][j] = heat_start;	
					m_tmp[i][j] = heat_start;	
				} else if (m[i][j] != heat_start) { // in case we should set both a col and row
					m[i][j] = cold_start;
					m_tmp[i][j] = cold_start;
				}
			}
			
			if (config.grid_rank % config.grid_dim[1] == 2) {
				//Set left to heat
				if(j == 1) {
					m[i][j] = heat_start;
					m_tmp[i][j] = heat_start;
				} else if (m[i][j] != heat_start) {
					m[i][j] = cold_start;
					m_tmp[i][j] = cold_start;
				}
			}
			
			if (config.grid_rank > (config.grid_dim[0] - 1) * (config.grid_dim[1])) {
				//Set lower row
				if (i == config.local_dims[0] - 1) {
					m[i][j] = heat_start; 
					m_tmp[i][j] = heat_start; 
				} else if (m[i][j] != heat_start) {
					m[i][j] = cold_start;
					m_tmp[i][j] = cold_start; 
				}
			}
			
			if (m[i][j] != heat_start) {
				m[i][j] = cold_start;
				m_tmp[i][j] = cold_start;
			}
		}
	}

	config.local_matrix = &m; 
	config.local_matrix_tmp = &m_tmp;
}




void print_matrix() {
	for (int i = 0; i < config.local_dims[0]; i++) {
		for (int j = 0; j < config.local_dims[1]; j++) {
			printf("%f ", config.local_matrix[i][j]);
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





void calc_heat() {
 
 for(double t = config.time; t > 0; t--) {
 	/* Broadcast values to close cells */
 	/* Calculate current theta */
 }

}
