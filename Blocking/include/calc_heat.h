#ifndef __CANNON_H__
#define __CANNON_H__

#include "mpi.h"

void init_solver(const double K, int size, const double heat_start, const double cold_start, const double time, char * outfile);
void calc_heat();
void cleanup();
void print_matrix(double * matrix, int row, int col);
void print_file();

#endif
