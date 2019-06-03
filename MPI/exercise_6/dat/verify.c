#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
	if (argc != 4) {
		printf("Error\n");
		exit(1);
	}

	FILE *file = fopen(argv[1], "rb");
	int dims [2];
	fread(dims, sizeof(int), 2, file);
	printf("(%d,%d)\n", dims[0], dims[1]);

	FILE *file1 = fopen(argv[3], "rb");
	fread(dims, sizeof(int), 2, file1);

	double *C = (double*)malloc(sizeof(double) * (dims[0] * dims[1]));
	double *C_ans = (double*)malloc(sizeof(double) * (dims[0] * dims[1]));
	
	double *A = (double*)malloc(sizeof(double) * (dims[0] * dims[1]));

	fread(C, sizeof(double), dims[0] * dims[1], file);
	fread(A, sizeof(double), dims[0] * dims[1], file1);
	int sizePerProc = 32;
	int offset = 0;
	//for(int i = 0; i < dims[1]; i++){
	//	printf("%d ", A[i]);

	//}
	//printf("\n");

	fclose(file);

	int ans_dims[2];
	file = fopen(argv[2], "rb");
	fread(ans_dims, sizeof(int), 2, file);
	fread(C_ans, sizeof(double), dims[0] * dims[1], file);

#pragma omp parallel for
	for (int i = 0; i < ans_dims[0] * ans_dims[1]; i++) {
		//printf("our: %f vs %f--------%d\n", C[i], C_ans[i], i);
		if (fabs(C[i] - C_ans[i]) > 1E-8) {
			printf("Error %f %f\n", C[i], C_ans[i]);
			exit(1);
		}
	}	
}
