#include "sum.h"

void omp_sum(double *sum_ret)
{
	double sum = 0; 
	#pragma omp parallel
	{
		for(int = 0; i < size; i++){
			sum += x[i];	
		}	
	}
	*sum_ret = sum;

}


void omp_critical_sum(double *sum_ret)
{
	double sum = 0;
	#pragma omp parallel
	{
		
		for(int i = 0; i < size; i++)
		{
			#pragma omp critical
			sum += x[i];
		}	
	}
	*sum_ret = sum;
	
}

void omp_atomic_sum(double *sum_ret)
{
  double sum = 0;
  #pragma omp parallel
  {
    for (int i = 0; i < size; i++){
	#pragma omp atomic
	 sum += x[i];
    }
  }	
  *sum_ret = sum;


}

void omp_local_sum(double *sum_ret)
{
 	double sum = 0; ii
	int num_threads = omp_get_num_threads();
	double *part_sum = (double*)calloc(num_threads, sizeof(double));
	#pragma omp parallel
	{
		int id = omp_get_thread_num();
		for(int i = id*size/num_threads; i < (id+1)*size/num_threads; i++)
		{
			part_sum[id] += x[i];
		}
	}

	double res_sum = 0;

	for(int i = 0; i < omp_get_num_threads(); i++)
		res_sum += part_sum[i];
		printf("SUM: %d\n", res_sum);

	*sum_ret = res_sum;
}

void omp_padded_sum(double *sum_ret)
{

}

void omp_private_sum(double *sum_ret)
{

}

void omp_reduction_sum(double *sum_ret)
{

}
