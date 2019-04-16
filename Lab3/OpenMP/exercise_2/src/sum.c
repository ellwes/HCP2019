#include "sum.h"

void omp_sum(double *sum_ret)
{
	double sum = 0;
	//double *part_sum = (double*)calloc(omp_get_num_threads(), sizeof(double));
	
	#pragma omp parallel
	{
	//	int id = omp_get_thread_num();
		int id, num_threads;
		id = omp_get_thread_num();
		num_threads = omp_get_num_threads();
		for(int i = id*size/num_threads; i < (id+1)*size/num_threads; i++)
		{
			sum += x[i];
		}
		//part_sum[id] = sum;
	}	
	
	//for(int i = 0; i < omp_get_num_threads(); i++){
	//	sum += part_sum[i];
	//}
	*sum_ret = sum;
}

void omp_critical_sum(double *sum_ret)
{
	double sum = 0;
	
	#pragma omp parallel
	{
		int id, num_threads;
		id = omp_get_thread_num();
		num_threads = omp_get_num_threads();

		
		for(int i = id; i < size; i+=num_threads)
		{
			#pragma omp critical
			sum += x[i];
		}	
	}
	
	*sum_ret = sum;
}

void omp_atomic_sum(double *sum_ret)
{

}

void omp_local_sum(double *sum_ret)
{
	//double sum = 0;
	double *part_sum = (double*)calloc(omp_get_num_threads(), sizeof(double));
	
	#pragma omp parallel
	{
		int num_threads = omp_get_num_threads();
		int id;
		id = omp_get_thread_num();
		for(int i = id*size/num_threads; i < (id+1)*size/num_threads; i++)
		{
			part_sum[id] += x[i];
		}
		//part_sum[id] = sum;
	}

	double res_sum = 0;

	for(int i = 0; i < omp_get_num_threads(); i++)
		res_sum += part_sum[i];

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
