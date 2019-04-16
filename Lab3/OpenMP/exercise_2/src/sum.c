#include "sum.h"

void omp_sum(double *sum_ret)
{
	double sum = 0; 
	#pragma omp parallel
	{
		#pragma omp for
		for(int i = 0; i < size; i++){
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
		#pragma omp for	
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
    #pragma omp for
    for (int i = 0; i < size; i++){
	#pragma omp atomic
	 sum += x[i];
    }
  }	
  *sum_ret = sum;


}

void omp_local_sum(double *sum_ret)
{
	double *part_sum = (double*)calloc(32, sizeof(double));
	#pragma omp parallel
	{
		int id = omp_get_thread_num();
		#pragma omp for
		for(int i = 0; i < size; i++)
		{
			part_sum[id] += x[i];
		}
	}

	double res_sum = 0;
	for(int i = 0; i < 32; i++)
		res_sum += part_sum[i];

	*sum_ret = res_sum;
}

void omp_padded_sum(double *sum_ret)
{

        double *part_sum = (double*)calloc(64, sizeof(double));
        #pragma omp parallel
        {
                int id = omp_get_thread_num();
		#pragma omp for
                for(int i = 0; i < size; i++)
                {
                        part_sum[id*3] += x[i];
                }
        }

        double res_sum = 0;
        for(int i = 0; i < 64; i++)
                res_sum += part_sum[i];

        *sum_ret = res_sum;
}

void omp_private_sum(double *sum_ret)
{

	double res_sum = 0;
	#pragma omp parallel
        {
		double sum = 0;
                #pragma omp for
                for(int i = 0; i < size; i++)
                {
                        sum += x[i];
                }
        	res_sum += sum;
	}
        *sum_ret = res_sum;

}

void omp_reduction_sum(double *sum_ret)
{
	double sum = 0; 
	#pragma omp parallel
	{
		#pragma omp for reduction(+: sum)
		for (int i = 0; i < size; i++) {
			sum += x[i];
		}
		
	}	
	*sum_ret = sum;

}
