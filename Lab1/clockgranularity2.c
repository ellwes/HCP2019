#include <stdint.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

double mysecond();

//Code from https://github.com/skuhl/sys-prog-examples/blob/master/simple-examples/rdtsc.c :
// Use the preprocessor so we know definitively that these are placed inline
#define RDTSC_START()            \
	__asm__ volatile("CPUID\n\t" \
	                 "RDTSC\n\t" \
	                 "mov %%edx, %0\n\t" \
	                 "mov %%eax, %1\n\t" \
	                 : "=r" (start_hi), "=r" (start_lo) \
	                 :: "%rax", "%rbx", "%rcx", "%rdx");

#define RDTSC_STOP()              \
	__asm__ volatile("RDTSCP\n\t" \
	                 "mov %%edx, %0\n\t" \
	                 "mov %%eax, %1\n\t" \
	                 "CPUID\n\t" \
	                 : "=r" (end_hi), "=r" (end_lo) \
	                 :: "%rax", "%rbx", "%rcx", "%rdx");

// Returns the elapsed time given the high and low bits of the start and stop time.
uint64_t elapsed(uint32_t start_hi, uint32_t start_lo,
                 uint32_t end_hi,   uint32_t end_lo)
{
	uint64_t start = (((uint64_t)start_hi) << 32) | start_lo;
	uint64_t end   = (((uint64_t)end_hi)   << 32) | end_lo;
	return end-start;
}


int main(int argc, char *argv[])
{
  double *ticks, mindist, maxdist, dist;
  int    i, n, maxn;

  if (argc > 0) n = atoi(argv[1]);

  ticks = (double *)malloc(n * sizeof(double));
  
  int sum;
  int e; 
  for ( e = 0; e < n; e++){
	sum = sum + mysecond();
  }
  printf("Cold start total cycles:%d\n", sum);

  
  sum = 0;
  for (j = 0; j<20; j++){

  	for (i=0; i<n; i++) {
     
     		ticks[i] = mysecond();
     sum = sum + ticks[i];
  }
}

  printf("Total cycles: %d\n", sum);


  
  mindist = 1e10;
  maxdist = 0;
  for (i=1; i<n; i++) {
    dist = ticks[i];
    if (dist > 0) {
      if (dist < mindist) mindist = dist;
      if (dist > maxdist) maxdist = dist;
    }
  }

  printf("# min dist = %.2e, max dist = %.2e",
	 mindist, maxdist);
  for (i=1; i<n; i++) {
    printf("%d\t%.2e\n", i, ticks[i]);
  }

  return 0;
}
// function with timer                                                             
double mysecond()
{
	uint32_t start_hi=0, start_lo=0; 
	uint32_t   end_hi=0,   end_lo=0;

	RDTSC_START();
	RDTSC_STOP();
	return(elapsed(start_hi, start_lo, end_hi, end_lo));
	 
}
