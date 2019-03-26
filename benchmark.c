#define N 100000
#include <stdio.h>
#include <sys/time.h>

double mysecond();

int main(){
  int i, j;
  int O_n = 20;
  double t1, t2; // timers                                                         
  double a[N], b[N], c[N]; // arrays  
                                             
  // init arrays                                                                   
  for (i = 0; i < N; i++){
    a[i] = 47.0;
    b[i] = 3.1415;
  }

  //Avoid coldstart
  for (i = 0; i < N; i++) { //TODO: HOw many iterations must coldstart be? 
     c[i] = a[i]*b[i];	
  }

  // measure performance                                                           
  double sum = 0;
  double var = 0; 
  for(j = 0; j < O_n; j++){
    
  t1 = mysecond();
    for(i = 0; i < N; i++){
      c[i] = a[i]*b[i];
     
    }     
  t2 = mysecond(); 

  for(i = 0; i < N; i++) {//Is this really necessary to avoid
    var +=  c[i];
  }
  sum += t2-t1;
}

	
 printf("For avoiding smart compiling %11.8f \n", var);
 double exe_t = sum/O_n;

  printf("Execution time: %11.8f s\n", (exe_t));
  return 0;
}

// function with timer                                                             
double mysecond(){
  struct timeval tp;
  struct timezone tzp;
  int i;

  i = gettimeofday(&tp,&tzp);
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}
