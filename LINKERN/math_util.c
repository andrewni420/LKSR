#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "util.h"


double sigmoid(double input){
    double output = (input>=0) ? 1/(1+exp(-input)) : exp(input)/(1+exp(input));
    return output;
}

double randUnif(double a, double b){
    return a+(b-a)*((double) rand () / (double) RAND_MAX);
}

double randNormal(double mu, double sigma){
    double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;
 
  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }
 
  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);
 
  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
 
  call = !call;
 
  return (mu + sigma * (double) X1);
}

void fisherYates(int *arr, int n){
    for (int i = n-1; i >= 0; --i){
    //generate a random number [0, n-1]
    int j = rand() % (i+1);

    //swap the last element with element at random index
    int temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
    }
}

void shuffle(int *array, size_t n) {   
    struct timeval tv;
    gettimeofday(&tv, NULL);
    int usec = tv.tv_usec;
    srand48(usec);
    if (n > 1) {
        size_t i;
        for (i = n - 1; i > 0; i--) {
            size_t j = (unsigned int) (drand48()*(i+1));
            int t = array[j];
            array[j] = array[i];
            array[i] = t;
        }
    }
}

void range(int *arr, int n){
    for (int i=0;i<n;i++) arr[i]=i;
}

void logrange(double *arr, int n, double loglow, double loghigh){
    for (int i=0;i<n;i++) arr[i] = exp(loglow+(loghigh-loglow)*n);
}


#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

double quickselect(double *arr, int n, int k) {
  unsigned long i,ir,j,l,mid;
  double a,temp;

  l=0;
  ir=n-1;
  for(;;) {
    if (ir <= l+1) { 
      if (ir == l+1 && arr[ir] < arr[l]) {
	SWAP(arr[l],arr[ir]);
      }
      return arr[k];
    }
    else {
      mid=(l+ir) >> 1; 
      SWAP(arr[mid],arr[l+1]);
      if (arr[l] > arr[ir]) {
	SWAP(arr[l],arr[ir]);
      }
      if (arr[l+1] > arr[ir]) {
	SWAP(arr[l+1],arr[ir]);
      }
      if (arr[l] > arr[l+1]) {
	SWAP(arr[l],arr[l+1]);
      }
      i=l+1; 
      j=ir;
      a=arr[l+1]; 
      for (;;) { 
	do i++; while (arr[i] < a); 
	do j--; while (arr[j] > a); 
	if (j < i) break; 
	SWAP(arr[i],arr[j]);
      } 
      arr[l+1]=arr[j]; 
      arr[j]=a;
      if (j >= k) ir=j-1; 
      if (j <= k) l=i;
    }
  }
}

double median(double *arr, int n){
    return quickselect(arr, n, n/2-1);
}

double qad(double *arr, double med, int n, double q){
    double *ad = CC_SAFE_MALLOC(n,double);
    for (int i=0;i<n;i++) ad[i] = fabs(arr[i]-med);
    double median_abs_diff = quickselect(ad, n, (int)(n*q)-1);
    CC_IFFREE(ad, double);
    return median_abs_diff;
}

void normalizeMAD(double *arr, int n){
    double *copy = CC_SAFE_MALLOC(n, double);
    for (int i=0;i<n;i++) copy[i] = arr[i];
    double med = median(copy,n);
    CC_IFFREE(copy, double);
    double median_abs_diff = qad(arr, med, n, 0.5);
    // printf("median %lf mad %lf n %d n/2 %d\n",med,median_abs_diff,n,n/2);
    for (int i=0;i<n;i++) arr[i] = (arr[i]-med)/median_abs_diff;
}

void scaleQuartile(double *arr, int n, double q){
    int idx = (int)(n*q);
    double *copy = CC_SAFE_MALLOC(n, double);
    for (int i=0;i<n;i++) copy[i] = fabs(arr[i]);
    double med = qad(arr, med, n, q);
    for (int i=0;i<n;i++) arr[i] = arr[i]/med;
}

void printIntArray(int *arr, int n){
    printf("[");
    for (int i=0;i<n-1;i++) printf("%d, ",arr[i]);
    printf("%d]",arr[n-1]);
}

void printBoolArray(bool *arr, int n){
    printf("[");
    for (int i=0;i<n-1;i++) printf("%d, ",arr[i]);
    printf("%d]",arr[n-1]);
}

void printDoubleArray(double *arr, int n){
    printf("[");
    for (int i=0;i<n-1;i++) printf("%lf, ",arr[i]);
    printf("%lf]",arr[n-1]);
}

void copyDouble(double *src, double *dest, int n){
  for (int i=0;i<n;i++) dest[i]=src[i];
}

void copyInt(int *src, int *dest, int n){
  for (int i=0;i<n;i++) dest[i]=src[i];
}

void reverseInt(int *src, int n){
  int low = 0;
  int high = n-1;
  while (high>low){
    int temp = src[low];
    src[low] = src[high];
    src[high] = temp;
    low++;
    high--;
  }
}

void copyBool(bool *src, bool *dest, int n){
    for (int i=0;i<n;i++) dest[i]=src[i];
}

int compareInt (const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}