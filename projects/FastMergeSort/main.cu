#include <stdio.h>
#include <time.h>
#include <windows.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#define N 10240
#define MIN(a,b) (b<a)?b:a

void mergeNumbers(int* numbers, int start, int middle, int ending){
    int* merged = (int*) malloc((ending-start)*sizeof(int));
    int left = 0, right = 0, i = 0;
    while(left < middle - start && right < ending - middle){
        merged[i++] = numbers[start + left] < numbers[middle + right]? numbers[start + left++]: numbers[middle + right++];
    }
    while(right < ending - middle)
        merged[i++] = numbers[middle + right++];

    while(left < middle - start)
        merged[i++] = numbers[start + left++];
    memcpy(numbers+start, merged, i*sizeof(int));
    free(merged);
}

void mergeSort(int* numbers, int* sorted, int numbersLength){
    memcpy(sorted, numbers, numbersLength*sizeof(int));
    for (int i = 1; i <= numbersLength / 2 + 1; i *= 2) {
        for (int j = i; j < numbersLength; j += 2 * i) {
            mergeNumbers(sorted, j-i, j, MIN(j+i,numbersLength));
        }
    }
}

__device__ void acceleratedMergeNumbers(int* numbers, int start, int middle, int ending){
    int* merged = (int*) malloc((ending-start)*sizeof(int));
    int left = 0, right = 0, i = 0;
    while(left < middle - start && right < ending - middle){
        merged[i++] = numbers[start + left] < numbers[middle + right]? numbers[start + left++]: numbers[middle + right++];
    }
    while(right < ending - middle)
        merged[i++] = numbers[middle + right++];

    while(left < middle - start)
        merged[i++] = numbers[start + left++];
    memcpy(numbers+start, merged, i*sizeof(int));
    free(merged);
}

__global__ void acceleratedMergeSort(int* numbers, int* sorted, int numbersLength){
    memcpy(sorted, numbers, numbersLength*sizeof(int));
    for (int i = 1; i <= numbersLength / 2 + 1; i *= 2) {
        for (int j = i; j < numbersLength; j += 2 * i) {
            acceleratedMergeNumbers(sorted, j-i, j, MIN(j+i,numbersLength));
        }
    }
}


int main(int argc, char* argv[]){
    int numbers[N], sorted[N], acceleratedSorted[N];
    //int* sorted = (int*) malloc(N*sizeof(int));
    //assert(sorted != NULL);
    srand(time(NULL));
    printf("%d", RAND_MAX);
    for(int i=0;i<N;i++){
        numbers[i] = rand()%10240;
    }
    printf("\nUnsorted first 10: ");
    for(int i=0;i<10;i++) printf("%d, ", numbers[i]);

    mergeSort(numbers, sorted, N);

    printf("\nUnsorted first 10 after sorting: ");
    for(int i=0;i<10;i++) printf("%d, ", numbers[i]);
    printf("\nSorted first 10: ");
    for(int i=0;i<10;i++) printf("%d, ", sorted[i]);

    int* dev_numbers;
    int* dev_sorted;

    cudaMalloc( (void**)&dev_numbers, N * sizeof(int) );
    cudaMalloc( (void**)&dev_sorted, N * sizeof(int) );

    cudaMemcpy( dev_numbers, numbers, N*sizeof(int), cudaMemcpyHostToDevice);

    cudaDeviceSetLimit(cudaLimitMallocHeapSize, 128*1024*1024);
    acceleratedMergeSort<<<1,1>>>(dev_numbers, dev_sorted, N);

    cudaMemcpy( sorted, dev_sorted, N*sizeof(int), cudaMemcpyHostToDevice);

    printf("\nUnsorted first 10 after sorting: ");
    for(int i=0;i<10;i++) printf("%d, ", numbers[i]);
    printf("\nSorted first 10: ");
    for(int i=0;i<10;i++) printf("%d, ", acceleratedSorted[i]);
    return 0;
}
