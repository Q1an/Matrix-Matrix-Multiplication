#include "mm_ip.h"
#include "Strassen_cache.h"
#include "time.h"
const int Next_2_power = 1024;
#define N 160 //when use the 1000*1000 matrix as a test, please change this to 1000
int main(){
	int i;
	int j;



	const int a[N][N]={
    #include "in_a_medium.txt"
	};//change the file as well

	const int b[N][N]={
    #include "in_b_medium.txt"
	};

	const int out_expected[N][N]={
    #include "out_c_medium.txt"
	};
    
    int *X, *Y, *Z;
    X = (int*) malloc(Next_2_power*Next_2_power*sizeof(int));
    Y = (int*) malloc(Next_2_power*Next_2_power*sizeof(int));
    Z = (int*) malloc(Next_2_power*Next_2_power*sizeof(int));
	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
            X[i*Next_2_power + j] = a[i][j];
            Y[i*Next_2_power + j] = b[i][j];
            Z[i*Next_2_power+j]=0;
		}
	}//initialization the X,Y,Z matrix as array here

	// HW computation
    int out[N][N];
    
    clock_t start=clock();
    // use clock in case we need to test the time on other platform
    
    // The naive version has been copied here, you can test it if you want
    // naive(Next_2_power, Next_2_power, X, Next_2_power, Y, Next_2_power, Z);
    
    strassen(Next_2_power, Next_2_power, X, Next_2_power, Y, Next_2_power, Z);
    
    clock_t end=clock();
    
    printf("%fms",(double)(end-start)/CLOCKS_PER_SEC*1000);
    
    for(i=0; i<N; i++){
        for(j=0; j<N; j++){
            out[i][j]=Z[i*Next_2_power + j];
        }
    }//output to the output 2D-array 
    
	// Verification
	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
			if(out[i][j]!=out_expected[i][j]){
				printf("Verification failed! out[%d][%d]=%d, out_expected[%d][%d]=%d\n", i, j, out[i][j], i, j, out_expected[i][j]);
				return -1;
			}
		}
	}
	printf("Verification passed!\n");
    free(X);
    free(Y);
    free(Z);
	return 0;

}
