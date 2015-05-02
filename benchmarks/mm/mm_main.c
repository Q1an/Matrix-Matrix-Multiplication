#include "mm_ip.h"
#include "Strassen_cache.h"
int main(){
	int i;
	int j;
	//int id;


	const int a[N][N]={
    #include "in_a_medium.txt"
	};

	const int b[N][N]={
    #include "in_b_medium.txt"
	};

	const int out_expected[N][N]={
    #include "out_c_medium.txt"
	};
    int *X, *Y, *Z;
    X = ( int*) malloc(256*256*sizeof( int));
    Y = ( int*) malloc(256*256*sizeof( int));
    Z = ( int*) malloc(256*256*sizeof( int));
	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
            X[i*256 + j] = a[i][j];
            Y[i*256 + j] = b[i][j];
		}
	}

	// HW computation
    int out[N][N];
    strassen(256, 256, X, 256, Y, 256, Z);
    for(i=0; i<N; i++){
        for(j=0; j<N; j++){
            out[i][j]=Z[i*256 + j];
        }
    }
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
