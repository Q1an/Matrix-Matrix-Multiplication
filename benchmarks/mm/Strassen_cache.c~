//
//  Strassen_cache.c
//  Strassen-cache
//
//  Created by Qian YingDong on 5/2/15.
//  Copyright (c) 2015 Qian YingDong. All rights reserved.
//

#include "Strassen_cache.h"


void naive(int size,int Xpitch, const int X[],int Ypitch, const int Y[],int Zpitch, int Z[]){
    int i;
    int j;
    int k;
    for(i=0; i<size; i++){
        for(j=0; j<size; j++){
            Z[i*Zpitch+j]=0;
            for(k=0; k<size; k++){
				Z[i*Zpitch+j]+=X[i*Xpitch+k]*Y[Ypitch*k+j];
            }
        }
    }
}



//Here used the cache tiling strategy to calculate the multiplication of certain size
//Variable block can be adjusted to optimize
// note that block here should be smaller than
void mmult(int size,int Xpitch, const int X[],int Ypitch, const int Y[],int Zpitch, int Z[]) {
    int i=0,j=0,k=0,jj=0,kk=0,Xi,Zi;
    for (jj=0;jj<size;jj=jj+block){
        for(kk=0;kk<size;kk=kk+block){
            for (i = 0; i < size; i++){
                Xi=i*Xpitch;
                Zi=Zpitch*i;
                for (j = jj; j < jj+block; j++) {
                    int sum = 0;
                    for (k = kk; k < kk+block; k++)
                        sum += X[Xi + k]*Y[k*Ypitch + j];
                    Z[Zi + j] += sum;
                }
            }
        }
    }
    
}


// S = X + Y
void madd(int size,int Xpitch, const int X[],int Ypitch, const int Y[],int Spitch, int S[]) {
    int i,j,Si,Xi,Yi;
    for (i = 0; i < size; i++){
        Si=i*Spitch;
        Xi=i*Xpitch;
        Yi=i*Ypitch;
        for (j = 0; j < size; j++)
            S[Si + j] = X[Xi + j] + Y[Yi + j];
    }
}


// S = X - Y
void msub(int size,int Xpitch, const int X[],int Ypitch, const int Y[], int Spitch, int S[]) {
    int i,j,Si,Xi,Yi;
    for (i = 0; i < size; i++){
        Si=i*Spitch;
        Xi=i*Xpitch;
        Yi=i*Ypitch;
        for (j = 0; j < size; j++)
            S[Si + j] = X[Xi + j] - Y[Yi + j];
    }
}


// Strassen algorithm for matrix multiplication
// Theoretical Runtime is O(N^2.807)
// But due to the cache size, we need to shift the size threshold to optimize
void strassen(int size,int Xpitch, const int X[],int Ypitch, const int Y[],int Zpitch, int Z[]) {

    
    if (size <= leafsize) { //change the threshold to optimize the performance
        int i,j;
        for (i = 0; i < size; i++){
            for (j = 0; j < size; j++){
                Z[i*Zpitch + j] = 0;
            }
        }
        
        mmult(size, Xpitch, X, Ypitch, Y, Zpitch, Z);
        return;
    }
    
    const int n = size/2;      // size of sub-matrices (size is always the power of 2)
    
    const int *A = X;    // Divide X into 4 sub-matrix
    const int *B = X + n;
    const int *C = X + n*Xpitch;
    const int *D = C + n;
    
    const int *E = Y;    // Divide Y into 4 sub-matrix
    const int *F = Y + n;
    const int *G = Y + n*Ypitch;
    const int *H = G + n;
    
    int *P[7];   // allocate temp matrices off heap
    const int sz = n*n*sizeof( int);
    int i=0;
    for (i = 0; i < 7; i++)
        P[i] = ( int *) malloc(sz);
    int *T = ( int *) malloc(sz);
    int *U = ( int *) malloc(sz);
    
    // P0 = A*(F - H);
    msub(n, Ypitch, F, Ypitch, H, n, T);
    strassen(n, Xpitch, A, n, T, n, P[0]);
    
    // P1 = (A + B)*H
    madd(n, Xpitch, A, Xpitch, B, n, T);
    strassen(n, n, T, Ypitch, H, n, P[1]);
    
    // P2 = (C + D)*E
    madd(n, Xpitch, C, Xpitch, D, n, T);
    strassen(n, n, T, Ypitch, E, n, P[2]);
    
    // P3 = D*(G - E);
    msub(n, Ypitch, G, Ypitch, E, n, T);
    strassen(n, Xpitch, D, n, T, n, P[3]);
    
    // P4 = (A + D)*(E + H)
    madd(n, Xpitch, A, Xpitch, D, n, T);
    madd(n, Ypitch, E, Ypitch, H, n, U);
    strassen(n, n, T, n, U, n, P[4]);
    
    // P5 = (B - D)*(G + H)
    msub(n, Xpitch, B, Xpitch, D, n, T);
    madd(n, Ypitch, G, Ypitch, H, n, U);
    strassen(n, n, T, n, U, n, P[5]);
    
    // P6 = (A - C)*(E + F)
    msub(n, Xpitch, A, Xpitch, C, n, T);
    madd(n, Ypitch, E, Ypitch, F, n, U);
    strassen(n, n, T, n, U, n, P[6]);
    
    // Z upper left = (P3 + P4) + (P5 - P1)
    madd(n, n, P[4], n, P[3], n, T);
    msub(n, n, P[5], n, P[1], n, U);
    madd(n, n, T, n, U, Zpitch, Z);
    
    // Z lower left = P2 + P3
    madd(n, n, P[2], n, P[3], Zpitch, Z + n*Zpitch);
    
    // Z upper right = P0 + P1
    madd(n, n, P[0], n, P[1], Zpitch, Z + n);
    
    // Z lower right = (P0 + P4) - (P2 + P6)
    madd(n, n, P[0], n, P[4], n, T);
    madd(n, n, P[2], n, P[6], n, U);
    msub(n, n, T, n, U, Zpitch, Z + n*(Zpitch + 1));
    
    free(U);  
    free(T);
    for (i = 6; i >= 0; i--)
        free(P[i]);
}


//I have written another version of cache tiling strategy, however based on the test on my own computer, I think 1000 is still not a big enough size to benefit from this strategy.


