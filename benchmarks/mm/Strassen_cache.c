//
//  Strassen_cache.c
//  Strassen-cache
//
//  Created by Qian YingDong on 5/2/15.
//  Copyright (c) 2015 Qian YingDong. All rights reserved.
//

#include "Strassen_cache.h"




//
// Classic O(N^3) square matrix multiplication.
// Z = X*Y
// All matrices are NxN and stored in row major order
// each with a specified pitch.
// The pitch is the distance (in  int's) between
// elements at (row,col) and (row+1,col).
//
void mmult(int size,int Xpitch, const int X[],int Ypitch, const int Y[],int Zpitch, int Z[]) {
    int i=0,j=0,k=0;
    for (i = 0; i < size; i++)
        for (j = 0; j < size; j++) {
            int sum = 0;
            for (k = 0; k < size; k++)
                sum += X[i*Xpitch + k]*Y[k*Ypitch + j];
            Z[i*Zpitch + j] = sum;
        }
}

//
// S = X + Y
//
void madd(int size,int Xpitch, const int X[],int Ypitch, const int Y[],int Spitch, int S[]) {
    int i,j;
    for (i = 0; i < size; i++)
        for (j = 0; j < size; j++)
            S[i*Spitch + j] = X[i*Xpitch + j] + Y[i*Ypitch + j];
}

//
// S = X - Y
//
void msub(int size,int Xpitch, const int X[],int Ypitch, const int Y[], int Spitch, int S[]) {
    int i,j;
    for (i = 0; i < size; i++)
        for (j = 0; j < size; j++)
            S[i*Spitch + j] = X[i*Xpitch + j] - Y[i*Ypitch + j];
}

//
// Volker Strassen algorithm for matrix multiplication.
// Theoretical Runtime is O(N^2.807).
// Assume NxN matrices where N is a power of two.
// Algorithm:
//   Matrices X and Y are split into four smaller
//   (N/2)x(N/2) matrices as follows:
//          _    _          _   _
//     X = | A  B |    Y = | E F |
//         | C  D |        | G H |
//          -    -          -   -
//   Then we build the following 7 matrices (requiring
//   seven (N/2)x(N/2) matrix multiplications -- this is
//   where the 2.807 = log2(7) improvement comes from):
//     P0 = A*(F - H);
//     P1 = (A + B)*H
//     P2 = (C + D)*E
//     P3 = D*(G - E);
//     P4 = (A + D)*(E + H)
//     P5 = (B - D)*(G + H)
//     P6 = (A - C)*(E + F)
//   The final result is
//        _                                            _
//   Z = | (P3 + P4) + (P5 - P1)   P0 + P1              |
//       | P2 + P3                 (P0 + P4) - (P2 + P6)|
//        -                                            -
//
void strassen(int size,int Xpitch, const int X[],int Ypitch, const int Y[],int Zpitch, int Z[]) {
    //
    // Recursive base case.
    // If matrices are 16x16 or smaller we just use
    // the conventional algorithm.
    // At what size we should switch will vary based
    // on hardware platform.
    //
    if (size <= 32) {
        mmult(size, Xpitch, X, Ypitch, Y, Zpitch, Z);
        return;
    }
    
    const int n = size/2;      // size of sub-matrices
    
    const int *A = X;    // A-D matrices embedded in X
    const int *B = X + n;
    const int *C = X + n*Xpitch;
    const int *D = C + n;
    
    const int *E = Y;    // E-H matrices embeded in Y
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
    
    free(U);  // deallocate temp matrices
    free(T);
    for (i = 6; i >= 0; i--)
        free(P[i]);
}

/*void mrand(int N, int pitch,  int M[]) {
    const int r = 10.0;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            M[i*pitch + j] = r*(2*drand48() - 1);
}*/




