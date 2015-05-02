//
//  Strassen_cache.h
//  Strassen-cache
//
//  Created by Qian YingDong on 5/2/15.
//  Copyright (c) 2015 Qian YingDong. All rights reserved.
//

#ifndef __Strassen_cache__Strassen_cache__
#define __Strassen_cache__Strassen_cache__

#include <stdio.h>
#include <stdlib.h>

//const int leafsize = 16;

void strassen(int size,int Xpitch, const int X[],int Ypitch, const int Y[],int Zpitch, int Z[]);
void msub(int size,int Xpitch, const int X[],int Ypitch, const int Y[],int Spitch, int S[]);
void madd(int size,int Xpitch, const int X[],int Ypitch, const int Y[],int Spitch, int S[]);
void mmult(int size,int Xpitch, const int X[], int Ypitch, const int Y[],int Zpitch, int Z[]);



#endif /* defined(__Strassen_cache__Strassen_cache__) */
