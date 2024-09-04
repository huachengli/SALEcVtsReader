//
// Created by huachengli on 9/1/24.
//

#ifndef SALECVTSREADER_LMATH_H
#define SALECVTSREADER_LMATH_H

#include <assert.h>
#include <sys/time.h>
#include <stdio.h>
void VecAddF(float *x, float *y, float p, int n);
void VecZeroF(float *x, int n);
int VecMaxArgF(const float *x, int n);
int VecMinArgF(const float *x, int n);
float VecMaxF(const float *x, int n);
float VecMinF(const float *x, int n);
double Clock(int i);


#endif //SALECVTSREADER_LMATH_H
