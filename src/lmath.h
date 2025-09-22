//
// Created by huachengli on 9/1/24.
//

#ifndef SALECVTSREADER_LMATH_H
#define SALECVTSREADER_LMATH_H
#include "lmath3d8.h"
#include <assert.h>
#include <sys/time.h>
#include <stdio.h>
void VecAddF(float *x, float *y, float p, int n);
void VecZeroF(float *x, int n);
void VecZero(double *x, int n);
int VecMaxArgF(const float *x, int n);
int VecMinArgF(const float *x, int n);
float VecMaxF(const float *x, int n);
float VecMinF(const float *x, int n);
float VecDisF(const float *x, const float * y, int n);
double VecDis(const double *x, const double *y, int n);
float VecLenF(const float *x, int n);
double Clock(int i);

double VecLen(const double *x, int n);
double VecDot(const double *x, const double *y, int n);
float VecDotF(const float *x, const float *y, int n);
void VecCross(double *z, const double *x, const double *y, int n);
void VecNormalize(double *x, int n);
void VecLinear(double *z, const double *x, double px, const double *y, double py, int n);
void VecAdd(double *x , const double *y, double p,int n);
void VecScale(double *x, double px, int n);
void VecCopy(double * y, double * x, int n);
void VecD2F(float *y, const double * x, int n);
void VecF2D(double *y, const float * x, int n);
void VecRotate(double * x, double ro, double fo, int n);
void Contour(double **d, int ilb, int iub, int jlb, int jub,
             double *x, double *y, double z,
             double ** con_pts, int * num_pts);
void Contour_pts_sort(double * con_pts, int * num_pts, int ** pseg, double tol);
#endif //SALECVTSREADER_LMATH_H
