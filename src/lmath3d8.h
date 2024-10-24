//
// Created by huacheng on 2021/1/8.
//

#ifndef SALECVTSREADER_LMATH3D8_H
#define SALECVTSREADER_LMATH3D8_H

#include <math.h>
#include <assert.h>
#define BETA 0.57735026918962576451

#define NGI3d 8
#define NIpV 8
#define NIpB 4

#define DIM 3
#define CARTESIAN 1
#if CARTESIAN
#define detJV detJVc
#else
#define detJV detJVs
#endif

// #define NSPE   9    // number of sub-element in every big element
// value of 1/sqrt(3) used in position of gaussian integral
//extern const double BETA;
// points of Gaussian integral for 2*2*2 in 3d
extern const double GIPS3d[NGI3d][3];
// weights of Gaussian integral for 2*2*2 in 3d
extern const double GIWS3d[NGI3d];
// definition for the interpolation function of boundary/volume
// interpolation function for boundary within 4 points
typedef double (*Interpolation)(const double *_x);
// interpolation function for Volume with in 8 points
extern Interpolation IpV[NIpV];
extern Interpolation IpV_X[NIpV][3];
void XgIpV(double Xg[], const double Xi[][DIM],const double xl[]);
void XgIpV_X(double Xg[], const double Xi[][DIM],const double xl[], const int _d);
double detJVc(const double Xi[][DIM],const double xl[]);
double DeriveVolume(const double Xi[][DIM]);
double DeriveVolumeAverage(const double Di[], const double Xi[][DIM]);
double DataIpV(const double Di[],const double xl[]);
double detJVs(const double Xi[][DIM],const double xl[]);

// manipulate points function
void Cross(const double * _a,const double * _b, double * _c);
void CrossAddition(const double * _a,const double * _b,double * _c);
double Dot(const double * _a,const double * _b);
void Addition(double * _a, const double * _b);
void ScalerAddition(double * _a,const double * _b, double _p);
void Copy(const double _a[], double _b[]);
void Scaler(double _a[], double _p);
void Zero(double _a[]);
double Determinant(const double tM[][DIM]);
double MatrixTime(double _a[], const double _m[][DIM], const double _b[]);
void Normalization(double p[]);
double Length(const double p[]);
double Contraction(const double _a[][DIM], const double _b[][DIM]);
double Max(double a, double b);
double Min(double a, double b);
double Wind(double x,double a, double b);

/// some function on surface
extern Interpolation IpB[NIpB];
extern Interpolation IpB_X[NIpB][2];
#define NGI2d 4
extern const double GIPS2d[NGI2d][2];
extern const double GIWS2d[NGI2d];
void DeriveArea(double Xi[][DIM], double a[]);
void XgIpB(double Xg[], double Xi[][DIM], const double xl[]);

#endif //SALECVTSREADER_LMATH3D8_H

