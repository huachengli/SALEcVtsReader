//
// Created by huacheng on 2021/1/8.
//

#include "lmath3d8.h"
const int X = 0;
const int Y = 1;
const int Z = 2;

const double GIPS3d[NGI3d][3] = {
        {-BETA, -BETA, -BETA} , { -BETA, -BETA, BETA}, {-BETA, BETA, -BETA} , { -BETA, BETA, BETA},
        { BETA, -BETA, -BETA} , {  BETA, -BETA, BETA}, { BETA, BETA, -BETA} , {  BETA, BETA, BETA}};
const double GIWS3d[NGI3d] = {1.0, 1.0, 1.0, 1.0,1.0, 1.0, 1.0, 1.0};

void XgIpV(double Xg[], const double Xi[][DIM],const double xl[])
{
    // Xg is the coordinate global
    // Xi is the Vertexes of Volume
    // xl is the local coordinate
    Zero(Xg);
    for(int i=0;i<NIpV;i++)
        ScalerAddition(Xg,Xi[i],IpV[i](xl));
}
double DataIpV(const double Di[],const double xl[])
{
    double tmp = 0.0;
    for(int i=0;i<NIpV;i++)
        tmp += Di[i] * IpV[i](xl);
    return tmp;
}

//double __xIpV(const double data[], , const double Di[NIpV][],const double xl[])

void XgIpV_X(double Xg[], const double Xi[][DIM],const double xl[], const int _d)
{
    // _d is in {X,Y,Z}
    Zero(Xg);
    for(int i=0;i<NIpV;i++)
        ScalerAddition(Xg,Xi[i],IpV_X[i][_d](xl));
}
double detJVc(const double Xi[][DIM],const double xl[])
{
    // the Determinant of jacobi matrix
    // Xi is the Vertexes of Volume
    // xl is the local coordinate
    // tM is the transpose of the jacobi matrix
    double tM[DIM][DIM] = {0.0};
    for(int k=0;k<DIM;k++)
    {
        XgIpV_X(tM[k],Xi,xl,k);
    }
    return Determinant(tM);
}

double detJVs(const double Xi[][DIM],const double xl[])
{
    // in citcom
    // x[DIM] = {colatitude,longitude,radius};
    double tM[DIM][DIM] = {0.0};
    double tXg[DIM];
    for(int k=0;k<DIM;k++)
    {
        XgIpV_X(tM[k],Xi,xl,k);
    }
    XgIpV(tXg,Xi,xl);
    return Determinant(tM) * tXg[2] * tXg[2] * fabs(cos(tXg[0]));
    //    return Determinant(tM);
}

double DeriveVolume(const double Xi[][DIM])
{
    double tV = 0.0;
    for(int k=0;k<NGI3d;k++)
    {
        tV += GIWS3d[k] * detJV(Xi,GIPS3d[k]);
    }
    return tV;
}

double DeriveVolumeAverage(const double Di[], const double Xi[][DIM])
{
    double tVxData = 0.0;
    for(int k=0;k<NIpV;k++)
        tVxData += GIWS3d[k] * detJV(Xi,GIPS3d[k]) * DataIpV(Di,GIPS3d[k]);
    return tVxData;
}



// define some function to mainipulate points
void Cross(const double * _a,const double * _b, double * _c)
{
    // _c = _a x _b
    _c[X] = _a[Y]*_b[Z] - _a[Z]*_b[Y];
    _c[Y] = _a[Z]*_b[X] - _a[X]*_b[Z];
    _c[Z] = _a[X]*_b[Y] - _a[Y]*_b[X];
}
void CrossAddition(const double * _a,const double * _b,double * _c)
{
    // _c += _a x _b
    _c[X] += _a[Y]*_b[Z] - _a[Z]*_b[Y];
    _c[Y] += _a[Z]*_b[X] - _a[X]*_b[Z];
    _c[Z] += _a[X]*_b[Y] - _a[Y]*_b[X];
}

double Dot(const double * _a,const double * _b)
{
    return _a[X]*_b[X] + _a[Y]*_b[Y] + _a[Z]*_b[Z];
}
void Addition(double * _a, const double * _b)
{
    // _a = _a + _b
    for(int k=0;k<DIM;k++)
        _a[k] = _a[k] + _b[k];
}
void ScalerAddition(double * _a,const double * _b, double _p)
{
    // _a = _a + p * _b
    for(int k=0;k<DIM;k++)
        _a[k] += _p*_b[k];
}
void Copy(const double _a[], double _b[])
{
    // _b = _a
    for(int k=0;k<DIM;k++)
        _b[k] = _a[k];
}
void Scaler(double _a[], double _p)
{
    for(int k=0;k<DIM;k++)
        _a[k] *= _p;
}
void Zero(double _a[])
{
    for(int k=0;k<DIM;k++)
        _a[k] = 0.0;
}
double Determinant(const double tM[][DIM])
{
    // operator det of 3*3 matrix
    double tmp = tM[X][X] * (tM[Y][Y]*tM[Z][Z] - tM[Z][Y]*tM[Y][Z])
                 - tM[Y][X] * (tM[X][Y]*tM[Z][Z] - tM[X][Z]*tM[Z][Y])
                 + tM[Z][X] * (tM[X][Y]*tM[Y][Z] - tM[X][Z]*tM[Y][Y]);
    return tmp;
}
double MatrixTime(double _a[], const double _m[][DIM], const double _b[])
{
    // _a = _m * _b
    for(int k=0;k<DIM;k++)
    {
        _a[k] = Dot(_m[k],_b);
    }
}

void Normalization(double p[])
{
    double tmp = 0;
    for(int k=0;k<DIM;k++)
    {
        tmp += p[k] * p[k];
    }
    for(int k=0;k<DIM;k++)
    {
        p[k] = p[k]/sqrt(tmp);
    }
}

double Length(const double p[])
{
    double tmp = 0.0;
    for(int k=0;k<DIM;k++)
    {
        tmp += p[k] * p[k];
    }
    return sqrt(tmp);
}

double Contraction(const double _a[][DIM], const double _b[][DIM])
{
    double tmp = 0.0;
    for(int k=0;k<DIM;k++)
        for(int j=0;j<DIM;j++)
        {
            tmp += _a[k][j] * _b[k][j];
        }
    return tmp;
}

double Max(double a, double b)
{
    if(a > b)
        return a;
    else
        return b;
}

double Min(double a, double b)
{
    if(a < b)
        return a;
    else
        return b;
}

double Wind(double x,double a, double b)
{
    if(x < a)
        return a;
    else if(x > b)
        return b;
    else
        return x;
}

// functions used in 8 points interpolate
// see the function array IpV and IpV_B
double Ip3d0(const double * _x)
{
    return 0.125*(1 + _x[0])*(1 + _x[1])*(1 - _x[2]);
}
double Ip3d1(const double * _x)
{
    return 0.125*(1 - _x[0])*(1 + _x[1])*(1 - _x[2]);
}
double Ip3d2(const double * _x)
{
    return 0.125*(1 - _x[0])*(1 - _x[1])*(1 - _x[2]);
}
double Ip3d3(const double * _x)
{
    return 0.125*(1 + _x[0])*(1 - _x[1])*(1 - _x[2]);
}
double Ip3d4(const double * _x)
{
    return 0.125*(1 + _x[0])*(1 + _x[1])*(1 + _x[2]);
}
double Ip3d5(const double * _x)
{
    return 0.125*(1 - _x[0])*(1 + _x[1])*(1 + _x[2]);
}
double Ip3d6(const double * _x)
{
    return 0.125*(1 - _x[0])*(1 - _x[1])*(1 + _x[2]);
}
double Ip3d7(const double * _x)
{
    return 0.125*(1 + _x[0])*(1 - _x[1])*(1 + _x[2]);
}


double Ip3d0_X1(const double * _x)
{
    return 0.125*(1 + _x[1])*(1 - _x[2]);
}
double Ip3d1_X1(const double * _x)
{
    return -0.125*(1 + _x[1])*(1 - _x[2]);
}
double Ip3d2_X1(const double * _x)
{
    return -0.125*(1 - _x[1])*(1 - _x[2]);
}
double Ip3d3_X1(const double * _x)
{
    return 0.125*(1 - _x[1])*(1 - _x[2]);
}
double Ip3d4_X1(const double * _x)
{
    return 0.125*(1 + _x[1])*(1 + _x[2]);
}
double Ip3d5_X1(const double * _x)
{
    return -0.125*(1 + _x[1])*(1 + _x[2]);
}
double Ip3d6_X1(const double * _x)
{
    return -0.125*(1 - _x[1])*(1 + _x[2]);
}
double Ip3d7_X1(const double * _x)
{
    return 0.125*(1 - _x[1])*(1 + _x[2]);
}

double Ip3d0_X2(const double * _x)
{
    return 0.125*(1 + _x[0])*(1 - _x[2]);
}
double Ip3d1_X2(const double * _x)
{
    return 0.125*(1 - _x[0])*(1 - _x[2]);
}
double Ip3d2_X2(const double * _x)
{
    return -0.125*(1 - _x[0])*(1 - _x[2]);
}
double Ip3d3_X2(const double * _x)
{
    return -0.125*(1 + _x[0])*(1 - _x[2]);
}
double Ip3d4_X2(const double * _x)
{
    return 0.125*(1 + _x[0])*(1 + _x[2]);
}
double Ip3d5_X2(const double * _x)
{
    return 0.125*(1 - _x[0])*(1 + _x[2]);
}
double Ip3d6_X2(const double * _x)
{
    return -0.125*(1 - _x[0])*(1 + _x[2]);
}
double Ip3d7_X2(const double * _x)
{
    return -0.125*(1 + _x[0])*(1 + _x[2]);
}

double Ip3d0_X3(const double * _x)
{
    return -0.125*(1 + _x[0])*(1 + _x[1]);
}
double Ip3d1_X3(const double * _x)
{
    return -0.125*(1 - _x[0])*(1 + _x[1]);
}
double Ip3d2_X3(const double * _x)
{
    return -0.125*(1 - _x[0])*(1 - _x[1]);
}
double Ip3d3_X3(const double * _x)
{
    return -0.125*(1 + _x[0])*(1 - _x[1]);
}
double Ip3d4_X3(const double * _x)
{
    return 0.125*(1 + _x[0])*(1 + _x[1]);
}
double Ip3d5_X3(const double * _x)
{
    return 0.125*(1 - _x[0])*(1 + _x[1]);
}
double Ip3d6_X3(const double * _x)
{
    return 0.125*(1 - _x[0])*(1 - _x[1]);
}
double Ip3d7_X3(const double * _x)
{
    return 0.125*(1 + _x[0])*(1 - _x[1]);
}

Interpolation IpV[NIpV] = {Ip3d0,Ip3d1,Ip3d2,Ip3d3,Ip3d4,Ip3d5,Ip3d6,Ip3d7};
Interpolation IpV_X[NIpV][3] = {
        { Ip3d0_X1, Ip3d0_X2, Ip3d0_X3},
        { Ip3d1_X1, Ip3d1_X2, Ip3d1_X3},
        { Ip3d2_X1, Ip3d2_X2, Ip3d2_X3},
        { Ip3d3_X1, Ip3d3_X2, Ip3d3_X3},
        { Ip3d4_X1, Ip3d4_X2, Ip3d4_X3},
        { Ip3d5_X1, Ip3d5_X2, Ip3d5_X3},
        { Ip3d6_X1, Ip3d6_X2, Ip3d6_X3},
        { Ip3d7_X1, Ip3d7_X2, Ip3d7_X3},
};

/// interpolation function on surface
double Ip2d0(const double *_x)
{
    return 0.25 * (1 + _x[0]) * (1 + _x[1]);
}

double Ip2d0_X1(const double *_x)
{
    return 0.25 * (1 + _x[1]);
}

double Ip2d0_X2(const double *_x)
{
    return 0.25 * (1 + _x[0]);
}

double Ip2d1(const double *_x)
{
    return 0.25 * (1 - _x[0]) * (1 + _x[1]);
}

double Ip2d1_X1(const double *_x)
{
    return -0.25 * (1 + _x[1]);
}

double Ip2d1_X2(const double *_x)
{
    return 0.25 * (1 - _x[0]);
}

double Ip2d2(const double *_x)
{
    return 0.25 * (1 - _x[0]) * (1 - _x[1]);
}

double Ip2d2_X1(const double *_x)
{
    return -0.25 * (1 - _x[1]);
}

double Ip2d2_X2(const double *_x)
{
    return -0.25 * (1 - _x[0]);
}

double Ip2d3(const double *_x)
{
    return 0.25 * (1 + _x[0]) * (1 - _x[1]);
}

double Ip2d3_X1(const double *_x)
{
    return 0.25 * (1 - _x[1]);
}

double Ip2d3_X2(const double *_x)
{
    return -0.25 * (1 + _x[0]);
}

Interpolation IpB[NIpB] = {Ip2d0, Ip2d1, Ip2d2, Ip2d3};
Interpolation IpB_X[NIpB][2] = {
        {Ip2d0_X1, Ip2d0_X2},
        {Ip2d1_X1, Ip2d1_X2},
        {Ip2d2_X1, Ip2d2_X2},
        {Ip2d3_X1, Ip2d3_X2},
};

void XgIpB(double Xg[], double Xi[][DIM], const double xl[])
{
    // Xg is the coordinate global
    // Xi is the Vertexes of Volume
    // xl is the local coordinate
    Zero(Xg);
    for (int i = 0; i < NIpB; i++)
    {
        ScalerAddition(Xg, Xi[i], IpB[i](xl));
    }
}

void XgIpB_X(double Xg[], double Xi[][DIM], const double xl[], const int _d)
{
    Zero(Xg);
    for (int i = 0; i < NIpB; i++)
    {
        ScalerAddition(Xg, Xi[i], IpB_X[i][_d](xl));
    }
}

void DeriveArea(double Xi[][DIM], double a[])
{
    // calculate the norm vector of a face
    double tArea[DIM] = {0.0, 0.0, 0.0};
    for (int k = 0; k < NGI2d; k++)
    {
        double tA[3] = {0.0, 0.0, 0.0};
        double tX1[3] = {0.0, 0.0, 0.0};
        double tX2[3] = {0.0, 0.0, 0.0};

        XgIpB_X(tX1, Xi, GIPS2d[k], X);
        XgIpB_X(tX2, Xi, GIPS2d[k], Y);

        Cross(tX1, tX2, tA);
        ScalerAddition(tArea, tA, GIWS2d[k]);
    }
    Copy(tArea, a);
}

const double GIPS2d[NGI2d][2] = {
        {-BETA, -BETA},
        {-BETA, BETA},
        {BETA,  -BETA},
        {BETA,  BETA}};
const double GIWS2d[NGI2d] = {1.0, 1.0, 1.0, 1.0};