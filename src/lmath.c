//
// Created by huachengli on 9/1/24.
//

#include "lmath.h"
#include <math.h>

void VecAddF(float *x, float *y, float p, int n)
{
    assert(n >= 1);
    for(int k=0;k<n;++k) x[k] += y[k]*p;
}

void VecZeroF(float *x, int n)
{
    assert(n >= 1);
    for(int k=0;k<n;++k) x[k] = 0.f;
}

int VecMaxArgF(const float *x, int n)
{
    assert(n >= 1);
    int rst = 0;
    for(int k=0;k<n;++k)
    {
        if(x[rst] < x[k]) rst = k;
    }
    return rst;
}

int VecMinArgF(const float *x, int n)
{
    assert(n >= 1);
    int rst = 0;
    for(int k=0;k<n;++k)
    {
        if(x[rst] > x[k]) rst = k;
    }
    return rst;
}

float VecMaxF(const float *x, int n)
{
    assert(n >= 1);
    float rst = x[0];
    for(int k=0;k<n;++k)
    {
        if(rst < x[k]) rst = x[k];
    }
    return rst;
}

float VecMinF(const float *x, int n)
{
    assert(n >= 1);
    float rst = x[0];
    for(int k=0;k<n;++k)
    {
        if(rst > x[k]) rst = x[k];
    }
    return rst;
}

float VecDisF(const float *x, const float * y, int n)
{
    assert(n>=1);
    float rst = 0;
    for(int k=0;k<n;++k)
    {
        rst += (x[k] - y[k])*(x[k] - y[k]);
    }
    return sqrtf(rst);
}



double VecScaler(double *x, double p,int n)
{
    for(int k=0;k<n;++k) x[k] *= p;
}

double VecLen(const double *x, int n)
{
    double rst = 0.;
    for(int k=0; k<n; ++k) rst += x[k]*x[k];
    return sqrt(rst);
}

float VecLenF(const float *x, int n)
{
    float rst = 0.f;
    for(int k=0; k<n; ++k) rst += x[k]*x[k];
    return sqrtf(rst);
}

double VecDot(const double *x, const double *y, int n)
{
    double rst = 0.;
    for(int k=0; k<n; ++k) rst += x[k]*y[k];
    return rst;
}

float VecDotF(const float *x, const float *y, int n)
{
    float rst = 0.f;
    for(int k=0; k<n; ++k) rst += x[k]*y[k];
    return rst;
}

void VecCross(double *z, const double *x, const double *y, int n) {
    assert(3 == n);
    z[0] = x[1] * y[2] - x[2] * y[1];
    z[1] = x[2] * y[0] - x[0] * y[2];
    z[2] = x[0] * y[1] - x[1] * y[0];
}

void VecLinear(double *z, const double *x, double px, const double *y, double py, int n)
{
    for(int k=0;k<n;++k) z[k] = x[k]*px + y[k]*py;
}

void VecAdd(double *x , const double *y, double p,int n)
{
    assert(n >= 1);
    for(int k=0;k<n;++k) x[k] += y[k]*p;
}

void VecNormalize(double *x, int n)
{
    double Ln = VecLen(x,n);
    VecScaler(x, 1.0/Ln, n);
}


void VecD2F(float *y, const double * x, int n)
{
    for(int k=0;k<n;++k)
        y[k] = (float) x[k];
}


double Clock(int i){
    static struct timeval start = {0, 0};

    if(i == 0)
    {
        gettimeofday(&start, NULL);
        return 0.;
    }
    else
    {
        struct timeval current;
        gettimeofday(&current,NULL);
        double elapsed = ( - start.tv_sec + current.tv_sec) + ( - start.tv_usec + current.tv_usec)/1000000.0;
        return elapsed;
    }
}