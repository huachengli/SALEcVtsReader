//
// Created by huachengli on 9/1/24.
//

#include "lmath.h"

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