//
// Created by huach on 26/11/2023.
//

#ifndef SALECVTSREADER_VTPTRACER_H
#define SALECVTSREADER_VTPTRACER_H

// let's map tracer info back to grid
typedef struct GridTracerDef
{
    float * ipos; // initial position
    float * cpos; // current position
    float * vel;
    float * pre;
    float * tem;
    float * den;
    float * mpre;
    float * mtem;
    float * ejecta_X;
    float * ejecta_T;
    float * id;
    float * matid;
    float * mask;
    int nx;
    int ny;
    int nx2;
    int ny2;
    int len;
    int stripe;
    int nvtp;
    int step;
    int maxstep;
    double dtsave;
    double t0;
    double gz;
    char prefix[MaxStrLen];
} GridTracer;

typedef struct TracerFileCollectDef
{
    VtpFile ** vtp;
    int NoF;
    const char Name[MaxStrLen];
} VtpTracerCollect;

typedef struct name2data
{
    char name[MaxStrLen];
    VtpData * data;
} Name2VtpData;

VtpTracerCollect * OpenVtpTracerCollect(const char * _prefix, int _nof);
int CloseVtpTracerCollect(VtpTracerCollect * _vtc);
int InitGridTracer(GridTracer * gtf,InputFile * ifp);
int LoadGridTxtFile(GridTracer * gtf,const char * fname);
int FlushGridTracerFromVtp(GridTracer * gtf, VtpFile * vfp);
int FlushGridTracerFromVtpCollect(GridTracer * gtf, VtpTracerCollect * tvtcp);
int WriteGridTracer(GridTracer * gtf, const char * vts_name);
#endif //SALECVTSREADER_VTPTRACER_H
