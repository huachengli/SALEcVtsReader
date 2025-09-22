//
// Created by huach on 26/11/2023.
//

#ifndef SALECVTSREADER_VTPTRACER_H
#define SALECVTSREADER_VTPTRACER_H
#include "Utility.h"
#include "VtpReader.h"
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
    float * ejecta_U;
    float * ejecta_V;
    float * ejecta_t;
    float * ejecta_x;
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
    int ejecta_num; // num of ejecta recorded in GridTracer
    float z_low;
    float z_up; // the threshold for detect ejectas
} GridTracer;

typedef struct TracerFileCollectDef
{
    VtpFile ** vtp;
    int NoF;
    char Name[MaxStrLen];
    char prefix[MaxStrLen];
    int step;
    int ne;
    int nex;
    int ney;
    int nez;
} VtpTracerCollect;

typedef struct name2data
{
    char name[MaxStrLen];
    VtpData * data;
} Name2VtpData;

VtpTracerCollect * OpenVtpTracerCollect(const char * _prefix, int _nof);
VtpTracerCollect * OpenSALEcTracerCollect(InputFile * ifp, int step);
VtpTracerCollect * FlushVtpTracerCollect(GridTracer * gtf,const char * _prefix, int _nof);
int CloseVtpTracerCollect(VtpTracerCollect * _vtc);
int ShowBriefVtpFile(VtpFile * vfp, FILE * fp);
int ShowBriefVtpColleect(VtpTracerCollect * vtc, FILE * fp);
VtpFile * SALEcVtpCollectMatFilter(VtpTracerCollect * vtc);
VtpFile * SALEcVtpCollectPosFilter(VtpTracerCollect * vtc);
VtpFile * SALEcVtpCollectPosFuncFilter(VtpTracerCollect * vtc, int (*tf)(const double*,const double*),const double * ctx);
VtpFile * SALEcVtpCollectPosFuncFilter2(VtpTracerCollect * vtc, int (*tf)(const double*,void *),void * ctx);
VtpFile * VtpGetMelting(VtpFile * in);
VtpFile * VtpGetConnect(VtpFile * in, VtpTracerCollect * vtc,SALEcData * ref);
int InitGridTracer(GridTracer * gtf,InputFile * ifp);
int LoadGridTxtFile(GridTracer * gtf,const char * fname);
int FlushGridTracerFromVtp(GridTracer * gtf, VtpFile * vfp);
int FlushGridTracerFromVtpCollect(GridTracer * gtf, VtpTracerCollect * tvtcp);
int WriteGridTracer(GridTracer * gtf, const char * vts_name);
int ExportGridTracerF32Bin(GridTracer * gtf, const char * binprefix);
unsigned long find_vtpfield(const char * _src, VtpFile * vfp);
VtpFile * duplicate_vtp(VtpFile * in, unsigned long n);
int vtp_add_field(VtpFile * in, const char * name, int noc);
int copy_vtp_k(VtpFile * x, unsigned long px, VtpFile * y, unsigned long  py);

typedef struct Ring
{
    double * pts;
    int npts;
    int * chords;
    int nchord;
} ring;
ring * load_ring_txt(const char * fname);
int check_in_ring(double * p, ring * _r);
int check_in_polygon(double * x, double * pl, int n);
void clean_ring(ring * _r);
#endif //SALECVTSREADER_VTPTRACER_H
