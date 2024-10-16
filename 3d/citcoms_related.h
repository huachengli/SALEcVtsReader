//
// Created by huachengli on 8/26/24.
//

#ifndef SALECVTSREADER_CITCOMS_RELATED_H
#define SALECVTSREADER_CITCOMS_RELATED_H

#include "Utility.h"
#include "VtkWriter.h"

typedef struct CitcomsTemImpl
{
    int nox;
    int noy;
    int noz;
    int nno;
    int ncaps;
    double * x;
    double * y;
    double * z;
    double * X[3];
    double * data;

    double * crust;
    double * dump;
} citcoms_temp_dump;

typedef struct CitcomsTracerImpl
{
    int ncaps;
    int num_basic_q;
    int num_extra_q;
    int nflavors;
    int itc;
    int * ntracers;

    double ** basicq;
    double ** extraq;
    int ** ielement;
} citcoms_tracer_dump;

typedef struct TracerMixedImpl
{
    int id[3]; // (procsid, capid, index)
    int mask[4]; // (blockid, xi, yi, zi)
    VTSDATAFLOAT vof[4];
    VTSDATAFLOAT d0;
    int sid;
    float flavor;
} tracer_mixed;

typedef struct CitcomsTracerMixedImpl
{
    tracer_mixed * data;
    unsigned int len;
    unsigned int len_alloc;
} citcoms_tracer_mixed;

typedef struct CitcomsDumpImpl
{
    int nproc;
    int nproc_surf;
    int nprocx;
    int nprocy;
    int nprocz;
    char temp_prefix[4096];
    char tracer_prefix[4096];
    double TransformR;
    double pad_step;
    citcoms_temp_dump * temp;
    citcoms_tracer_dump * tracer;
    citcoms_tracer_mixed * mtracer;
} citcoms_dump;

typedef struct TracerIJKFinder
{
    int nox;
    int noy;
    int noz;
    int elx;
    int ely;
    int elz;

    double  A[3];
    double  B[3];
    double  C[3];
    double  D[3];

    double P0[3];
    double P1[3];
    double P2[3];

    double Q0[3];
    double Q1[3];
    double Q2[3];

    float P0f[3];
    float P1f[3];
    float P2f[3];

    float Q0f[3];
    float Q1f[3];
    float Q2f[3];
} tracer_finder;

int load_citcoms_temp_dump(citcoms_temp_dump * _ctd, const char * _fname);
int load_citcoms_tracer_dump(citcoms_tracer_dump * _ctd, const char * _fname);
int load_citcoms_dump(citcoms_dump * _cd, InputFile * ifp);

int clean_citcoms_temp_dump(citcoms_temp_dump * _ctd);
int clean_citcoms_tracer_dump(citcoms_tracer_dump * _ctd);
int clean_citcoms_dump(citcoms_dump * x);
citcoms_dump * InitCitcomsDump(InputFile * ifp);
int CloseCitcomsDump(citcoms_dump * x);
SALEcData * CrInitSALEcData(InputFile * ifp);
SALEcData * CrInitSALEcData_ref(InputFile * ifp);
void CrCloseSALEcData(SALEcData * _sdata);

int UpdateCitcomsTempDump(citcoms_dump * _cd, SALEcData * _sdata, SALEcData * _rdata);
int CheckCitcomsTracerDump(citcoms_dump * _cd);
int UpdateCitcomsTracerDump(citcoms_dump * _cd, SALEcData * _sdata);
int UpdateCitcomsDump(citcoms_dump * _cdp, SALEcData * _sdata, SALEcData * _rdata);
int SALEcGetCData(SALEcData * _sdata, int fId, VTSDATAFLOAT * _pos, VTSDATAFLOAT * _data);
int SALEcGetCDataN(SALEcData * _sdata, int *fId, int length, VTSDATAFLOAT * _pos, VTSDATAFLOAT * _data, int * mask);
int SALEcGetCDataMask(SALEcData * _sdata, VTSDATAFLOAT * _pos, int * _id);
int WriteCitcomsDump(citcoms_dump * _cdp);

void citcoms_tracer_dump_vtp(citcoms_tracer_dump * _ctd, const char * name);
void citcoms_tracer_dump_pvtp(citcoms_dump * _cdp, const char * name);
int write_citcoms_temp_dump(citcoms_temp_dump * _ctd, const char * fname);
int write_citcoms_tracer_dump(citcoms_tracer_dump * _ctd, const char * fname);
int citcoms_tracer_mixed_init(citcoms_tracer_mixed * _ctm);
void citcoms_tracer_mixed_clean(citcoms_tracer_mixed * _ctm);
int citcoms_tracer_mixed_push(citcoms_tracer_mixed * _ctm, tracer_mixed * x);
int tracer_mixed_cmp(const void * _a, const void * _b);
int tracer_mixed_vofcmp(const void * _a, const void * _b);
void citcoms_tracer_mixed_export(citcoms_tracer_mixed * _ctm, citcoms_dump * _cd, int sid, int len,const char * name);

int citcoms_offset(int i, int j, int k, int nx, int ny, int nz);
void tracer_finder_init(tracer_finder * _tf,citcoms_dump * _cd,int p[4]);
int citcoms_check_tracer_element(citcoms_dump * _cd);
double solve_local(double * x,double * v2, double * v1, double * v0);
float solve_local_f(float * x,float * v2, float * v1, float * v0);
void set_projection_axis(double * n2, double * n1, double * n0, double * A, double * B, double *C, double *D);
#endif //SALECVTSREADER_CITCOMS_RELATED_H
