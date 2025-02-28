//
// Created by huachengli on 8/26/24.
//

#ifndef SALECVTSREADER_CITCOMS_RELATED_H
#define SALECVTSREADER_CITCOMS_RELATED_H

#include "Utility.h"
#include "VtkWriter.h"

typedef struct
{
    int nox;
    int noy;
    int noz;
    int nno;
    int nel;
    int nproc;
    int nproc_surf;
    int nprocx;
    int nprocy;
    int nprocz;
    char VtsPrefix[4096];
    char OutPrefix[4096];
    char datafile[4096];
    char datapath[4096];
    int step0;
    int step1;
    int step_inc;
    VtsInfo * VSF;
    char attach[200][4096];
    int len_attach;
    char sol_liq_file[4096];

    double * gr;
    double * sol;
    double * liq;

    int ncomp;
    double * tscomp_ff; // sol/liq shift according to composition
    char melt_post_dir[4096];
} CitcomsData;


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

typedef struct CitcomsSphereDumpImpl
{
    int nox;
    int noy;
    int noz;
    int nno;
    float * pos;
    float * data; // data of cell
    float * pdata; // data of point
    int noc;
    int nel;
    float * marker;
    float * pmarker;
    float * area;
    float * vstat;
} citcoms_sphere_dump;

typedef struct CitcomsTracerDumpImpl
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

typedef struct CitcomsSphereImpl
{
    int nproc;
    int nproc_surf;
    int nprocx;
    int nprocy;
    int nprocz;
    citcoms_sphere_dump * cap;
} citcoms_sphere;


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

citcoms_sphere * init_citcoms_sphere(citcoms_dump * _cd, int noc);
int clean_citcoms_sphere(citcoms_sphere * _cs);
int write_citcoms_sphere(citcoms_sphere * _cs, CitcomsData * _cd,const char * _name);

citcoms_dump * InitCitcomsDump(InputFile * ifp);
int CloseCitcomsDump(citcoms_dump * x);
SALEcData * CrInitSALEcData(InputFile * ifp);
SALEcData * CrInitSALEcData_ref(InputFile * ifp);
void CrCloseSALEcData(SALEcData * _sdata);

int StructedGridIntf(double ** X, double * _data, int *shape, int *eid, double * res, double (*f)(double *, double *), double * ctx);
int StructedGridIntf2(float * X, float * _data, int *shape, int *eid, double * res, double (*f)(double *, double *), double * ctx);
int VIntCitcomsTempDump(citcoms_temp_dump * _ctd, float * dump, int len_dump,double (*f[])(double *,double *),  double * ctx);
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
int citcoms_eid(int n, int eid[], int nx, int ny, int nz);
void tracer_finder_init(tracer_finder * _tf,citcoms_dump * _cd,int p[4]);
int citcoms_check_tracer_element(citcoms_dump * _cd);
double solve_local(double * x,double * v2, double * v1, double * v0);
float solve_local_f(float * x,float * v2, float * v1, float * v0);
void set_projection_axis(double * n2, double * n1, double * n0, double * A, double * B, double *C, double *D);

int SphereIntegrateCitcomsDump(citcoms_dump * _cd, const char * _prefix);
int SphereIntegrateCitcomsDump2(CitcomsData * _cd, const char * _prefix);
CitcomsData * init_citcoms_data(const char * input);
int load_citcoms_step(CitcomsData * _cdata, int step);
int write_citcoms_step(CitcomsData * _cdata, int step);
int clean_citcoms_data(CitcomsData * _cdata);
int close_citcoms_data(CitcomsData * _cdata);
int init_sol_liq_prof(CitcomsData * _cdata);
int update_melting(CitcomsData * _cdata);

citcoms_sphere * init_citcoms_sphere2(CitcomsData * _cd, int noc);
citcoms_sphere * init_citcoms_sphere3(int npx, int nx, int noc);
int set_citcoms_sphere_coord(citcoms_sphere * _cs, int mcap);

int set_ring_scope(citcoms_sphere * _cs, const char * fname);
int vts_find_solidify_thickness(citcoms_sphere * _cs, CitcomsData * _cd);
int calculate_effective_depth(citcoms_sphere * _cs, CitcomsData * _cd, const char * tname);
#endif //SALECVTSREADER_CITCOMS_RELATED_H
