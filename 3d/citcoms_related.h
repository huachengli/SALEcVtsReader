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


typedef struct CitcomsDumpImpl
{
    int nproc;
    char temp_prefix[4096];
    char tracer_prefix[4096];
    double TransformR;
    citcoms_temp_dump * temp;
    citcoms_tracer_dump * tracer;
} citcoms_dump;

int load_citcoms_temp_dump(citcoms_temp_dump * _ctd, const char * _fname);
int load_citcoms_tracer_dump(citcoms_tracer_dump * _ctd, const char * _fname);
int load_citcoms_dump(citcoms_dump * _cd, InputFile * ifp);

int clean_citcoms_temp_dump(citcoms_temp_dump * _ctd);
int clean_citcoms_tracer_dump(citcoms_tracer_dump * _ctd);
int clean_citcoms_dump(citcoms_dump * x);
citcoms_dump * InitCitcomsDump(InputFile * ifp);
int CloseCitcomsDump(citcoms_dump * x);
SALEcData * CrInitSALEcData(InputFile * ifp);
void CrCloseSALEcData(SALEcData * _sdata);

int UpdateCitcomsTempDump(citcoms_dump * _cd, SALEcData * _sdata);
int UpdateCitcomsTracerDump(citcoms_dump * _cd, SALEcData * _sdtat);
int UpdateCitcomsDump(citcoms_dump * _cdp, SALEcData * _sdata);
int SALEcGetCData(SALEcData * _sdata, int fId, VTSDATAFLOAT * _pos, VTSDATAFLOAT * _data);

void citcoms_tracer_dump_vtp(citcoms_tracer_dump * _ctd, const char * name);
#endif //SALECVTSREADER_CITCOMS_RELATED_H
