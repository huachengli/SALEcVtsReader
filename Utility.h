//
// Created by huacheng on 10/19/21.
//

#ifndef SALECVTSREADER_UTILITY_H
#define SALECVTSREADER_UTILITY_H
#include "VtsReader.h"

typedef struct
{
    unsigned long Npx[VTSDIM];
    unsigned long Npgx[VTSDIM];
    unsigned long Noffset;
    unsigned long VtsBlockNum;
    char VtsPrefix[200];
    VTSDATAFLOAT * GCLV[VTSDIM]; // global coordinate line (vertex)
    VTSDATAFLOAT * GCLC[VTSDIM]; // global coordinate line (cell)
    unsigned long nGCLV[VTSDIM];
    unsigned long nGCLC[VTSDIM];

    VTSDATAFLOAT * BCLV[VTSDIM]; // block coordinate line (vertex)
    VTSDATAFLOAT * BCLC[VTSDIM]; // block coordinate line (cell)
    unsigned long nBCLV[VTSDIM];
    unsigned long nBCLC[VTSDIM];


    VtsInfo * VSF;
} SALEcData;

typedef struct
{
    VTSDATAFLOAT n[VTSDIM];
    VTSDATAFLOAT d;

    unsigned long shape[VTSDIM-1];
    unsigned long NoC;
    unsigned long Id;
    char Name[MaxNameLen];
    VTSDATAFLOAT nX[VTSDIM-1][VTSDIM];
    VTSDATAFLOAT X0[VTSDIM];
    VTSDATAFLOAT * CL[VTSDIM-1];
    VTSDATAFLOAT ** data;
    int ** mask;
    unsigned long nCL[VTSDIM-1];

    unsigned long * scores;
} Plane;

void LoadInpInfo(SALEcData *,const char*);
void LoadVtsData(SALEcData *,const char *);
void CleanSALEcData(SALEcData *);
void WriteGCL(SALEcData * _sdata,unsigned int _c,FILE *fp);

int cbsearch(VTSDATAFLOAT * _cl, VTSDATAFLOAT _x,unsigned long _len);
void IndexSearchC(SALEcData *,VTSDATAFLOAT const *,int*);
void IndexSearchV(SALEcData *,VTSDATAFLOAT const *,int*);
int GetBlockIndexC(SALEcData *, const int *);
int GetBlockIndexV(SALEcData *, const int *);
int BlockSearch(SALEcData * _sdata,VTSDATAFLOAT * _x);
#define vdot3(a,b) (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])
#define vzero(a) a[0] = a[1] = a[2] = 0.
void SetPlaneMeshV(SALEcData * _sdata, Plane * _out);
void SetPlaneMask(SALEcData * _sdata, Plane * _out);
#endif //SALECVTSREADER_UTILITY_H
