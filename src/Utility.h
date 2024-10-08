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
    char Name[MaxNameLen];
    char Vacuum[MaxNameLen];

    unsigned long shape[VTSDIM-1];
    unsigned long NoC;
    unsigned long Id;
    unsigned long VacuumId;

    VTSDATAFLOAT nX[VTSDIM-1][VTSDIM];
    VTSDATAFLOAT X0[VTSDIM];
    VTSDATAFLOAT * CL[VTSDIM-1];
    VTSDATAFLOAT *** data;
    int ** mask;
    unsigned long nCL[VTSDIM-1];

    unsigned long * scores;
    int *** offset;
    VTSDATAFLOAT *** weight;
    VTSDATAFLOAT *** coord;
} Plane;

typedef struct
{
    int CacheState;
    char Vacuum[MaxNameLen];
    unsigned long VacuumId;
    VTSDATAFLOAT n[VTSDIM];
    VTSDATAFLOAT d;
    VTSDATAFLOAT tol;

    unsigned long nCL[VTSDIM-1];
    int ** mask;
    int ** kz;
    VTSDATAFLOAT ** Lambda;
} ProfileCache;


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

int OffsetSerchC(VtsInfo * _vsf, VTSDATAFLOAT * _x, int * _indices,VTSDATAFLOAT *);
int OffsetSerchV(VtsInfo * _vsf, VTSDATAFLOAT * _x, int * _indices,VTSDATAFLOAT *);
void SetPlaneMask(SALEcData * _sdata, Plane * _out,int (*_search)(VtsInfo *, VTSDATAFLOAT *, int *,VTSDATAFLOAT*));
#define SetPlaneMaskV(_sdata,_out) SetPlaneMask(_sdata,_out,OffsetSerchV)
#define SetPlaneMaskC(_sdata,_out) SetPlaneMask(_sdata,_out,OffsetSerchC)

void SetMeshCLV(SALEcData * _sdata,Plane * _out,unsigned int px,unsigned int py);
void SetMeshCLC(SALEcData * _sdata,Plane * _out,unsigned int px,unsigned int py);
void SetPlaneMesh(SALEcData * _sdata, Plane * _out, void (*_set)(SALEcData * ,Plane * ,unsigned int ,unsigned int));
#define SetPlaneMeshV(_sdata,_out) SetPlaneMesh(_sdata,_out,SetMeshCLV)
#define SetPlaneMeshC(_sdata,_out) SetPlaneMesh(_sdata,_out,SetMeshCLC)

void GetPlaneDataC(SALEcData * _sdata, Plane * _out);
void WritePlaneData(Plane * _out,int compoent,const char * fname);
void WritePlaneDataAll(Plane * _out,const char * fname);
void WritePlaneCoord(Plane * _out, const char * fname);
void CleanPlane(Plane * _out);
void CleanCache(ProfileCache * _cache,int _n);

VTSDATAFLOAT* VtmGetCellData(SALEcData * _sdata, unsigned long k, unsigned long _i,unsigned long _j, unsigned long _k);

#define v_normalize(x) do {\
    VTSDATAFLOAT tmp = sqrt((x)[0]*(x)[0] + (x)[1]*(x)[1] + (x)[2]*(x)[2]); \
    (x)[0] /= tmp;(x)[1] /= tmp;(x)[2] /= tmp;                 \
    }while(0);
#define v_liner_op(x,pa,a,pb,b) do \
    {                              \
        for(int _vk=0;_vk<3;_vk++) (x)[_vk] = (pa)*(a)[_vk] + (pb)*(b)[_vk];\
    } while(0);
#define v_copy(x,a) do \
    {\
        for(int _vk=0;_vk<3;_vk++) (x)[_vk] = (a)[_vk];\
    }while(0);
#define v_add_liner(x,pa,a,pb,b) do \
    {                              \
        for(int _vk=0;_vk<3;_vk++) (x)[_vk] += (pa)*(a)[_vk] + (pb)*(b)[_vk];\
    } while(0);
#define vdot3(a,b) ((a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2])
#define vzero(a) (a)[0] = (a)[1] = (a)[2] = 0.
#define lerp(_l,_r,_t) ((_l)*(1.-_t) + (_r)*(_t))
#define vlerp(_x,_n,_l,_r,_t) do{\
    for(int _iv=0;_iv<(_n);_iv++)\
    {                        \
        (_x)[_iv] = lerp((_l)[_iv],(_r)[_iv],(_t));\
    }}while(0);
#define _lId3(x,y,z,nx,ny,nz) ((x)+(nx)*((y)+(ny)*(z)))
#define _lId2(x,y,nx,ny) ((x)+(nx)*(y))

#define print_vec(x,nx,t) do{\
    fprintf(stdout,"\n###START###\n");   \
    for(int _vk=0;_vk<(nx);_vk++) \
    {\
        fprintf(stdout,"%10.6f,",(x)[_vk]);\
        if(_vk%(t)==(t-1)) fprintf(stdout,"\n");\
    }                                            \
    fprintf(stdout,"\n###END###\n");     \
    }while(0);

char * head__strcasestr(const char * src, const char *tgt);
unsigned long find_cellfield(const char _src[], VtsInfo * _vsf);

// function using cache to search profile
void GetProfileLim(SALEcData * _sdata, Plane * _out, VTSDATAFLOAT _tol);
void GetProfileVOFLim(SALEcData * _sdata, Plane * _out, VTSDATAFLOAT _tol);

#define GetProfile(_sdata,_out) GetProfileLim(_sdata,_out,0.95)
#define GetProfileVOF(_sdata,_out) GetProfileVOFLim(_sdata,_out,0.95)

void GetProfileVOFReadCache(SALEcData * _sdata, Plane * _out,ProfileCache * _cache);
void GetProfileVOFWriteCache(SALEcData * _sdata, Plane * _out, ProfileCache * _cache, VTSDATAFLOAT _tol);
void GetProfileVOFWithCache__(SALEcData * _sdata, Plane * _out, ProfileCache * _cache, VTSDATAFLOAT _tol);
#define GetProfileVOFWithCache(_sdata,_out,_cache) GetProfileVOFWithCache__(_sdata,_out,_cache,0.95)

void GetProfileWithCache__(SALEcData * _sdata, Plane * _out, ProfileCache * _cache, VTSDATAFLOAT _tol);
#define GetProfileWithCache(_sdata,_out,_cache) GetProfileWithCache__(_sdata,_out,_cache,0.05)

void GetRemnantLim(SALEcData * _sdata, Plane * _out, VTSDATAFLOAT _tol);
#define GetRemant(_sdata,_out) GetRemnantLim(_sdata,_out,0.05)
#define SumRemant(_sdata,_out) GetRemnantLim(_sdata,_out,-1.0)

#endif //SALECVTSREADER_UTILITY_H
