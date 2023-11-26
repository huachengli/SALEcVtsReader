//
// Created by huacheng on 10/14/21.
//

#ifndef SALECVTSREADER_VTSREADER_H
#define SALECVTSREADER_VTSREADER_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>
#include "InputParser.h"


#define VTSDIM 3
#define MaxNameLen 50
#define MaxStackDepth 50
#define VTSDATAFLOAT float

#define SALEC_VTS_TAG_TYPES 14
#define SALEC_VTS_HEADER 1
#define SALEC_VTS_VTKFILE 5
#define SALEC_VTS_STRUCTUREDGRID 9
#define SALEC_VTS_PIECE 13
#define SALEC_VTS_POINTDATA 17
#define SALEC_VTS_DATAARRAY 21
#define SALEC_VTS_DATAARRAY_OPTS 4
#define SALEC_VTS_CELLDATA 25
#define SALEC_VTS_POINTS 29
#define SALEC_VTS_NONE 0

typedef struct
{
    char Name[MaxNameLen];
    char Type[MaxNameLen];
    char Format[MaxNameLen];
    unsigned int NoC; /* NumberOfComponent */
    VTSDATAFLOAT * Data;
    unsigned long DataLen;
} VtsData;

typedef struct
{
    unsigned int Tag;
    char Name[MaxNameLen];
} VtsStackFrame;

typedef struct
{
    unsigned long PieceExtent[VTSDIM][2];
    unsigned long WholeExtent[VTSDIM][2];
    unsigned long Nxp[VTSDIM];
    VTSDATAFLOAT **** Point;
    VTSDATAFLOAT **** Cell;
    VTSDATAFLOAT * CLV[VTSDIM];
    VTSDATAFLOAT * CLC[VTSDIM];

    unsigned int CellNoF; /*NumberOfField*/
    VtsData * CellField;

    unsigned int PointNoF;
    VtsData * PointField;

    VtsStackFrame * VtsStack;
    unsigned int StackPos;
    unsigned int DataNodeType;
    VtsData * ActiveVtsData;
} VtsInfo;


/*
 * Base64 encoding/uncoding
 */
void Base64Encode(unsigned char * _str, unsigned _slen, unsigned char * _code);
int Base64Decode(unsigned char * _code, unsigned _clen, unsigned char * _str);
void RandomStr(unsigned char * _str,int _slen);


void ReadVtsBinaryF32(float ** _data,unsigned long * _dlen,FILE * fp);
void VtsLoad(VtsInfo * _vfp,FILE * fp);
int VtsFrameHeadLoad(VtsInfo * _vsp,FILE *fp);
int VtsFrameLoad(VtsInfo * _vsf,FILE *fp);
int VtsCoordinateReshape(VtsInfo * _vsf);
void VtsInfoClean(VtsInfo * _vsf);
VTSDATAFLOAT * VtsGetPoint(VtsInfo * _vsf,unsigned long _i, unsigned long _j, unsigned long _k);
VTSDATAFLOAT * VtsGetCellData(VtsInfo * _vsf,unsigned long k,unsigned long _i, unsigned long _j, unsigned long _k);
VTSDATAFLOAT * VtsGetPointData(VtsInfo * _vsf,unsigned long k,unsigned long _i, unsigned long _j, unsigned long _k);
void VtsSetCoordLine(VtsInfo * _vsf);

/*
 * Some Tag processing function in VtsTag
 */
extern void (*ReadDataArrayProperty[])(const char * _vbuffer, unsigned int * _start,VtsData * _vdp);
extern char DataArrayProperty[][100];
extern char TagName[][100];
extern void (*TagNameP[])(const char*,VtsInfo *);
#endif //SALECVTSREADER_VTSREADER_H
