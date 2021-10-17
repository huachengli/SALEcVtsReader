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

#define VTSDIM 3
#define MaxNameLen 50
#define MaxStackDepth 512
#define VTSDATAFLOAT double

#define SALEC_VTS_TAG_TYPES 14
#define SALEC_VTS_HEADER 1
#define SALEC_VTS_VTKFILE 5
#define SALEC_VTS_STRUCTUREDGRID 9

typedef struct
{
    char Name[MaxNameLen];
    unsigned int NoC; /* NumberOfComponent */
    VTSDATAFLOAT *** Data;
} VtsData;

typedef struct
{
    unsigned int Tag;
    char Name[MaxNameLen];
} VtsStackFrame;

typedef struct
{
    unsigned int PieceExtent[VTSDIM][2];
    unsigned int WholeExtent[VTSDIM][2];
    double *** Point;
    double *** Cell;

    unsigned int NoF; /*NumberOfField*/
    VtsData * Field;

    VtsStackFrame * VtsStack;
    unsigned int StackPos;
} VtsInfo;


/*
 * Base64 encoding/uncoding
 */
void Base64Encode(unsigned char * _str, unsigned _slen, unsigned char * _code);
int Base64Decode(unsigned char * _code, unsigned _clen, unsigned char * _str);
void RandomStr(unsigned char * _str,int _slen);

/*
 * just copy from citcoms, never used.
 */
void zlibcompress(unsigned char* in, int nn, unsigned char** out, int *nn2);
void IntToUnsignedChar(int * _iarray, int _inn, unsigned char * _carray);
void UnsignedCharToInt(unsigned char * _carray, int _cnn, int * _iarray);
void FloatToUnsignedChar(float * _farray, int _fnn, unsigned char * _carray);

void ReadVtsBinaryF32(float * _data,unsigned long * _dlen,FILE * fp);
void VtsLoad(VtsInfo * _vfp,FILE * fp);
int VtsFrameHeadLoad(VtsInfo * _vsp,FILE *fp);
int VtsFrameLoad(VtsStackFrame * _vsf,FILE *fp);

/*
 * basic string processing
 */
char *trim(char *str);
int Strok(const char _str[],const char _delim[], char value[]);
int InSubset(char _c,const char _set[]);
int ReadLineTrim(unsigned char _buffer[],FILE *fp);
#endif //SALECVTSREADER_VTSREADER_H
