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

#define VTSDIM 3
#define MaxNameLen 50
#define VTSDATAFLOAT double
typedef struct
{
    char Name[MaxNameLen];
    unsigned int NoC; /* NumberOfComponent */
    VTSDATAFLOAT *** Data;
} VtsData;

typedef struct
{
    unsigned int PieceExtent[VTSDIM][2];
    unsigned int WholeExtent[VTSDIM][2];
    double *** Point;
    double *** Cell;

    unsigned int NoF; /*NumberOfField*/
    VtsData * Field;
} VtsInfo;

typedef struct
{
    unsigned int Tag;
    char Name[MaxNameLen];
} VtsStackFrame;

void VtsLoad(VtsInfo * _vfp,FILE * fp);

void Base64Encode(unsigned char * _str, unsigned _slen, unsigned char * _code);
int Base64Decode(unsigned char * _code, unsigned _clen, unsigned char * _str);
void RandomStr(unsigned char * _str,int _slen);
void zlibcompress(unsigned char* in, int nn, unsigned char** out, int *nn2);
void test_compress(Byte*,uLong,Byte*,uLong);

void IntToUnsignedChar(int * _iarray, int _inn, unsigned char * _carray);
void UnsignedCharToInt(unsigned char * _carray, int _cnn, int * _iarray);
#endif //SALECVTSREADER_VTSREADER_H
