//
// Created by huacheng on 3/23/22.
//

#ifndef SALECVTSREADER_TRACERREADER_H
#define SALECVTSREADER_TRACERREADER_H

#include "InputParser.h"
#include "AVLTree.h"

typedef struct _TracerInfo
{
    int nproc;
    int maxstep;

    char type[200];
    char prefix[200];
    char path[200];
} TracerInfo;

typedef struct  _TracerData
{
    double position[3];
    double data[3];
    int Id;
    int eId;
} TracerData;

typedef struct _Tracer
{
    int num;
    TracerData * data;
    AVLNode * node; // memory allocate for avlnode
    AVLNode * root; // the root used to query
} Tracer;


void LoadTracerInfo(TracerInfo * _tinfo,const char *);
void TracerName(TracerInfo * _tinfo,int rank,int step,char * _tname);
int LoadTracerFile(const char * fname, Tracer * _tracer);
Tracer * OpenTracerFile(const char * fname);

inline static int TracerIdCmp(AVLNode * a,AVLNode * b)
{
    int Id_a = ((TracerData *)a->attach)->Id;
    int Id_b = ((TracerData *)b->attach)->Id;
    return Id_b - Id_a;
}

void CloseTracerFile(Tracer * _tracer);
double CompareTracerFile(Tracer * a,Tracer * b);
#endif //SALECVTSREADER_TRACERREADER_H
