//
// Created by lhc on 27/7/2024.
//

#ifndef SALECVTSREADER_SLICEFILTER_H
#define SALECVTSREADER_SLICEFILTER_H
#include "Utility.h"

typedef struct SliceFilerImpl
{
    // section for plane data, do not change it!
    VTSDATAFLOAT n[VTSDIM];
    VTSDATAFLOAT d;
    char Name[MaxNameLen];   // unused;
    char Vacuum[MaxNameLen]; // unused;

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
    // end of plane struct data
    VTSDATAFLOAT * data_coord;
    VTSDATAFLOAT * data_vars;
} SliceFilter;

typedef SliceFilter slice_filter;

int GetSliceDataC(SALEcData * _sdata, SliceFilter * _out, unsigned long Id);
void CleanSlice(SliceFilter * _out);
void SetSliceMask(SALEcData * _sdata, SliceFilter * _out,int (*_search)(VtsInfo *, VTSDATAFLOAT *, int *,VTSDATAFLOAT*));


#endif //SALECVTSREADER_SLICEFILTER_H