//
// Created by huach on 26/11/2023.
//

#ifndef SALECVTSREADER_VTUREADER_H
#define SALECVTSREADER_VTUREADER_H

#include "VtsReader.h"
#define VtpMaxNameLen 50
#define VtpMaxStackDepth 50
#define VTPDATAFLOAT float

#define SALEC_VTP_TAG_TYPES 14
#define SALEC_VTP_HEADER 1
#define SALEC_VTP_VTKFile 5
#define SALEC_VTP_Piece 13
#define SALEC_VTP_PointData 17
#define SALEC_VTP_DataArray 21
#define SALEC_VTP_DATAARRAY_OPTS 4
#define SALEC_VTP_Points 29
#define SALEC_VTP_PolyData 31
#define SALEC_VTP_NONE 0
#define VTPDIM 3
typedef VtsData VtpData;
typedef struct vtp_file
{
    unsigned int NoP; // number of points
    VTPDATAFLOAT * Point;

    unsigned int PointNoF;
    VtpData PointField[MaxStackDepth];

    // data for process vtk file
    VtsStackFrame VtpStack[MaxStackDepth];
    unsigned int StackPos;
    unsigned int DataNodeType;
    VtpData * ActiveData;
    char name[4096];
} VtpFile;

VtpFile * OpenVtpFile(const char vtp_name[]);
int CloseVtpFile(VtpFile * fp);
int VtpFrameHeadLoad(VtpFile * _vfp,FILE *fp);
int VtpFrameLoad(VtpFile * _vsf,FILE *fp);
int ShowVtpFileInfo(VtpFile * vfp);
int VtpCoordinateReshape(VtpFile * _vsf);

typedef struct {
    char TagName[MaxStrLen];
    void (*P)(const char*,VtpFile *);
} VtpTagPair;
extern VtpTagPair VtpTag2P[];
#endif //SALECVTSREADER_VTUREADER_H
