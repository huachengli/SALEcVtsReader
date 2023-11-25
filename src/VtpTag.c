//
// Created by huach on 26/11/2023.
//

#include "VtpReader.h"

#define VtpBeginFunc(name) void VtpTag##name##BeginFunc(const char* _values,VtpFile * _vfp) \
        { \
            VtsStackFrame * _vsf = _vfp->VtpStack + _vfp->StackPos; \
            _vsf->Tag = SALEC_VTP_##name; \
            strcpy(_vsf->Name,__func__); \
            _vfp->StackPos ++;                                                              \
            VtpTag##name##Process(_values,_vfp);\
        }

#define VtpEndFunc(name) void VtpTag##name##EndFunc(const char* _values,VtpFile * _vfp) \
        {\
            VtsStackFrame * _vsf = _vfp->VtpStack + _vfp->StackPos - 1; \
            if(_vsf->Tag!= SALEC_VTP_##name) \
            { \
                fprintf(stdout,#name" Tag does not match!\n"); \
                exit(0); \
            } else \
            { \
                _vfp->StackPos --; \
            } \
        }
#define DeclareVtpFunc(name) VtpBeginFunc(name) VtpEndFunc(name)
#define DeclareVtpSkip(name) void VtpTag##name##Process(const char* _values,VtpFile * _vfp){}

/*
 *  generate example DeclareVtpFunc(VTKFile)
void VtpTagVTKFileBeginFunc(const char* _values,VtpFile * _vfp)
{
    VtsStackFrame * _vsf = _vfp->VtpStack + _vfp->StackPos;
    _vsf->Tag = SALEC_VTP_VTKFILE;
    strcpy(_vsf->Name,__func__);
    _vfp->StackPos ++;
}

void VtpTagVTKFileEndFunc(const char* _values,VtpFile * _vfp)
{
    VtsStackFrame * _vsf = _vfp->VtpStack + _vfp->StackPos - 1;
    if(_vsf->Tag!= SALEC_VTP_VTKFILE)
    {
        fprintf(stdout,"VTKFile Tag does not match!\n");
        exit(0);
    } else
    {
        _vfp->StackPos --;
    }
}
 */

DeclareVtpSkip(VTKFile)
void VtpTagDataArrayProcess(const char* _values,VtpFile * _vfp)
{
    VtsData * _vdp = NULL;
    switch (_vfp->DataNodeType) {
        case SALEC_VTP_Points:
            _vdp = _vfp->PointField + _vfp->PointNoF++;
            break;
        case SALEC_VTP_PointData:
            _vdp = _vfp->PointField + _vfp->PointNoF++;
            break;
        case SALEC_VTP_NONE:
            _vdp = _vfp->PointField + _vfp->PointNoF++;
            break;
        default:
            _vdp = NULL;
            fprintf(stdout,"unknown DataNodeTye!\n");
            exit(0);
            break;
    }

    _vfp->ActiveData = _vdp;
    _vdp->NoC = 1; // set the default values for compoent;

    unsigned int r = 0;
    unsigned int dr = 0;
    char OptName[200];
    while((dr=Strok(_values+r," \"<>=",OptName))>0)
    {
        r+= dr;
        int k =0;
        for(;k<SALEC_VTS_DATAARRAY_OPTS;k++)
        {
            if(0==strcasecmp(DataArrayProperty[k],OptName)) break;
        }
        if(SALEC_VTS_DATAARRAY_OPTS==k)
        {
            fprintf(stdout,"DataArray:undefined properity:%s\n",OptName);
            exit(0);
        }
        ReadDataArrayProperty[k](_values,&r,_vdp);
    }
}

void VtpTagPieceProcess(const char* _values,VtpFile * _vfp)
{
    VtsStackFrame * _vsf = _vfp->VtpStack + _vfp->StackPos;
    char _vbuffer[50];
    int r = Strok(_values," =\"",_vbuffer);
    if(strcasecmp("Extent",_vbuffer)==0)
    {
        r += Strok(_values+r," =\"",_vbuffer);
        int ex0 = atoi(_vbuffer);
        r += Strok(_values+r," =\"",_vbuffer);
        int ex1 = atoi(_vbuffer);
        r += Strok(_values+r," =\"",_vbuffer);
        int ey0 = atoi(_vbuffer);
        r += Strok(_values+r," =\"",_vbuffer);
        int ey1 = atoi(_vbuffer);
        r += Strok(_values+r," =\"",_vbuffer);
        int ez0 = atoi(_vbuffer);
        r += Strok(_values+r," =\"",_vbuffer);
        int ez1 = atoi(_vbuffer);
    } else
    {
        fprintf(stdout,"TagPiece undefined property:%s!",_vbuffer);
    }
}

DeclareVtpFunc(VTKFile)
DeclareVtpFunc(Piece)


VtpTagPair VtpTag2P[] = {
        {.TagName = "VTKFile",  .P = VtpTagVTKFileBeginFunc},
        {.TagName = "/VTKFile", .P = VtpTagVTKFileEndFunc},
        {.TagName = "XXXXXXXX", .P = NULL}
};
