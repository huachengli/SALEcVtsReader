//
// Created by huacheng on 26/11/2023.
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
#define DeclareVtpSkip(name) void VtpTag##name##Process(const char* _values,VtpFile * _vfp) \
                            {_vfp->DataNodeType = SALEC_VTP_##name;}
#define SETVtpTag2P(name) {.TagName = #name,  .P = VtpTag##name##BeginFunc},\
                        {.TagName = "/"#name, .P = VtpTag##name##EndFunc},

DeclareVtpSkip(VTKFile)
DeclareVtpSkip(Points)
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
    // just process the NumberOfPoints property
    if(strcasecmp("NumberOfPoints",_vbuffer)==0)
    {
        r += Strok(_values+r," =\"",_vbuffer);
        _vfp->NoP = atoi(_vbuffer);
    } else
    {
        fprintf(stdout,"TagPiece undefined property:%s!",_vbuffer);
    }
}

DeclareVtpSkip(PolyData)
DeclareVtpSkip(PointData)

DeclareVtpFunc(VTKFile)
DeclareVtpFunc(Piece)
DeclareVtpFunc(DataArray)
DeclareVtpFunc(Points)
DeclareVtpFunc(PolyData)
DeclareVtpFunc(PointData)
VtpTagPair VtpTag2P[] = {
        SETVtpTag2P(VTKFile)
        SETVtpTag2P(Piece)
        SETVtpTag2P(DataArray)
        SETVtpTag2P(Points)
        SETVtpTag2P(PolyData)
        SETVtpTag2P(PointData)
        {.TagName = "unknown", .P = NULL}
};
