//
// Created by huacheng on 10/18/21.
//

#include "VtsReader.h"

void TagVTKFileBeginFunc(const char* _values,VtsInfo * _vfp)
{
    VtsStackFrame * _vsf = _vfp->VtsStack + _vfp->StackPos;
    _vsf->Tag = SALEC_VTS_VTKFILE;
    strcpy(_vsf->Name,"VTKFileBegin");
    _vfp->StackPos ++;
}
void TagVTKFileEndFunc(const char* _values,VtsInfo * _vfp)
{
    VtsStackFrame * _vsf = _vfp->VtsStack + _vfp->StackPos - 1;
    if(_vsf->Tag!=SALEC_VTS_VTKFILE)
    {
        fprintf(stdout,"VTKFile Tag does not match!\n");
        exit(0);
    } else
    {
        _vfp->StackPos --;
    }
}

void TagStructuredGridBeginFunc(const char* _values,VtsInfo * _vfp)
{
    VtsStackFrame * _vsf = _vfp->VtsStack + _vfp->StackPos;
    _vsf->Tag = SALEC_VTS_STRUCTUREDGRID;
    strcpy(_vsf->Name,"StructuredGrid");
    _vfp->StackPos ++;

    char _vbuffer[50];
    int r = Strok(_values," =\"",_vbuffer);
    if(strcasecmp("WholeExtent",_vbuffer)==0)
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

        _vfp->WholeExtent[0][0] = ex0;
        _vfp->WholeExtent[0][1] = ex1;

        _vfp->WholeExtent[1][0] = ey0;
        _vfp->WholeExtent[1][1] = ey1;

        _vfp->WholeExtent[2][0] = ez0;
        _vfp->WholeExtent[2][1] = ez1;
    } else
    {
        fprintf(stdout,"StructuredGrid undefined property:%s!",_vbuffer);
    }
}
void TagStructuredGridEndFunc(const char* _values,VtsInfo * _vfp)
{
    VtsStackFrame * _vsf = _vfp->VtsStack + _vfp->StackPos - 1;
    if(_vsf->Tag!=SALEC_VTS_STRUCTUREDGRID)
    {
        fprintf(stdout,"StructuredGrid Tag does not match!\n");
        exit(0);
    } else
    {
        _vfp->StackPos --;
    }
}
void TagPieceBeginFunc(const char* _values,VtsInfo * _vfp)
{
    VtsStackFrame * _vsf = _vfp->VtsStack + _vfp->StackPos;
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

        _vfp->PieceExtent[0][0] = ex0;
        _vfp->PieceExtent[0][1] = ex1;

        _vfp->PieceExtent[1][0] = ey0;
        _vfp->PieceExtent[1][1] = ey1;

        _vfp->PieceExtent[2][0] = ez0;
        _vfp->PieceExtent[2][1] = ez1;
    } else
    {
        fprintf(stdout,"TagPiece undefined property:%s!",_vbuffer);
    }

    _vsf->Tag = SALEC_VTS_PIECE;
    strcpy(_vsf->Name,"Piece");
    _vfp->StackPos ++;
}

void TagPieceEndFunc(const char* _values,VtsInfo * _vfp)
{
    VtsStackFrame * _vsf = _vfp->VtsStack + _vfp->StackPos - 1;
    if(_vsf->Tag!=SALEC_VTS_PIECE)
    {
        fprintf(stdout,"StructuredGrid Tag does not match!\n");
        exit(0);
    } else
    {
        _vfp->StackPos --;
    }
}

void TagPointDataBeginFunc(const char* _values,VtsInfo * _vfp)
{
    VtsStackFrame * _vsf = _vfp->VtsStack + _vfp->StackPos;
    _vsf->Tag = SALEC_VTS_POINTDATA;
    strcpy(_vsf->Name,"PointData");
    _vfp->StackPos ++;
    _vfp->DataNodeType = SALEC_VTS_POINTDATA;
}

void TagPointDataEndFunc(const char* _values,VtsInfo * _vfp)
{
    VtsStackFrame * _vsf = _vfp->VtsStack + _vfp->StackPos - 1;
    if(_vsf->Tag!=SALEC_VTS_POINTDATA)
    {
        fprintf(stdout,"StructuredGrid Tag does not match!\n");
        exit(0);
    } else
    {
        _vfp->StackPos --;
    }
    _vfp->DataNodeType = SALEC_VTS_NONE;
}

void TagDataArrayBeginFunc(const char* _values,VtsInfo * _vfp)
{
    VtsStackFrame * _vsf = _vfp->VtsStack + _vfp->StackPos;
    _vsf->Tag = SALEC_VTS_DATAARRAY;
    strcpy(_vsf->Name,"DataArray");
    _vfp->StackPos ++;
    VtsData * _vdp = NULL;

    switch (_vfp->DataNodeType) {
        case SALEC_VTS_POINTS:
            _vdp = _vfp->PointField + _vfp->PointNoF++;
            break;
        case SALEC_VTS_CELLDATA:
            _vdp = _vfp->CellField + _vfp->CellNoF++;
            break;
        case SALEC_VTS_POINTDATA:
            _vdp = _vfp->PointField + _vfp->PointNoF++;
            break;
        case SALEC_VTS_NONE:
            _vdp = _vfp->PointField + _vfp->PointNoF++;
            break;
        default:
            _vdp = NULL;
            fprintf(stdout,"unknown DataNodeTye!\n");
            exit(0);
            break;
    }

    _vfp->ActiveVtsData = _vdp;
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

void TagDataArrayEndFunc(const char* _values,VtsInfo * _vfp)
{
    VtsStackFrame * _vsf = _vfp->VtsStack + _vfp->StackPos - 1;
    if(_vsf->Tag!=SALEC_VTS_DATAARRAY)
    {
        fprintf(stdout,"StructuredGrid Tag does not match!\n");
        exit(0);
    } else
    {
        _vfp->StackPos --;
    }
}

void TagCellDataBegin(const char* _values,VtsInfo * _vfp)
{
    VtsStackFrame * _vsf = _vfp->VtsStack + _vfp->StackPos;
    _vsf->Tag = SALEC_VTS_CELLDATA;
    strcpy(_vsf->Name,"CellData");
    _vfp->StackPos ++;
    _vfp->DataNodeType = SALEC_VTS_CELLDATA;
}

void TagCellDataEndFunc(const char* _values,VtsInfo * _vfp)
{
    VtsStackFrame * _vsf = _vfp->VtsStack + _vfp->StackPos - 1;
    if(_vsf->Tag!=SALEC_VTS_CELLDATA)
    {
        fprintf(stdout,"CellData Tag does not match!\n");
        exit(0);
    } else
    {
        _vfp->StackPos --;
    }
    _vfp->DataNodeType = SALEC_VTS_NONE;
}

void TagPointsBegin(const char* _values,VtsInfo * _vfp)
{
    VtsStackFrame * _vsf = _vfp->VtsStack + _vfp->StackPos;
    _vsf->Tag = SALEC_VTS_POINTS;
    strcpy(_vsf->Name,"Points");
    _vfp->StackPos ++;
    _vfp->DataNodeType = SALEC_VTS_POINTS;
}

void TagPointsEndFunc(const char* _values,VtsInfo * _vfp)
{
    VtsStackFrame * _vsf = _vfp->VtsStack + _vfp->StackPos - 1;
    if(_vsf->Tag!=SALEC_VTS_POINTS)
    {
        fprintf(stdout,"Points Tag does not match!\n");
        exit(0);
    } else
    {
        _vfp->StackPos --;
    }
    _vfp->DataNodeType = SALEC_VTS_NONE;
}


void ReadDataArrayType(const char * _vbuffer, unsigned int * _start,VtsData * _vdp)
{
    int t = *_start;
    *_start += Strok(_vbuffer+*_start," \"<>=",_vdp->Type);
}

void ReadDataArrayName(const char * _vbuffer, unsigned int * _start,VtsData * _vdp)
{
    *_start += Strok(_vbuffer+*_start," \"<>=",_vdp->Name);
}

void ReadDataArrayNoC(const char * _vbuffer, unsigned int * _start,VtsData * _vdp)
{
    char tmp[200];
    *_start += Strok(_vbuffer+*_start," \"<>=",tmp);
    _vdp->NoC = atoi(tmp);
}

void ReadDataArrayFormat(const char * _vbuffer, unsigned int * _start,VtsData * _vdp)
{
    *_start += Strok(_vbuffer+*_start," \"<>=",_vdp->Format);
}

/*
 * property label in DataArray
 */

char DataArrayProperty[][100] = {
        "Type",
        "Name",
        "NumberOfComponents",
        "format"
};


void (*ReadDataArrayProperty[])(const char * _vbuffer, unsigned int * _start,VtsData * _vdp) = {
        ReadDataArrayType,
        ReadDataArrayName,
        ReadDataArrayNoC,
        ReadDataArrayFormat
};

char TagName[][100] = {
        "VTKFile",
        "/VTKFile",
        "StructuredGrid",
        "/StructuredGrid",
        "Piece",
        "/Piece",
        "PointData",
        "/PointData",
        "DataArray",
        "/DataArray",
        "Points",
        "/Points",
        "CellData",
        "/CellData"
};

void (*TagNameP[])(const char*,VtsInfo *) = {
        TagVTKFileBeginFunc,
        TagVTKFileEndFunc,
        TagStructuredGridBeginFunc,
        TagStructuredGridEndFunc,
        TagPieceBeginFunc,
        TagPieceEndFunc,
        TagPointDataBeginFunc,
        TagPointDataEndFunc,
        TagDataArrayBeginFunc,
        TagDataArrayEndFunc,
        TagPointsBegin,
        TagPointsEndFunc,
        TagCellDataBegin,
        TagCellDataEndFunc
};