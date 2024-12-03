//
// Created by huach on 26/11/2023.
//

#include "VtpReader.h"
#include "VtkWriter.h"

void VtpLoad(VtpFile * _vfp, FILE *fp)
{
    _vfp->StackPos = 0;
    _vfp->PointNoF = 0;
    VtpFrameHeadLoad(_vfp,fp);
    while(VtpFrameLoad(_vfp,fp))
    {
        VtsStackFrame * _vsf = _vfp->StackPos - 1 + _vfp->VtpStack;
        if(_vsf->Tag == SALEC_VTP_DataArray)
        {
            if(strcasecmp(_vfp->ActiveData->Format,"binary")==0)
            {
                ReadVtsBinaryF32(&(_vfp->ActiveData->Data), &(_vfp->ActiveData->DataLen), fp);
            }
        }

        if(_vsf->Tag == SALEC_VTP_Piece && _vfp->NoP == 0)
        {
            break;
        }
    }
    VtpCoordinateReshape(_vfp);
}

VtpFile * OpenVtpFile(const char vtp_name[])
{
    FILE * fp = fopen(vtp_name,"r");
    if(fp==NULL)
    {
        fprintf(stdout,"cannot open %s\n",vtp_name);
        exit(0);
    }
    VtpFile * tmp = malloc(sizeof(VtpFile));
    snprintf(tmp->name,4096,"%s",vtp_name);
    VtpLoad(tmp,fp);
    fclose(fp);
    return tmp;
}

int CloseVtpFile(VtpFile * vfp)
{
    if(NULL == vfp) return 0;
    // clear PointField
    for(int k=0;k<vfp->PointNoF;++k)
    {
        free(vfp->PointField[k].Data);
    }
    // clear Points
    // Points has been cleared in PointsData
    return 1;
}

int ShowVtpFileInfo(VtpFile * vfp)
{
    fprintf(stdout,"%s\n",vfp->name);
    for(int k=0;k<vfp->PointNoF;++k)
    {
        fprintf(stdout,"PointFiles:%s,len=%ld\n",vfp->PointField[k].Name,vfp->PointField[k].DataLen);
    }
    return 0;
}


int VtpFrameHeadLoad(VtpFile * _vfp,FILE *fp)
{
    /*
     * Read the first line of vts file
     */

    unsigned char LineBuffer[1024];
    ReadLineTrim(LineBuffer,fp);
    VtsStackFrame * _vsf = _vfp->VtpStack + _vfp->StackPos;

    // check header
    unsigned char SALEcVtsHead[] = "<?xml version=\"1.0\"?>";
    if(0!= strcmp(SALEcVtsHead,LineBuffer))
    {
        fprintf(stdout,"Warning/the header of vts is not consistent with SALEc!\n");
        exit(0);
    }

    // set tag for header
    strcpy(_vsf->Name,LineBuffer);
    _vsf->Tag = SALEC_VTS_HEADER;
    _vfp->StackPos ++;
    return 1;
}

int VtpFrameLoad(VtpFile * _vsf,FILE *fp)
{
    unsigned char LineBuffer[1024];
    ReadLineTrim(LineBuffer,fp);

    char tKey[100];
    int r = Strok(LineBuffer+1," <>",tKey)+1;
    int k = 0;
    for(;VtpTag2P[k].P != NULL;++k)
    {
        if(0==strcmp(tKey,VtpTag2P[k].TagName)) break;
    }
    if(VtpTag2P[k].P == NULL)
    {
        fprintf(stdout,"Undefined TagName:%s\n",tKey);
        exit(0);
    } else
    {
        VtpTag2P[k].P(LineBuffer+r,_vsf);
    }
    return _vsf->StackPos - 1;
}

int VtpCoordinateReshape(VtpFile * _vsf)
{
    int k = 0;
    for(;k<_vsf->PointNoF;++k)
    {
        if(0== strcasecmp("coordinate",_vsf->PointField[k].Name)) break;
    }

    if(k==_vsf->PointNoF)
    {
        if(_vsf->NoP == 0) return 0;
        fprintf(stdout,"No coordinate information in the vtp!\n");
        return -1;
    }
    else
    {
        VTSDATAFLOAT * PointData = _vsf->PointField[k].Data;
        unsigned long PointDataLen = _vsf->PointField[k].DataLen;

        if(_vsf->NoP*VTPDIM!=PointDataLen)
        {
            fprintf(stdout,"number of coordinates is %ld, but %d is wanted\n",PointDataLen,_vsf->NoP*VTPDIM);
            exit(0);
        }
        // set Points ptr
        _vsf->Point = PointData;
        return (int) PointDataLen;
    }
}

int WriteVtpFile(VtpFile * _vsf)
{
    FILE * fp = fopen(_vsf->name,"w");
    vtp_file_header(fp,_vsf->NoP);
    vtk_point_data_header(fp);

    for(int k=0;k<_vsf->PointNoF;++k)
    {
        VtpData * _vdk = _vsf->PointField + k;
        if(strcasecmp(_vdk->Name,"coordinate") == 0)
            continue;
        vtk_dataarray_vec_f2(fp,_vdk->Name,_vdk->Format,_vdk->Type,_vdk->Data,_vdk->DataLen/_vdk->NoC,_vdk->NoC);
    }
    vtk_point_data_trailer(fp);
    vtk_point_header(fp);
    vtk_dataarray_vec_f(fp,"coordinate","binary",_vsf->Point,_vsf->NoP,3);
    vtk_point_trailer(fp);
    vtp_file_trailer(fp);
    fclose(fp);

    return _vsf->NoP;
}


