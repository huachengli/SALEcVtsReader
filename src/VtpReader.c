//
// Created by huach on 26/11/2023.
//

#include "VtpReader.h"

void VtpLoad(VtpFile * _vfp, FILE *fp)
{
    _vfp->VtpStack = (VtsStackFrame *) malloc(sizeof(VtsStackFrame) * MaxStackDepth);
    _vfp->StackPos = 0;
    _vfp->PointField = (VtsData *) malloc(sizeof(VtsData)*MaxStackDepth);
    _vfp->PointNoF = 0;
    VtpFrameHeadLoad(_vfp,fp);
    while(VtpFrameLoad(_vfp,fp))
    {
        VtsStackFrame * _vsf = _vfp->StackPos - 1 + _vfp->VtpStack;
        if(_vsf->Tag == SALEC_VTP_DataArray)
        {
            ReadVtsBinaryF32(&(_vfp->ActiveData->Data), &(_vfp->ActiveData->DataLen), fp);
        }
    }
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
    VtpLoad(tmp,fp);
    fclose(fp);
    return tmp;
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
        fprintf(stdout,"Wranning/the header of vts is not consistent with SALEc!\n");
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


