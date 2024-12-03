//
// Created by huachengli on 8/30/24.
//
#include "VtpReader.h"
#include "VtkWriter.h"
#include "InputParser.h"
#include "VtpTracer.h"
#include "Utility.h"
#include "lmath.h"

int check_depth(const double * pos, const double * ctx);

int main(int argc,char * argv[])
{
    char inp_file[200];
    switch(argc)
    {
        case 1:
            strcpy(inp_file,"planes.plot");
            break;
        case 2:
            strcpy(inp_file,argv[1]);
            break;
        default:
            fprintf(stdout,"error in command line\n");
    }
    InputFile * ifp = OpenInputFile(inp_file);

    // load position vts reference
    char salec_inp_name[4096];
    char salec_tracer_ref[4096];
    GetValueS(ifp,"SALEc.input",salec_inp_name,"SALEc.inp");
    GetValueS(ifp,"Tracer.ref",salec_tracer_ref,"ParaTest.proc%d.0.vts");

    SALEcData * ref = InitSALEcData(salec_inp_name,salec_tracer_ref);

    int step = 650;
    VtpTracerCollect * vtc = OpenSALEcTracerCollect(ifp,step);

    // ShowBriefVtpColleect(vtc,stdout);
    // VtpFile * tmp = SALEcVtpCollectMatFilter(vtc);
    // VtpFile * tmp = SALEcVtpCollectPosFilter(vtc);
    double ctx[3] = {1e6,1.74e6,0};
    VtpFile * tmp = SALEcVtpCollectPosFuncFilter(vtc,check_depth,ctx);
    VtpGetMelting(tmp);

    VtpGetConnect(tmp, vtc, ref);

    strcpy(tmp->name,"test.vtp");
    WriteVtpFile(tmp);
    CloseVtpFile(tmp);
    CloseVtpTracerCollect(vtc);
    CloseInputFile(ifp);
    CleanSALEcData(ref);
    return 0;
}

int check_depth(const double * pos, const double * ctx)
{
    double maxdepth = ctx[0];
    double Rm = ctx[1];
    if((1.0 - VecLen(pos,3)) < maxdepth/Rm)
    {
        return 1;
    }
    else
        return 0;
}