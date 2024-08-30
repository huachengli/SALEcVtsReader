//
// Created by huachengli on 8/30/24.
//
#include "VtpReader.h"
#include "VtkWriter.h"
#include "InputParser.h"
#include "VtpTracer.h"
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
    int step = 950;
    VtpTracerCollect * vtc = OpenSALEcTracerCollect(ifp,step);

    ShowBriefVtpColleect(vtc,stdout);
    VtpFile * tmp = SALEcVtpCollectFilter(vtc);
    strcpy(tmp->name,"test.vtp");
    WriteVtpFile(tmp);
    CloseVtpFile(tmp);
    CloseVtpTracerCollect(vtc);
    CloseInputFile(ifp);
    return 0;
}