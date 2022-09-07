//
// Created by huacheng on 3/23/22.
//

#include "InputParser.h"
#include "TracerReader.h"

int main(int argc,char * argv[])
{
    char get_plane[200];
    switch (argc) {
        case 1:
            strcpy(get_plane,"planes.plot");
            break;
        case 2:
            strcpy(get_plane,argv[1]);
            break;
        default:
            fprintf(stdout,"error in command line\n");
    }


    char SALEcInpName[200];
    TracerInfo pTracerInfo;

    InputFile * ifp = OpenInputFile(get_plane);
    GetValueS(ifp,"SALEc.input",SALEcInpName,"SALEc.inp");
    GetValueS(ifp,"tracer.path",pTracerInfo.path,"unknown");
    CloseInputFile(ifp);
    LoadTracerInfo(&pTracerInfo,SALEcInpName);


    Tracer ** pTracer = (Tracer **) malloc(sizeof(Tracer*)*pTracerInfo.nproc);
    Tracer ** cTracer = (Tracer **) malloc(sizeof(Tracer*)*pTracerInfo.nproc);
    for(int k=0;k<pTracerInfo.nproc;++k)
    {
        char tname[200];
        TracerName(&pTracerInfo,k,0,tname);
        pTracer[k] = OpenTracerFile(tname);
        TracerName(&pTracerInfo,k,8,tname);
        cTracer[k] = OpenTracerFile(tname);

        fprintf(stdout,"[%d]:%f\n",k, CompareTracerFile(pTracer[k],cTracer[k]));
        if(k>0)
        fprintf(stdout,"[%d+]:%f\n",k, CompareTracerFile(pTracer[k-1],cTracer[k]));

    }

    // check the instance






    for(int k=0;k<pTracerInfo.nproc;++k)
    {
        CloseTracerFile(pTracer[k]);
        CloseTracerFile(cTracer[k]);

    }
    free(pTracer);
    free(cTracer);
    return 0;
}





