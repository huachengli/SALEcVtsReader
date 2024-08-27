//
// Created by huachengli on 8/26/24.
//

#include "citcoms_related.h"

int main(int argc,char * argv[])
{
    char inp_file[4096] = "set_dump.inp";
    switch(argc)
    {
        case 1:
            strcpy(inp_file,"set_dump.inp");
            break;
        case 2:
            strcpy(inp_file,argv[1]);
            break;
        default:
            fprintf(stdout,"error in command line\n");
    }

    InputFile * ifp = OpenInputFile(inp_file);
    citcoms_dump * cdp = InitCitcomsDump(ifp);
    SALEcData * sdp = CrInitSALEcData(ifp);
    UpdateCitcomsDump(cdp,sdp);
    CrCloseSALEcData(sdp);
    CloseCitcomsDump(cdp);
    CloseInputFile(ifp);
    return 0;
}