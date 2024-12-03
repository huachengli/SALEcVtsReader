//
// Created by huachengli on 10/22/24.
//

#include "citcoms_related.h"

int main(int argc,char * argv[])
{
    char inp_file[4096] = "post.inp";
    switch(argc)
    {
        case 1:
            break;
        case 2:
            strcpy(inp_file,argv[1]);
            break;
        default:
            fprintf(stdout,"error in command line\n");
    }

    CitcomsData * _cdata = init_citcoms_data(inp_file);

    for(int k=_cdata->step0;k<_cdata->step1;k+=_cdata->step_inc)
    {
        load_citcoms_step(_cdata,k);
        char tmp_prefix[4096];
        snprintf(tmp_prefix, 4096, "%s.%04d", _cdata->OutPrefix, k);
        SphereIntegrateCitcomsDump2(_cdata,tmp_prefix);
        clean_citcoms_data(_cdata);
        fprintf(stdout, "write sphere data to %s.vtm\n", tmp_prefix);
    }
    close_citcoms_data(_cdata);

    return 0;
}