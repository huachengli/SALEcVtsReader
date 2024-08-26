//
// Created by li on 8/20/24.
//
// convert ejecta txt files to vtps


#include "ejecta_analysis.h"
#include "InputParser.h"
#include <stdio.h>
#include <assert.h>

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

    ejecta_collect EC;
    // ejecta_collect_test_init(&EC);
    ejecta_collect_init(&EC,ifp);

    char ExportStepOpt[4096];
    GetValueSk(ifp,"Tracer.step",ExportStepOpt,0,"Range");
    int * export_steps;
    int num_export_steps;
    if(strcasecmp(ExportStepOpt,"range") == 0)
    {
        int s_step = GetValueIk(ifp,"Tracer.step",1,"0");
        int e_step = GetValueIk(ifp,"Tracer.step",2,"0");
        int interval = GetValueIk(ifp,"Tracer.step",3,"0");
        assert(s_step < e_step && interval > 0);
        num_export_steps = (e_step - s_step)/interval + 1;
        export_steps = malloc(sizeof(int)*num_export_steps);
        for(int k=0;k<num_export_steps;++k) export_steps[k] = e_step + k*interval;
    }
    else
    {
        fprintf(stdout,"Unimplemented range opt:%s\n",ExportStepOpt);
        exit(1);
    }
    for(int k=0;k<num_export_steps;++k)
    {
        int new_ejecta_num = load_ejecta_collect(&EC,k);
        fprintf(stdout,"%d step: %d ejecta\n",k,new_ejecta_num);
    }

    ejecta_collect_to_vtp(&EC,EC.output);

    char predictLoc[4096];
    GetValueS(ifp,"Tracer.loops",predictLoc,"none");
    if(strcasecmp("none",predictLoc) != 0)
    {
        for(int k=0;k<200;++k)
        {
            char _vtp_name[1025];
            snprintf(_vtp_name,1024,"%s.%04d.vtp",predictLoc,k+1);
            numerical_ejecta_orbit_moon(&EC,300.0);
            ejecta_collect_to_vtp(&EC,_vtp_name);
            fprintf(stdout,"#");
            fflush(stdout);
        }
    }
    fprintf(stdout,"\n %d ejecta detected\n",EC.len);
    ejecta_collect_test_clean(&EC);
    CloseInputFile(ifp);
    return 1;
}