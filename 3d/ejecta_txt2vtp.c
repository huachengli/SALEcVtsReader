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
    ejecta_collect_init(&EC,ifp);

    char ExportStepOpt[4096];
    GetValueSk(ifp,"Ejecta.step",ExportStepOpt,0,"Range");
    double Rm = GetValueD(ifp, "Ejecta.Rm", "1.74e6");
    int * export_steps;
    int num_export_steps;
    if(strcasecmp(ExportStepOpt,"range") == 0)
    {
        int s_step = GetValueIk(ifp,"Ejecta.step",1,"0");
        int e_step = GetValueIk(ifp,"Ejecta.step",2,"0");
        int interval = GetValueIk(ifp,"Ejecta.step",3,"0");
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
    GetValueS(ifp,"Ejecta.predictLoc",predictLoc,"none");
    if(strcasecmp("none",predictLoc) == 0)
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
    else if (strcasecmp("analytical",predictLoc)==0)
    {
        char _vtp_name[1025];
        snprintf(_vtp_name,1024,"%s.%04d.vtp",predictLoc,1);
        fprintf(stdout, "use analytical method to calculate landing site\n");
        analytical_ejecta_orbit_moon(&EC, 1.740e6, 1.622);
        ejecta_collect_to_vtp(&EC,_vtp_name);
    }

    citcoms_sphere * _cs = init_citcoms_sphere3(12, 128, 4);
    const char _prefix[] = "ejecta_txt";
    set_ring_scope(_cs,"inring.txt");
    calculate_ejecta_thickness(_cs, &EC, 1.740e6);
    write_citcoms_sphere(_cs,NULL, _prefix);

    fprintf(stdout,"\n %d ejecta detected\n",EC.len);
    ejecta_collect_test_clean(&EC);
    clean_citcoms_sphere(_cs);

    free(export_steps);
    CloseInputFile(ifp);
    return 1;
}