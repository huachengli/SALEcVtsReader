//
// Created by li on 8/20/24.
//

#include "ejecta_analysis.h"
#include "InputParser.h"
#include <stdio.h>

int main(int argc,char * argv[])
{
//    char get_plane[200];
//    switch(argc)
//    {
//        case 1:
//            strcpy(get_plane,"planes.plot");
//            break;
//        case 2:
//            strcpy(get_plane,argv[1]);
//            break;
//        default:
//            fprintf(stdout,"error in command line\n");
//    }
//
//    InputFile * ifp = OpenInputFile(get_plane);

    ejecta_collect EC;
    ejecta_collect_test_init(&EC);

    for(int k=0;k<400;++k)
    {
        int new_ejecta_num = load_ejecta_collect(&EC,k);
        fprintf(stdout,"%d step: %d ejecta\n",k,new_ejecta_num);
    }

    ejecta_collect_to_vtp(&EC,"ejecta_vtp/test.0000.vtp");
    for(int k=0;k<200;++k)
    {
        char _vtp_name[1025];
        snprintf(_vtp_name,1024,"ejecta_vtp/test.%04d.vtp",k+1);
        numerical_ejecta_orbit_moon(&EC,300.0);
        ejecta_collect_to_vtp(&EC,_vtp_name);
        fprintf(stdout,"#");
        fflush(stdout);
    }
    fprintf(stdout,"\n %d ejecta detected\n",EC.len);
    ejecta_collect_test_clean(&EC);

//    CloseInputFile(ifp);

    return 1;
}