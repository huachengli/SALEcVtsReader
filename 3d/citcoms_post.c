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
    load_citcoms_step(_cdata,0);

    SphereIntegrateCitcomsDump2(_cdata,_cdata->OutPrefix);

    clean_citcoms_data(_cdata);
    close_citcoms_data(_cdata);

    // FILE * fp = fopen("/public/home/huachengli/exec-citcoms/V2VD0IC350H100-job26/impact_spa_high_alumina/a.proc95.980.vts","r");
    // unsigned char LineBuffer[1024];
    // int k = 0;
    // while(1)
    // {
    //     k++;
    //     int l = ReadLineTrim(LineBuffer,fp);
    //     if(l<=0)
    //         break;
    //     char * r0 = strstr(LineBuffer,"DataArray");
    //     char * r1 = strstr(LineBuffer,"/DataArray");
    //     if(r0 != NULL && r1 == NULL)
    //     {
    //         float * data;
    //         unsigned long ndata;
    //         ReadVtsAsciiF32(&data,&ndata,fp);
    //         fprintf(stdout,"DATA(%d):%f,%f,...,%f,%f\n",ndata, data[0], data[1], data[ndata-2],data[ndata-1]);
    //         free(data);
    //     }
    //
    //     if(r0 != NULL && r1 != NULL)
    //     {
    //         fprintf(stdout,"%d:%s\n",k,LineBuffer);
    //     }
    // }
    // fclose(fp);

    // InputFile * ifp = OpenInputFile(inp_file);
    // citcoms_dump * cdp = InitCitcomsDump(ifp);
    //
    // citcoms_check_tracer_element(cdp);
    //
    // CloseCitcomsDump(cdp);
    // CloseInputFile(ifp);
    return 0;
}