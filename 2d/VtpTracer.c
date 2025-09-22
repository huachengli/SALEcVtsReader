//
// Created by huachengli on 8/30/24.
//
#include "VtpReader.h"
#include "VtpTracer.h"
#include "InputParser.h"
#include "VtkWriter.h"
#include <omp.h>
#include <unistd.h>
#define MAXNAMELEN 4096

int main(int argc,char * argv[])
{
    char data_inp[MAXNAMELEN] = "sale2d.inp";
    char data_dir[MAXNAMELEN] = ".";
    int maxstep = 5;

    int c;
    int write_vts = 0;
    opterr = 0;
    while ((c = getopt (argc, argv, "n:v:f:d:")) != -1)
    {
        switch (c)
        {
            case 'n':
                maxstep = atoi(optarg);
                break;
            case 'v':
                write_vts = atoi(optarg);
                break;
            case 'f':
                strcpy(data_inp,optarg);
                break;
            case 'd':
                strcpy(data_dir,optarg);
                break;
            case '?':
                if(optopt == 'c')
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                else if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf (stderr,"Unknown option character `\\x%x'.\n",optopt);
                return 1;
            default:
                abort ();
        }
    }

    char inp_path[MAXNAMELEN*2] = "";
    snprintf(inp_path,MAXNAMELEN*2,"%s/%s",data_dir,data_inp);
    GridTracer gTracer;
    GridTracer * gTracer_ptr = &gTracer;
    InputFile * data_ifp = OpenInputFile(inp_path);
    InitGridTracer(gTracer_ptr,data_ifp);
    CloseInputFile(data_ifp);

    char data_name[MAXNAMELEN*2];
    // load initial position
    snprintf(data_name,MAXNAMELEN*2,"%s/txt/grid.txt",data_dir);
    LoadGridTxtFile(gTracer_ptr,data_name);

    for(int step=0;step<maxstep;++step)
    {
        fprintf(stdout,"processing step %d (",step);
        gTracer.step = step;
        snprintf(data_name, MAXNAMELEN*2, "%s/vtp/%s.tracer.proc%%04d.%04d.vtp", data_dir, gTracer.prefix, step);
        VtpTracerCollect * tfcp = FlushVtpTracerCollect(gTracer_ptr,data_name, gTracer.nvtp);
        char post_vts_name[MAXNAMELEN*2];
        snprintf(post_vts_name, MAXNAMELEN*2, "%s/post/%s.post.%04d.vts", data_dir, gTracer.prefix, step);

        if(write_vts > 0 && step%write_vts == 0)
        {
            WriteGridTracer(gTracer_ptr,post_vts_name);
            // fprintf(stdout," post/vts ");
        }

        snprintf(post_vts_name, MAXNAMELEN*2, "%s/post/%s", data_dir, gTracer.prefix);
        // fprintf(stdout," post/bin ");

        ExportGridTracerF32Bin(gTracer_ptr,post_vts_name);
        CloseVtpTracerCollect(tfcp);
        fprintf(stdout,")\n");
    }

    return 0;
}