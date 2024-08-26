//
// Created by li on 8/22/24.
//

#include "VtkWriter.h"
#include <assert.h>
#include <stdio.h>

int main(int argc,char * argv[])
{
    char * txt_name = "./tif2.txt";
    char * vts_name = "./surface2.vts";
    char * ref_nameOuter = "./SPAOuterRing.txt";
    char * ref_nameInner = "./SPAInnerRing.txt";
    FILE * tfp = fopen(txt_name,"r");
    FILE * vfp = fopen(vts_name,"w");
    FILE * rfp1 = fopen(ref_nameOuter, "r");
    FILE * rfp2 = fopen(ref_nameInner, "r");
    assert(tfp != NULL);
    assert(vfp != NULL);
    assert(rfp1 != NULL);
    assert(rfp2 != NULL);

    int nx,ny,nd;
    fscanf(tfp,"%d %d %d",&nx,&ny,&nd);
    float * xt = malloc(sizeof(float)*nx*ny*3);
    float * xtC = malloc(sizeof(float)*nx*ny*3);
    float * dt = malloc(sizeof(float)*nx*ny);
    float Rm = 1.74e6;

    for(int k=0;k<nx;++k)
    {
        for(int j=0;j<ny;++j)
        {
            int lid = j + k*ny;
            fscanf(tfp,"%f %f %f",xt+3*lid+0,xt+3*lid+1,dt+lid);
            xt[3*lid+2] = 0.0;

            float lon = xt[3*lid + 0] * M_PI / 180.0;
            float lat = xt[3*lid + 1] * M_PI / 180.0;

            xtC[3*lid + 0] = Rm*  cos(lat) * cos(lon);
            xtC[3*lid + 1] = Rm*  cos(lat) * sin(lon);
            xtC[3*lid + 2] = Rm* sin(lat) - Rm;

        }
    }

    int ref_n;
    fscanf(rfp1, "%d", &ref_n);
    float * ref_xt = malloc(sizeof(float)*ref_n*3);
    float * ref_xtC = malloc(sizeof(float)*ref_n*3);
    for(int k=0;k<ref_n;++k)
    {
        int id1,id2,id3,id4;
        fscanf(rfp1, "%d,%d,%d,%d,%f,%f", &id1, &id2, &id3, &id4, ref_xt + k * 3, ref_xt + k * 3 + 1);
        ref_xt[k*3+2] = 0.;

        float lon = ref_xt[3*k + 0] * M_PI / 180.0;
        float lat = ref_xt[3*k + 1] * M_PI / 180.0;

        ref_xtC[3*k + 0] = Rm*  cos(lat) * cos(lon);
        ref_xtC[3*k + 1] = Rm*  cos(lat) * sin(lon);
        ref_xtC[3*k + 2] = Rm* sin(lat) - Rm;
    }

    int ref_n2;
    fscanf(rfp1, "%d", &ref_n2);
    float * ref_xt2 = malloc(sizeof(float)*ref_n*3);
    float * ref_xtC2 = malloc(sizeof(float)*ref_n*3);
    for(int k=0;k<ref_n2;++k)
    {
        int id1,id2,id3,id4;
        fscanf(rfp1, "%d,%d,%d,%d,%f,%f", &id1, &id2, &id3, &id4, ref_xt2 + k * 3, ref_xt2 + k * 3 + 1);
        ref_xt2[k*3+2] = 0.;

        float lon = ref_xt2[3*k + 0] * M_PI / 180.0;
        float lat = ref_xt2[3*k + 1] * M_PI / 180.0;

        ref_xtC2[3*k + 0] = Rm*  cos(lat) * cos(lon);
        ref_xtC2[3*k + 1] = Rm*  cos(lat) * sin(lon);
        ref_xtC2[3*k + 2] = Rm* sin(lat) - Rm;
    }


    for(int k=0;k<nx*ny;++k)
    {
        xt[3*k] = 2.0*M_PI*Rm;
        xt[3*k + 1] = 2.0*M_PI*Rm;
        for(int j=0;j<ref_n;++j)
        {
            float *p1 = ref_xtC + 3*j;
            float *p2 = xtC + 3*k;
            float disP12 = 0.0;
            for(int i=0;i<3;++i)
                disP12 += (p1[i]-p2[i])*(p1[i]-p2[i]);
            disP12 = sqrt(disP12);
            if(disP12 < xt[3*k]) xt[3*k] = disP12;
        }

        for(int j=0;j<ref_n2;++j)
        {
            float *p1 = ref_xtC2 + 3*j;
            float *p2 = xtC + 3*k;
            float disP12 = 0.0;
            for(int i=0;i<3;++i)
                disP12 += (p1[i]-p2[i])*(p1[i]-p2[i]);
            disP12 = sqrt(disP12);
            if(disP12 < xt[3*k + 1]) xt[3*k+1] = disP12;
        }
    }


    // write data
    char whole_extent[4096], piece_extent[4096];
    snprintf(whole_extent,4096,"%d %d %d %d 0 0",1,nx,1,ny);
    snprintf(piece_extent,4096,"%d %d %d %d 0 0",1,nx,1,ny);
    vts_file_header(vfp,whole_extent,piece_extent);
    char TimeValueAttr[4096];
    vtk_point_data_header(vfp);
    vtk_dataarrayf(vfp,"h","binary",dt,nx*ny);
    vtk_dataarray_vec_f(vfp,"s","binary",xt,nx*ny,3);
    vtk_point_data_trailer(vfp);
    vtk_cell_data_header(vfp);
    vtk_cell_data_trailer(vfp);
    vtk_output_coordf(vfp,"binary",xtC,nx*ny);
    vts_file_trailer(vfp);

    free(xt);
    free(dt);
    free(xtC);
    free(ref_xt);
    free(ref_xtC);

    fclose(tfp);
    fclose(vfp);
    fclose(rfp1);
    fclose(rfp2);
    return 0;
}