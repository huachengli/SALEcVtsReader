//
// Created by huach on 26/11/2023.
//
// process *.vtp format tracer (ascii)


#include "VtpReader.h"
#include "VtpTracer.h"
#include "InputParser.h"
#include "VtkWriter.h"
#define MAXNAMELEN 4096

int main(int argc,char * argv[])
{
    char data_inp[MAXNAMELEN] = "sale2d.inp";
    char data_dir[MAXNAMELEN] = ".";
    int maxstep = 1;

    if(argc >= 2)
    {
        maxstep = atoi(argv[1]);
    }

    if(argc >= 3)
    {
        strcpy(data_dir,argv[2]);
    }

    if(argc >= 4)
    {
        strcpy(data_inp,argv[3]);
    }

    char inp_path[MAXNAMELEN] = "";
    snprintf(inp_path,MAXNAMELEN,"%s/%s",data_dir,data_inp);
    GridTracer gTracer;
    GridTracer * gTracer_ptr = &gTracer;
    InputFile * data_ifp = OpenInputFile(inp_path);
    InitGridTracer(gTracer_ptr,data_ifp);
    CloseInputFile(data_ifp);

    char data_name[MAXNAMELEN];
    // load initial position
    snprintf(data_name,MAXNAMELEN,"%s/txt/grid.txt",data_dir);
    LoadGridTxtFile(gTracer_ptr,data_name);

    for(int step=0;step<maxstep;++step)
    {
        gTracer.step = step;
        snprintf(data_name, MAXNAMELEN, "%s/vtp/%s.tracer.proc%%04d.%04d.vtp", data_dir, gTracer.prefix, step);
        VtpTracerCollect * tfcp = OpenVtpTracerCollect(data_name, gTracer.nvtp);
        FlushGridTracerFromVtpCollect(gTracer_ptr,```````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````tfcp```````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````);
        char post_vts_name[MAXNAMELEN];
        snprintf(post_vts_name, MAXNAMELEN, "%s/%s.post.proc%%04d.%04d.vts", data_dir, gTracer.prefix, step);
        WriteGridTracer(gTracer_ptr,post_vts_name);
        CloseVtpTracerCollect(tfcp);
        fprintf(stdout,"writing %s.post.proc%%04d.%04d.vts\n",gTracer.prefix, step);
    }

    return 0;
}

VtpTracerCollect * OpenVtpTracerCollect(const char * _prefix, int _nof)
{
    VtpTracerCollect * tmp = malloc(sizeof(VtpTracerCollect));
    tmp->vtp = malloc(sizeof(VtpFile*)*_nof);
    char _vtpname[MAXNAMELEN];
    tmp->NoF = _nof;
    for(int k=0;k<_nof;++k)
    {
        snprintf(_vtpname,MAXNAMELEN,_prefix,k);
        tmp->vtp[k] = OpenVtpFile(_vtpname);
    }
    return tmp;
}

int CloseVtpTracerCollect(VtpTracerCollect * _vtc)
{
    for(int k=0;k<_vtc->NoF;++k){
        CloseVtpFile(_vtc->vtp[k]);
    }
    free(_vtc->vtp);
    _vtc->vtp = NULL;
    free(_vtc);
    return 1;
}

int InitGridTracer(GridTracer * gtf,InputFile * ifp)
{
    int npgx = GetValueI(ifp,"processor.npgx","1");
    int npgy = GetValueI(ifp,"processor.npgy","1");
    int npx = GetValueI(ifp,"mesh.npx","1");
    int npy = GetValueI(ifp,"mesh.npy","1");
    int pad = 2;
    gtf->nvtp = npgx*npgy;
    gtf->stripe = GetValueI(ifp,"numerical.tarcer_strip","2");
    if(gtf->stripe < 0) gtf->stripe *= -1;
    gtf->nx = npgx*npx + 2*pad;
    gtf->ny = npgy*npy + 2*pad;
    gtf->dtsave = GetValueD(ifp,"output.time","0.1");
    gtf->t0 =  GetValueD(ifp,"cycle.time0","0.0");
    gtf->gz = GetValueDk(ifp,"condition.gravity",1,"-9.8");

    gtf->nx2 = (gtf->nx+gtf->stripe-1)/gtf->stripe;
    gtf->ny2 = (gtf->ny+gtf->stripe-1)/gtf->stripe;

    // malloc memory for Grid
    int ngrid = gtf->nx2*gtf->nx2;
    gtf->ipos = calloc(ngrid*3,sizeof(float));
    gtf->cpos = calloc(ngrid*3,sizeof(float));
    gtf->vel  = calloc(ngrid*3,sizeof(float));
    gtf->den       = calloc(ngrid,sizeof(float));
    gtf->pre       = calloc(ngrid,sizeof(float));
    gtf->mpre      = calloc(ngrid,sizeof(float));
    gtf->tem       = calloc(ngrid,sizeof(float));
    gtf->mtem      = calloc(ngrid,sizeof(float));
    gtf->ejecta_T  = calloc(ngrid,sizeof(float));
    gtf->ejecta_X  = calloc(ngrid,sizeof(float));
    gtf->matid     = calloc(ngrid,sizeof(float));
    gtf->id        = calloc(ngrid,sizeof(float));
    gtf->mask      = calloc(ngrid,sizeof(float));
    gtf->step = 0;
    gtf->len = ngrid;
    GetValueS(ifp,"output.prefix",gtf->prefix,"ParaTest");
    return 1;
}

int LoadGridTxtFile(GridTracer * gtf,const char * fname)
{
    FILE * fp = fopen(fname,"r");
    if(NULL == fp)
    {
        fprintf(stderr,"cannot open %s\n",fname);
        return 0;
    }
    int nx=0,ny=0;
    fscanf(fp,"%d %d\n",&nx,&ny);

    if(nx!=gtf->nx || ny!=gtf->ny)
    {
        fprintf(stdout,"incompatible grid.txt loaded!\n");
        exit(0);
    }

    float * x = calloc(nx, sizeof(float));
    float * y = calloc(ny, sizeof(float));


    for(int k=0;k<nx;++k)
    {
        fscanf(fp,"%f\n",x+k);
    }

    for(int k=0;k<ny;++k)
    {
        fscanf(fp,"%f\n",y+k);
    }

    for(int k=0;k<nx;++k)
    {
        for(int j=0;j<ny;++j)
        {
            if(k%gtf->stripe != 0 || j%gtf->stripe != 0) continue;
            int index = (k/gtf->stripe) + gtf->nx2*(j/gtf->stripe);
            gtf->ipos[index*3 + 0] = x[k];
            gtf->ipos[index*3 + 1] = y[j];
            gtf->ipos[index*3 + 2] = 0.;
        }
    }
    free(x);
    free(y);
    fclose(fp);
    return 1;
}

int FlushGridTracerFromVtp(GridTracer * gtf, VtpFile * vfp)
{

    // ptr need binded
    const int i_pos  = 0;
    const int i_gx   = 1;
    const int i_gy   = 2;
    const int i_id   = 3;
    const int i_matid= 4;
    const int i_vel  = 5;
    const int i_pre  = 6;
    const int i_tem  = 7;
    const int i_den  = 8;
    const int i_mpre = 9;
    const int i_mtem = 10;
    Name2VtpData NameBinded[] = {
            [ 0] = {.name = "coordinate", .data = NULL},
            [ 1] = {.name = "gx", .data = NULL},
            [ 2] = {.name = "gy", .data = NULL},
            [ 3] = {.name = "tag", .data = NULL},
            [ 4] = {.name = "matid", .data = NULL},
            [ 5] = {.name = "vel", .data = NULL},
            [ 6] = {.name = "pre", .data = NULL},
            [ 7] = {.name = "tem", .data = NULL},
            [ 8] = {.name = "den", .data = NULL},
            [ 9] = {.name = "mpre", .data = NULL},
            [10] = {.name = "mtem", .data = NULL},
            [11] = {.name = "unknown", .data = NULL},
    };

    // bind VtpData from VtpFile to named pointer
    for(int k=0;k<vfp->PointNoF;++k)
    {

        for(int j=0;strcasecmp(NameBinded[j].name,"unknown")!=0;++j)
        {
            if(strcasecmp(NameBinded[j].name,vfp->PointField[k].Name) == 0)
            {
                NameBinded[j].data = vfp->PointField + k;
            }
        }
    }

    // flush data into corresponding dataarray
    for(int j=0;strcasecmp(NameBinded[j].name,"unknown")!=0;++j)
    {
        if(NameBinded[j].data == NULL)
        {
            fprintf(stderr,"!%s unfound in vtpfile\n",NameBinded[j].name);
        }
    }

    if(NameBinded[i_mpre].data == NULL) NameBinded[i_mpre].data = NameBinded[i_pre].data;
    if(NameBinded[i_mtem].data == NULL) NameBinded[i_mtem].data = NameBinded[i_tem].data;
    int fresh_deteced = 0;

    for(int k=0;k<vfp->NoP;++k)
    {
        float * pos = NameBinded[i_pos].data->Data + 3*k;
        float * vel = NameBinded[i_vel].data->Data + 2*k;
        int gx  = (int)roundf(NameBinded[i_gx].data->Data[k]);
        int gy  = (int)roundf(NameBinded[i_gy].data->Data[k]);
        float id  = NameBinded[i_id].data->Data[k];
        float pre = NameBinded[i_pre].data->Data[k];
        float tem = NameBinded[i_tem].data->Data[k];
        float den = NameBinded[i_den].data->Data[k];
        float mpre = NameBinded[i_mpre].data->Data[k];
        float mtem = NameBinded[i_mtem].data->Data[k];
        float matid = NameBinded[i_matid].data->Data[k];

        if(gx < 0 || gy< 0)
        {
            // this is a dummy tracer
            continue;
        }

        // convert to index, colum first
        gx = gx/gtf->stripe;
        gy = gy/gtf->stripe;
        int index = gtf->nx2*gy + gx;
        // update variables in gtf
        gtf->cpos[index*3 + 0] = pos[0];
        gtf->cpos[index*3 + 1] = pos[1];
        gtf->cpos[index*3 + 2] = pos[2];

        gtf->vel[index*3 + 0] = vel[0];
        gtf->vel[index*3 + 1] = vel[1];
        gtf->vel[index*3 + 2] = 0.f;

        gtf->matid[index] = matid;
        gtf->id[index]    = id;
        gtf->pre[index] = pre;
        gtf->tem[index] = tem;
        gtf->den[index] = den;

        gtf->mpre[index] = mpre;
        gtf->mtem[index] = mtem;

        if(matid < 0.){
            if(gtf->mask[index] == 0) fresh_deteced++;
            gtf->mask[index] += 1.0f;
            float gz = (float)(-1.0*gtf->gz);
            float land_delay_time = (vel[1] + sqrtf(vel[1]*vel[1] + 2*gz*pos[1]))/gz;
            float eT = gtf->t0 + gtf->dtsave*gtf->step + land_delay_time;
            float eX = pos[0] + land_delay_time*vel[0];
            gtf->ejecta_X[index] = eX;
            gtf->ejecta_T[index] = eT;
        }
    }
    return fresh_deteced;
}

int FlushGridTracerFromVtpCollect(GridTracer * gtf, VtpTracerCollect * tvtcp)
{
    for(int k=0;k<tvtcp->NoF;++k)
    {
        FlushGridTracerFromVtp(gtf,tvtcp->vtp[k]);
    }
}

int WriteGridTracer(GridTracer * gtf, const char * vts_name)
{
    FILE * fp = fopen(vts_name,"w");
    char whole_extent[4096], piece_extent[4096];
    snprintf(whole_extent,4096,"%d %d %d %d 0 0",1,gtf->nx2,1,gtf->ny2);
    snprintf(piece_extent,4096,"%d %d %d %d 0 0",1,gtf->nx2,1,gtf->ny2);
    vts_file_header(fp,whole_extent,piece_extent);
    char TimeValueAttr[4096];
    snprintf(TimeValueAttr,4096,"TimestepValues = \"%f\"",gtf->t0 + gtf->dtsave*gtf->step);
    vtk_point_data_header_with_attr(fp,TimeValueAttr);
    vtk_dataarrayf(fp,"pre","binary",gtf->pre,gtf->len);
    vtk_dataarrayf(fp,"tem","binary",gtf->tem,gtf->len);
    vtk_dataarrayf(fp,"den","binary",gtf->den,gtf->len);
    vtk_dataarrayf(fp,"mpre","binary",gtf->mpre,gtf->len);
    vtk_dataarrayf(fp,"mtem","binary",gtf->mtem,gtf->len);
    vtk_dataarrayf(fp,"matid","binary",gtf->matid,gtf->len);
    vtk_dataarrayf(fp,"mask","binary",gtf->mask,gtf->len);
    vtk_dataarrayf(fp,"id","binary",gtf->id,gtf->len);
    vtk_dataarrayf(fp,"eX","binary",gtf->ejecta_X,gtf->len);
    vtk_dataarrayf(fp,"eT","binary",gtf->ejecta_T,gtf->len);
    vtk_dataarray_vecf(fp,"vel","binary",gtf->vel,gtf->len,3);
    vtk_point_data_trailer(fp);
    vtk_cell_data_header(fp);
    vtk_cell_data_trailer(fp);
    vtk_output_coordf(fp,"binary",gtf->ipos,gtf->len);
    vts_file_trailer(fp);
    fclose(fp);
}
