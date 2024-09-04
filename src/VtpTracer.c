//
// Created by huach on 26/11/2023.
//
// process *.vtp format tracer (ascii)


#include "VtpReader.h"
#include "VtpTracer.h"
#include "InputParser.h"
#include "VtkWriter.h"
#include <omp.h>
#include <unistd.h>
#define MAXNAMELEN 4096

VtpTracerCollect * OpenVtpTracerCollect(const char * _prefix, int _nof)
{
    VtpTracerCollect * tmp = malloc(sizeof(VtpTracerCollect));
    tmp->vtp = malloc(sizeof(VtpFile*)*_nof);
    tmp->NoF = _nof;
    #pragma omp parallel for num_threads(4) shared(_nof,_prefix,tmp) default(none)
    for(int k=0;k<_nof;++k)
    {
        char _vtpname[MAXNAMELEN];
        snprintf(_vtpname,MAXNAMELEN,_prefix,k);
        tmp->vtp[k] = OpenVtpFile(_vtpname);
    }
    return tmp;
}

VtpTracerCollect * OpenSALEcTracerCollect(InputFile * ifp, int step)
{
    char SALEcInp[4096];
    GetValueS(ifp,"SALEc.input",SALEcInp,"SALEc.inp");
    InputFile * sifp = OpenInputFile(SALEcInp);
    int npgx = GetValueI(sifp,"processor.npgx","2");
    int npgy = GetValueI(sifp,"processor.npgy","2");
    int npgz = GetValueI(sifp,"processor.npgz","1");
    int _nof = npgx*npgy*npgz;
    CloseInputFile(sifp);

    char _prefix[MAXNAMELEN];
    GetValueS(ifp,"tracer.data",_prefix,"bm");

    char tmp_prefix[MAXNAMELEN*2];
    snprintf(tmp_prefix,MAXNAMELEN*2-1,"%s.proc%%d.%d.vtp",_prefix,step);

    VtpTracerCollect * x =  OpenVtpTracerCollect(tmp_prefix,_nof);
    x->step = step;
    return x;
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

int ShowBriefVtpFile(VtpFile * vfp, FILE * fp)
{
    if(vfp==NULL)
    {
        return 0;
    }

    fprintf(fp,"Name:%s\n",vfp->name);
    fprintf(fp,"  points:%d\n",vfp->NoP);
    fprintf(fp,"Field list:\n");
    for(int k=0;k < vfp->PointNoF;++k)
    {
        fprintf(fp,"%d|%s|\n",k,vfp->PointField[k].Name);
    }
    return vfp->NoP;
}

int ShowBriefVtpColleect(VtpTracerCollect * vtc, FILE * fp)
{
    fprintf(fp,"<<<<[collect]%s\n",vtc->prefix);
    int sum_nop = 0;
    for(int k=0;k<vtc->NoF;++k)
    {
        sum_nop += ShowBriefVtpFile(vtc->vtp[k],fp);
    }
    fprintf(fp,"%d pts\n",sum_nop);
    fprintf(fp,">>>>[*]%s\n",vtc->prefix);
    return sum_nop;
}

VtpFile * SALEcVtpCollectFilter(VtpTracerCollect * vtc)
{

    int sum_nop_filtered = 0;
    unsigned long matFieldId = 100;
    for(int k=0;k<vtc->NoF;++k)
    {
        if(vtc->vtp[k]->NoP == 0)
            continue;
        if(matFieldId == 100)
            matFieldId = find_vtpfield("matid",vtc->vtp[k]);
        VtpData * vtpx = vtc->vtp[k]->PointField + matFieldId;
        for(int j=0;j<vtpx->DataLen;++j)
        {
            if(fabs(fabs(vtpx->Data[j]) - 1.0) < 0.1) sum_nop_filtered ++;
        }
    }

    VtpFile * tmp = NULL;
    unsigned long pindex = 0;

    for(int k=0;k<vtc->NoF;++k)
    {
        if(vtc->vtp[k]->NoP == 0)
            continue;

        if(NULL == tmp) tmp = duplicate_vtp(vtc->vtp[k],sum_nop_filtered);
        VtpFile * vtp = vtc->vtp[k];
        VtpData * vtpx = vtc->vtp[k]->PointField + matFieldId;
        for(int j=0;j<vtpx->DataLen;++j)
        {
            if(fabs(fabs(vtpx->Data[j]) - 1.0) < 0.1){
                copy_vtp_k(tmp,pindex,vtp,j);
                pindex++;
            }
        }
    }
    VtpCoordinateReshape(tmp);
    return tmp;
}

unsigned long find_vtpfield(const char * _src, VtpFile * vfp)
{
    unsigned long result= 100;
    for(int k=0;k<vfp->PointNoF;++k)
    {
        if(0== strcasecmp(_src,vfp->PointField[k].Name))
        {
            result = k;
            break;
        }
    }
    return result;
}

VtpFile * duplicate_vtp(VtpFile * in, unsigned long n)
{
    VtpFile * tmp = malloc(sizeof(VtpFile));
    tmp->NoP = n;
    tmp->PointNoF = in->PointNoF;
    if(n == 0) return tmp;

    for(int k=0;k<tmp->PointNoF;++k)
    {
        tmp->PointField[k].DataLen = in->PointField[k].NoC*n;
        tmp->PointField[k].NoC = in->PointField[k].NoC;
        tmp->PointField[k].Data =  malloc(sizeof(VTSDATAFLOAT)*tmp->PointField[k].DataLen);
        strcpy(tmp->PointField[k].Name,in->PointField[k].Name);
        strcpy(tmp->PointField[k].Type,in->PointField[k].Type);
        strcpy(tmp->PointField[k].Format,in->PointField[k].Format);
    }
    return tmp;
}

int copy_vtp_k(VtpFile * x, unsigned long px, VtpFile * y, unsigned long  py)
{
    if(x->PointNoF != y->PointNoF) return 0;
    for(int k=0;k<x->PointNoF;++k)
    {
        memcpy(x->PointField[k].Data + px*x->PointField[k].NoC,
               y->PointField[k].Data + py*y->PointField[k].NoC, sizeof(VTSDATAFLOAT)*x->PointField[k].NoC);
    }
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
    double projectile_radii = GetValueD(ifp,"projectile.radiu","1.0");
    double z_threshold =  GetValueD(ifp, "ejecta.threshold", "2.");
    gtf->z_up = 3.0f*z_threshold*projectile_radii;
    gtf->z_low = z_threshold*projectile_radii;

    gtf->nx2 = (gtf->nx+gtf->stripe-1)/gtf->stripe;
    gtf->ny2 = (gtf->ny+gtf->stripe-1)/gtf->stripe;

    // malloc memory for Grid
    int ngrid = gtf->nx2*gtf->ny2;
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
    gtf->ejecta_U  = calloc(ngrid,sizeof(float));
    gtf->ejecta_V  = calloc(ngrid,sizeof(float));
    gtf->ejecta_t  = calloc(ngrid,sizeof(float));
    gtf->ejecta_x  = calloc(ngrid,sizeof(float));
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
            fprintf(stderr,"!%s unknown in vtpfile\n",NameBinded[j].name);
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

        gtf->ejecta_num = 0;
        if(matid < 0.){
            gtf->ejecta_num ++;
            if(gtf->mask[index] == 0) fresh_deteced++;
            if(gtf->mask[index] >= 2.0) continue; // this tracer have been tracked with enough time step.
            gtf->mask[index] += 1.0f;
            float gz = (float)(-1.0*gtf->gz);
            float land_delay_time = (vel[1] + sqrtf(vel[1]*vel[1] + 2*gz*pos[1]))/gz;
            float launch_time = (vel[1] - sqrtf(vel[1]*vel[1] + 2*gz*pos[1]))/gz;
            float eT = gtf->t0 + gtf->dtsave*gtf->step + land_delay_time;
            float eX = pos[0] + land_delay_time*vel[0];
            gtf->ejecta_X[index] = eX;
            gtf->ejecta_T[index] = eT;
            gtf->ejecta_U[index] = vel[0];
            gtf->ejecta_V[index] = vel[1] - gz*launch_time;
            gtf->ejecta_t[index] = gtf->t0 + gtf->dtsave*gtf->step + launch_time;
            gtf->ejecta_x[index] = pos[0] + launch_time*vel[0];
        }
    }
    return fresh_deteced;
}

VtpTracerCollect * FlushVtpTracerCollect(GridTracer * gtf,const char * _prefix, int _nof)
{
    VtpTracerCollect * tmp = malloc(sizeof(VtpTracerCollect));
    tmp->vtp = malloc(sizeof(VtpFile*)*_nof);
    tmp->NoF = _nof;
#pragma omp parallel for num_threads(8) shared(_nof,_prefix,tmp,gtf) default(none)
    for(int k=0;k<_nof;++k)
    {
        char _vtpname[MAXNAMELEN];
        snprintf(_vtpname,MAXNAMELEN,_prefix,k);
        tmp->vtp[k] = OpenVtpFile(_vtpname);
        FlushGridTracerFromVtp(gtf,tmp->vtp[k]);
    }
    return tmp;
}

int FlushGridTracerFromVtpCollect(GridTracer * gtf, VtpTracerCollect * tvtcp)
{
    for(int k=0;k<tvtcp->NoF;++k)
    {
        FlushGridTracerFromVtp(gtf,tvtcp->vtp[k]);
    }
    return 0;
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
    return 0;
}

int f32ArrayTofile(float * arr, int n, const char fname[])
{
    FILE * fp = fopen(fname,"wb");
    if(NULL==fp) return -1;
    fwrite(arr, sizeof(float),n,fp);
    fclose(fp);
    return n;
}

int ExportGridTracerF32Bin(GridTracer * gtf, const char * binprefix)
{
    int ngrid = gtf->nx2*gtf->ny2;
    int ejecta_num_hold = 0;
    // scan all the recorded grid, some ejecta would be removed from computation
    for(int k=0;k<ngrid;++k)
    {
        if(gtf->mask[k] > 0) ejecta_num_hold++;
    }
    if(ejecta_num_hold < gtf->ejecta_num)
    {
        fprintf(stdout,"tracers(%d) recorded by mask[k] is less than ejecta_num(%d)\n",ejecta_num_hold,gtf->ejecta_num);
        exit(0);
    }

    if(ejecta_num_hold <= 0) return 0;
    float * tmp_eX = calloc(ejecta_num_hold,sizeof(float));
    float * tmp_eT = calloc(ejecta_num_hold,sizeof(float));
    float * tmp_ix = calloc(ejecta_num_hold,sizeof(float));
    float * tmp_iy = calloc(ejecta_num_hold,sizeof(float));
    float * tmp_eU  = calloc(ejecta_num_hold,sizeof(float));
    float * tmp_eV  = calloc(ejecta_num_hold,sizeof(float));
    float * tmp_et  = calloc(ejecta_num_hold,sizeof(float));
    float * tmp_ex  = calloc(ejecta_num_hold,sizeof(float));


    int i = 0;
    for(int k=0;k<ngrid;++k)
    {
        if(gtf->mask[k] <= 0) continue;
        tmp_eT[i] = gtf->ejecta_T[k];
        tmp_eX[i] = gtf->ejecta_X[k];
        tmp_ix[i] = gtf->ipos[3*k + 0];
        tmp_iy[i] = gtf->ipos[3*k + 1];
        tmp_eU[i] = gtf->ejecta_U[k];
        tmp_eV[i] = gtf->ejecta_V[k];
        tmp_et[i] = gtf->ejecta_t[k];
        tmp_ex[i] = gtf->ejecta_x[k];
        i++;
    }

    char tmp_name[MaxStrLen];
    snprintf(tmp_name,MaxStrLen,"%s.ix.%04d.bin",binprefix,gtf->step);
    f32ArrayTofile(tmp_ix,ejecta_num_hold,tmp_name);
    snprintf(tmp_name,MaxStrLen,"%s.iy.%04d.bin",binprefix,gtf->step);
    f32ArrayTofile(tmp_iy,ejecta_num_hold,tmp_name);
    snprintf(tmp_name,MaxStrLen,"%s.eX.%04d.bin",binprefix,gtf->step);
    f32ArrayTofile(tmp_eX,ejecta_num_hold,tmp_name);
    snprintf(tmp_name,MaxStrLen,"%s.eT.%04d.bin",binprefix,gtf->step);
    f32ArrayTofile(tmp_eT,ejecta_num_hold,tmp_name);
    snprintf(tmp_name,MaxStrLen,"%s.eU.%04d.bin",binprefix,gtf->step);
    f32ArrayTofile(tmp_eU,ejecta_num_hold,tmp_name);
    snprintf(tmp_name,MaxStrLen,"%s.eV.%04d.bin",binprefix,gtf->step);
    f32ArrayTofile(tmp_eV,ejecta_num_hold,tmp_name);
    snprintf(tmp_name,MaxStrLen,"%s.et.%04d.bin",binprefix,gtf->step);
    f32ArrayTofile(tmp_et,ejecta_num_hold,tmp_name);
    snprintf(tmp_name,MaxStrLen,"%s.ex.%04d.bin",binprefix,gtf->step);
    f32ArrayTofile(tmp_ex,ejecta_num_hold,tmp_name);


    free(tmp_eX);
    free(tmp_eT);
    free(tmp_ix);
    free(tmp_iy);
    free(tmp_eU);
    free(tmp_eV);
    free(tmp_et);
    free(tmp_ex);


    return ejecta_num_hold;
}
