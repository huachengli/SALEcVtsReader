//
// Created by lhc on 25/4/2024.
//
// read vts data and generate profile automatically
#include "Utility.h"
#include <stdio.h>
#include <float.h>
#include "Utility2d.h"
#include <unistd.h>
#include "lmath.h"
#include "VtpTracer.h"

void Vts2dProfile(SALEcData * _sdata, SALEc2dProfile * _spdata);
void Vts2dContour(SALEcData * _sdata, SALEc2dProfile * _spdata, const char * field_name, double level_z, int use_point);
void Write2dProfile(SALEc2dProfile * _spdata, const char * fname, const char * format);
void Clean2dProfile(SALEc2dProfile * _spdata);
int profile_vtp_filter(const double * pos, void * ctx);
int load_vtp_grid(VtpFile * vfp, const char * fname);
typedef struct {
    SALEc2dProfile * p;
    double t;
} VtpFilterCtx;
int main(int argc,char * argv[])
{

    char SALEcDataDir[256] = ".";
    char SALEcInpName[256] = "sale2d.inp";
    int Expstep = 1;

    int c;
    int f_set_flag = 0;
    char f_set_name[256];
    int f_bin_out = 0;
    int f_use_contour = 1;
    double p_expose_tracer = 0.;
    char p_field_name[256] = "e_den";
    double p_field_threshold = 200.0;
    char o_set_name[256] = "profile2d";
    opterr = 0; // defined in getopt*.h
    const char usage_info[] = "-OPTION [ARGS]        DEFAULT\n"
                              "-d [data directory]   $PWD   \n"
                              "-f [input name]       salec2d.inp\n"
                              "-e [exported step]    1\n"
                              "-o [output name]      profile_2d\n"
                              "-b                    no default, write binary output\n"
                              "-v [field name]       e_den\n"
                              "-t [field threshold]  200.0\n"
                              "-r [tracer filter scope]";
    while ((c = getopt (argc, argv, "d:f:e:o:bv:t:c:r:h")) != -1)
    {
        switch (c)
        {
            case 'd':
                strcpy(SALEcDataDir,optarg);
                break;
            case 'f':
                strcpy(f_set_name,optarg);
                f_set_flag = 1;
                break;
            case 'e':
                Expstep = atoi(optarg);
                break;
            case 'b':
                f_bin_out = 1;
                break;
            case 'v':
                strcpy(p_field_name,optarg);
                break;
            case 't':
                p_field_threshold = atof(optarg);
                break;
            case 'r':
                p_expose_tracer = atof(optarg);
                break;
            case 'c':
                f_use_contour = atoi(optarg);
                break;
            case 'o':
                strcpy(o_set_name,optarg);
                break;
            case 'h':
                fprintf(stdout,"USAGE:\n%s\n",usage_info);
                break;
            case '?':
                fprintf (stderr,"Unknown option character `\\x%x'.\n",optopt);
                fprintf(stdout,"USAGE:\n%s\n",usage_info);
                return 1;
            default:
                fprintf(stdout,"USAGE:\n%s\n",usage_info);
                abort();
        }
    }

    if(f_set_flag == 0)
    {
        sprintf(SALEcInpName,"%s/sale2d.inp",SALEcDataDir);
    } else
    {
        sprintf(SALEcInpName,"%s/%s",SALEcDataDir,f_set_name);
    }

    // load vts data
    // load the reference data
    SALEcData * SaleData = (SALEcData *) malloc(sizeof(SALEcData));
    SALEcPlanetInfo * SaleInfo = malloc(sizeof(SALEcPlanetInfo));
    Load2dInpInfo(SaleData, SaleInfo, SALEcInpName);
    char DataPath[256];
    snprintf(DataPath,256,"%s/vts/%s.proc%%04d.%04d.vts",SALEcDataDir,SaleInfo->prefix,Expstep);
    LoadVts2dData(SaleData, DataPath);

    SALEc2dProfile profile = {.n=0,.pos=NULL,.pos2=NULL};

    if(f_use_contour)
    {
        Vts2dContour(SaleData, &profile, p_field_name, p_field_threshold, 1);
        strcat(o_set_name,".csv");
    }
    else
    {
        Vts2dProfile(SaleData,&profile);
        strcat(o_set_name,"c0.csv");
    }

    // write profile to csv files
    Write2dProfile(&profile,o_set_name,"csv");
    Write2dProfile(&profile,o_set_name,"vtp");

    if(f_bin_out)
    {
        strcat(o_set_name,".bin");
        Write2dProfile(&profile,o_set_name,"bin");
    }

    // process tracer data, check the excavation depth
    if(p_expose_tracer > 0.0)
    {
        char VtpPath[256];
        snprintf(VtpPath,256,"%s/vtp/%s.tracer.proc%%04d.%04d.vtp",SALEcDataDir,SaleInfo->prefix,Expstep);
        VtpTracerCollect * Vtc = OpenVtpTracerCollect(VtpPath, SaleData->VtsBlockNum);
        VtpFilterCtx ctx = {.p = &profile, .t=p_expose_tracer};

        VtpFile * Vtpf = SALEcVtpCollectPosFuncFilter2(Vtc,profile_vtp_filter,&ctx);
        for(int s= strlen(o_set_name)-1; s>=0; s--)
        {
            if(o_set_name[s] == '.')
            {
                o_set_name[s] = '\0';
                break;
            }
        }
        strcat(o_set_name,".tracer.vtp");
        strcpy(Vtpf->name,o_set_name);

        // replace gx,gy to real coordinate
        char grid_txt_name[512];
        snprintf(grid_txt_name, 512, "%s/txt/grid.txt", SALEcDataDir);
        load_vtp_grid(Vtpf, grid_txt_name);

        WriteVtpFile(Vtpf);
        WriteVtpTxt(Vtpf);
        CloseVtpFile(Vtpf);
        CloseVtpTracerCollect(Vtc);
    }


    CleanSALEcData(SaleData);
    Clean2dProfile(&profile);
    free(SaleData);

    // logs
    fprintf(stdout,"\nwrite profile[%d] to %s\n",Expstep,o_set_name);
    return 0;
}

void Vts2dProfile(SALEcData * _sdata, SALEc2dProfile * _spdata)
{
    //allocate mem
    _spdata->n = _sdata->nGCLC[0];
    _spdata->pos = (VTSDATAFLOAT *) malloc(sizeof(VTSDATAFLOAT)*_spdata->n);
    _spdata->pos2= (VTSDATAFLOAT *) malloc(sizeof(VTSDATAFLOAT)*_spdata->n);
    _spdata->seg = NULL;

    // get density field id
    unsigned long VofId = find_cellfield("e_vof",_sdata->VSF);
    unsigned long DenId = find_cellfield("e_den",_sdata->VSF);
    const vts_float DenThreshold = 200.0f;
    const vts_float VofThreshold = 0.5f;

    #pragma omp parallel for num_threads(OMP2D_THREADS) default(shared)
    for(int i=0;i<_sdata->nGCLC[0];++i)
    {
        _spdata->pos[i] = _sdata->GCLC[0][i]; // x coordinate of contour points
        _spdata->pos2[i] = _sdata->GCLC[1][0]; // y coordinate of contour points

        vts_float vof_t=0, vof_b=0;
        vts_float den_t=0, den_b=0;
        for(int j = 0; j < _sdata->nGCLC[1]; ++j) {
            vts_float *vof_f = Vtm2dGetCellData(_sdata, VofId, i, j);
            vts_float *den_f = Vtm2dGetCellData(_sdata, DenId, i, j);

            vof_t = vof_f[0];
            den_t = den_f[0];

            if(j >= 2 && vof_t > VofThreshold)
            {
                vts_float w0 = (vof_t - VofThreshold)/(vof_t - vof_b);
                vts_float w1 = (VofThreshold - vof_b)/(vof_t - vof_b);
                _spdata->pos2[i] = _sdata->GCLC[1][j]*w1 + _sdata->GCLC[1][j-1]*w0;
                break;
            }

            if(j >= 2 && den_t < DenThreshold)
            {
                vts_float w0 = (den_t - DenThreshold)/(den_t - den_b);
                vts_float w1 = (DenThreshold - den_b)/(den_t - den_b);
                _spdata->pos2[i] = _sdata->GCLC[1][j]*w1 + _sdata->GCLC[1][j-1]*w0;
                break;
            }

            vof_b = vof_f[0];
            den_b = den_f[0];
            _spdata->pos2[i] = _sdata->GCLC[1][j];
        }
    }
}

void Vts2dContour(SALEcData * _sdata, SALEc2dProfile * _spdata, const char * field_name, double level_z, int use_point_data)
{
    // check data id (only in cell)
    unsigned fId = find_cellfield(field_name,_sdata->VSF);
    if(UNKNOWNFIELD == fId)
    {
        fId = find_cellfield("e_den",_sdata->VSF);
    }

    // prepare data for Contour(CONREC)
    unsigned num_x_cell = _sdata->nGCLC[0];
    unsigned num_y_cell = _sdata->nGCLC[1];
    unsigned num_cell = num_x_cell * num_y_cell;

    double * cell_data_mem = malloc(sizeof(double)*num_cell);
    double ** cell_data = malloc(sizeof(double*)*num_x_cell);


    // prepare coordinate of cell
    double * cell_x = malloc(sizeof(double)*num_x_cell);
    double * cell_y = malloc(sizeof(double)*num_y_cell);
    assert(cell_x!=NULL);
    assert(cell_y!=NULL);
    for(unsigned i=0; i<num_x_cell; ++i)
    {
        cell_x[i] = _sdata->GCLC[0][i];
    }
    for(unsigned j=0;j<num_y_cell;++j)
    {
        cell_y[j] = _sdata->GCLC[1][j];
    }
    // copy data into cell_data
    for(unsigned i=0; i<num_x_cell; ++i)
    {
        cell_data[i] = cell_data_mem + i*num_y_cell;
        for(unsigned j=0;j<num_y_cell;++j)
        {
            cell_data[i][j] = (double) Vtm2dGetCellData(_sdata, fId, i, j)[0];
        }
    }

    int num_pts[2] = {0};
    double  * con_pts = NULL;
    int * con_seg = NULL;

    if(use_point_data == 0)
    {
        Contour(cell_data,0,num_x_cell-1,0,num_y_cell-1,cell_x,cell_y,level_z,&con_pts,num_pts);
    }
    else
    {
        // interpolate to vertex points
        unsigned num_point = (num_x_cell + 1)*(num_y_cell + 1);
        double * point_data_mem = malloc(sizeof(double)*num_point);
        double ** point_data = malloc(sizeof(double*)*(num_x_cell+1));
        for(int i=0; i<num_x_cell+1;++i)
        {
            point_data[i] = point_data_mem + i*(num_y_cell + 1);
            int l = MAX(0,i-1);
            int r = MIN(i,num_x_cell-1);
            for(int j=0;j<num_y_cell+1;++j)
            {
                int t = MIN(j,num_y_cell-1);
                int b = MAX(0,j-1);
                point_data[i][j] = 0.25*(cell_data[l][b] + cell_data[l][t] + cell_data[r][b] + cell_data[r][t]);
            }
        }

        double * point_x = malloc(sizeof(double)*num_x_cell);
        double * point_y = malloc(sizeof(double)*num_y_cell);
        assert(point_x!=NULL && point_y!=NULL);
        for(unsigned i=0; i<num_x_cell+1; ++i)
        {
            point_x[i] = _sdata->GCLV[0][i];
        }
        for(unsigned j=0;j<num_y_cell+1;++j)
        {
            point_y[j] = _sdata->GCLV[1][j];
        }

        Contour(point_data,0,num_x_cell,0,num_y_cell,point_x,point_y,level_z,&con_pts,num_pts);

        free(point_x);
        free(point_y);
        free(point_data);
        free(point_data_mem);
    }

    double tol = fabs(cell_x[0] - cell_x[1])*1e-2;
    Contour_pts_sort(con_pts,num_pts,&con_seg,tol);

    _spdata->n = num_pts[0];
    _spdata->pos = (VTSDATAFLOAT *) malloc(sizeof(VTSDATAFLOAT)*_spdata->n);
    _spdata->pos2= (VTSDATAFLOAT *) malloc(sizeof(VTSDATAFLOAT)*_spdata->n);
    _spdata->seg = con_seg;

    for(int k=0; k<_spdata->n;++k)
    {
        _spdata->pos[k] = (VTSDATAFLOAT) con_pts[2*k + 0];
        _spdata->pos2[k] =(VTSDATAFLOAT) con_pts[2*k + 1];
    }

    free(con_pts);
    free(cell_x);
    free(cell_y);
    free(cell_data);
    free(cell_data_mem);
}


void WriteProfileVtk(SALEc2dProfile * _spdata, const char * fname)
{
    int n = _spdata->n;
    const vts_float * x = _spdata->pos;
    const vts_float * y = _spdata->pos2;

    int m = 0;
    const int * seg = NULL;
    if(NULL != _spdata->seg)
    {
        m = _spdata->seg[0];
        seg = _spdata->seg + 1;
    }


    if(m > 0)
    {
        int sum = 0;
        for (int i = 0; i < m; ++i) sum += seg[i];
        if (sum != n) {
            fprintf(stderr, "ERROR: sum(seg) != n\n");
            return;
        }
    }

    FILE * fp = fopen(fname, "w");
    assert(NULL!=fp);

    fprintf(fp,
            "<?xml version=\"1.0\"?>\n"
            "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
            "  <PolyData>\n"
            "    <Piece NumberOfPoints=\"%d\" NumberOfLines=\"%d\">\n"
            "      <Points>\n"
            "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n",
            n, m>0?m:1);

    for (int i = 0; i < n; ++i)
        fprintf(fp, "%g %g 0\n", x[i], y[i]);

    fprintf(fp,
            "        </DataArray>\n"
            "      </Points>\n"
            "      <Lines>\n");

    fprintf(fp,
            "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
    for (int i = 0; i < n; ++i) fprintf(fp, "%d ", i);
    fprintf(fp, "\n        </DataArray>\n");

    fprintf(fp,
            "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
    if(m>0)
    {
        int off = 0;
        for (int i = 0; i < m; ++i) {
            off += seg[i];
            fprintf(fp, "%d ", off);
        }
    }
    else
    {
        fprintf(fp, "%d ", n);
    }
    fprintf(fp, "\n        </DataArray>\n"
                "      </Lines>\n"
                "    </Piece>\n"
                "  </PolyData>\n"
                "</VTKFile>\n");

    fclose(fp);
    printf("VTK file written: %s\n", fname);
}


void Write2dProfile(SALEc2dProfile * _spdata, const char * fname, const char * format)
{
    if(strcasecmp(format,"csv") == 0)
    {
        FILE * fp = fopen(fname,"w");
        assert(fp!=NULL);

        if(NULL == _spdata->seg)
        {
            fprintf(fp,"%10d,%10d\n", 1, _spdata->n);
            fprintf(fp,"%10d,%10d\n", 0, _spdata->n);
        }
        else
        {
            fprintf(fp,"%10d,%10d\n", _spdata->seg[0], _spdata->n);
            for(int k=1;k<=_spdata->seg[0];++k)
            {
                fprintf(fp,"%10d,%10d\n", k-1, _spdata->seg[k]);
            }
        }

        for(int k=0;k<_spdata->n;++k)
        {
            fprintf(fp,"%10.5e, %10.5e\n",_spdata->pos[k],_spdata->pos2[k]);
        }
        fclose(fp);
    }

    if(strcasecmp(format,"bin") == 0)
    {
        FILE * fp = fopen(fname,"wb");
        fwrite(_spdata->pos, sizeof(VTSDATAFLOAT),_spdata->n, fp);
        fwrite(_spdata->pos2, sizeof(VTSDATAFLOAT),_spdata->n, fp);
        fclose(fp);
    }

    if(strcasecmp(format,"vtp")==0)
    {
        char vtk_name[512];
        strcpy(vtk_name,fname);
        for(size_t s=strlen(fname)-1; s>=0; --s)
        {
            if(fname[s] == '.')
            {
                vtk_name[s] = '\0';
                break;
            }
        }
        strcat(vtk_name,".vtp");
        WriteProfileVtk(_spdata, vtk_name);
    }
}

void Clean2dProfile(SALEc2dProfile * _spdata)
{
    // clean _spdata
    free(_spdata->pos2);
    free(_spdata->pos);
    if(NULL != _spdata->seg)
        free(_spdata->seg);
}

typedef struct { double y; } YPt;
int cmp_yp(const void *a, const void *b)
{
    double d = ((const YPt*)a)->y - ((const YPt*)b)->y;
    return (d>0)-(d<0);
}

double dist_point_seg(double px, double py, double x1, double y1, double x2, double y2)
{
    double dx = x2 - x1;
    double dy = y2 - y1;
    if (dx == 0 && dy == 0) return hypot(px-x1, py-y1);
    double t = ((px-x1)*dx + (py-y1)*dy) / (dx*dx + dy*dy);
    t = t < 0 ? 0 : (t > 1 ? 1 : t);
    double cx = x1 + t*dx;
    double cy = y1 + t*dy;
    return hypot(px-cx, py-cy);
}

int side_of_polyline(int n,const float *x, const float *y,double a, double b, double *dist)
{
    /*
     * -1 : lower side of polyline
     * 0  : on polyline
     * 1  : upper side of polyline
     */
    YPt ylist[1024];
    int cnt = 0;
    double min_dist = DBL_MAX;

    for (int i = 0; i < n-1; ++i) {
        double x1 = x[i], y1 = y[i];
        double x2 = x[i+1], y2 = y[i+1];

        double d = dist_point_seg(a, b, x1, y1, x2, y2);
        if (d < min_dist) min_dist = d;

        if (x1 == x2) {
            if (x1 == a) {
                ylist[cnt++].y = fmin(y1, y2);
                ylist[cnt++].y = fmax(y1, y2);
            }
            continue;
        }
        if ((x1-a)*(x2-a) > 0) continue;
        double t = (a - x1) / (x2 - x1);
        double y0 = y1 + t*(y2 - y1);
        ylist[cnt++].y = y0;
    }
    *dist = min_dist;
    if (cnt == 0) return 0;

    qsort(ylist, cnt, sizeof(YPt), cmp_yp);

    int m = 1;
    for (int i = 1; i < cnt; ++i)
        if (fabs(ylist[i].y - ylist[i-1].y) > 1e-12)
            ylist[m++].y = ylist[i].y;
    cnt = m;

    int k = 0;
    for(int i = 0; i < cnt; ++i)
        if (b > ylist[i].y + 1e-12) ++k;

    for (int i = 0; i < cnt; ++i)
        if (fabs(b - ylist[i].y) < 1e-9) return 0;

    return (k & 1) ? 1 : -1;
}

int profile_vtp_filter(const double * pos, void * ctx)
{
    VtpFilterCtx * _c = (VtpFilterCtx *) ctx;
    SALEc2dProfile * _prof = _c->p;
    double x = pos[0], y = pos[1];
    double tol = _c->t;

    float * px=NULL, *py=NULL;
    int np = 0;

    if(_prof->seg == NULL)
    {
        px = _prof->pos;
        py = _prof->pos2;
        np = _prof->n;
    }
    else
    {
        // find the longest segment
        int i = 0;
        for(int k=1;k<=_prof->seg[0];++k)
        {
            if(_prof->seg[k] > np)
            {
                np = _prof->seg[k];
                px = _prof->pos  + i;
                py = _prof->pos2 + i;
            }
            i += _prof->seg[k];
        }
    }

    // calculate the nearest distance
    double min_dist;
    int side = side_of_polyline(np, px, py, x, y, &min_dist);
    return (-1 == side && min_dist<= tol)?1:0;
}

int load_vtp_grid(VtpFile * vfp, const char * fname)
{
    FILE * fp = fopen(fname,"r");
    if(NULL == fp)
    {
        fprintf(stderr,"cannot open %s\n",fname);
        return 0;
    }
    int nx=0,ny=0;
    fscanf(fp,"%d %d\n",&nx,&ny);
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

    int i_gx = -1, i_gy = -1;
    for(int k=0;k<vfp->PointNoF;++k)
    {
        if(strcasecmp("gx",vfp->PointField[k].Name) == 0)
        {
            i_gx = k;
        }

        if(strcasecmp("gy",vfp->PointField[k].Name) == 0)
        {
            i_gy = k;
        }
    }

    assert(i_gx!=-1 && i_gy!=-1);

    for(int k=0;k<vfp->NoP;++k)
    {
        int gx = (int) roundf(vfp->PointField[i_gx].Data[k]);
        int gy = (int) roundf(vfp->PointField[i_gy].Data[k]);

        assert(0<=gx && gx<nx);
        assert(0<=gy && gy<ny);

        vfp->PointField[i_gx].Data[k] = 0.5f*(x[gx] + x[gx+1]);
        vfp->PointField[i_gy].Data[k] = 0.5f*(y[gy] + y[gy+1]);
    }

    free(x);
    free(y);
    return vfp->NoP;
}
