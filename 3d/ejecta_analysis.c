//
// Created by li on 8/20/24.
//

#include "ejecta_analysis.h"
#include "omp.h"
int load_ejecta_collect(ejecta_collect * _ec, int step)
{
    int new_ejecta_num = 0;
    for(int k=0;k<_ec->nproc;++k)
    {
        char _tmp_name[4097];
        snprintf(_tmp_name,4096,"%s.proc%d.%04d.ejecta",_ec->prefix,k,step);
        FILE * fp = fopen(_tmp_name,"r");
        if(NULL == fp){
            fprintf(stdout,"cannot open %s\n",_tmp_name);
            continue;
        }
        ejecta_t tmp_e;
        int NeE;
        while(15 == fscanf(fp,"%d, %d, %d, %d, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf",
                         &(tmp_e.id),&(tmp_e.rank), &(tmp_e.matid), &(NeE),&(tmp_e.t),
                         tmp_e.pos,tmp_e.pos+1,tmp_e.pos+2,
                         &(tmp_e.maxpre),&(tmp_e.maxtem),
                         tmp_e.vel, tmp_e.vel+1, tmp_e.vel+2,
                         &(tmp_e.pre),&(tmp_e.tem)
                         ))
        // while(14 == fscanf(fp,"%d, %d, %d,  %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf",
        //                    &(tmp_e.id),&(tmp_e.rank), &(tmp_e.matid), &(tmp_e.t),
        //                    tmp_e.pos,tmp_e.pos+1,tmp_e.pos+2,
        //                    &(tmp_e.maxpre),&(tmp_e.maxtem),
        //                    tmp_e.vel, tmp_e.vel+1, tmp_e.vel+2,
        //                    &(tmp_e.pre),&(tmp_e.tem)
        // ))
        {
            ejecta_collect_push(_ec,&tmp_e);
            new_ejecta_num++;
        }

        fclose(fp);
    }
    return new_ejecta_num;
}

int load_ejecta_collect_single_file(ejecta_collect * _ec, const char * _tmp_name)
{
    int new_ejecta_num = 0;

    FILE * fp = fopen(_tmp_name,"r");
    if(NULL == fp){
        fprintf(stdout,"cannot open %s\n",_tmp_name);
        exit(0);
    }
    ejecta_t tmp_e;
    while(13 == fscanf(fp,"%d, %d, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf",
                       &(tmp_e.id),&(tmp_e.rank),&(tmp_e.t),
                       tmp_e.pos,tmp_e.pos+1,tmp_e.pos+2,
                       &(tmp_e.maxpre),&(tmp_e.maxtem),
                       tmp_e.vel, tmp_e.vel+1, tmp_e.vel+2,
                       &(tmp_e.pre),&(tmp_e.tem)
    ))
    {
        ejecta_collect_push(_ec,&tmp_e);
        new_ejecta_num++;
    }
    fclose(fp);

    return new_ejecta_num;
}

int ejecta_collect_push(ejecta_collect * _ec, ejecta_t * _e)
{
    if(_ec->len < 0) return 0;

    if(_ec->len+2 > _ec->len_allocated)
    {
        _ec->data = realloc(_ec->data,2*_ec->len_allocated* sizeof(ejecta_t));
        _ec->len_allocated *= 2;

        if(_ec->data == NULL)
        {
            fprintf(stdout,"error in realloc\n");
            exit(0);
        }
    }

     memcpy(_ec->data + _ec->len ,_e, sizeof(ejecta_t));
    _ec->len++;

    return 1;
}

int ejecta_collect_test_init(ejecta_collect * _ec)
{
    _ec->nproc = 490;
    strcpy(_ec->prefix,"./ejecta/bm");
    _ec->len_allocated = 256*256;
    _ec->data = malloc(sizeof(ejecta_t)*_ec->len_allocated);
    _ec->len = 0;
    return _ec->len_allocated;
}

int ejecta_collect_init(ejecta_collect * _ec, InputFile * ifp)
{
    char SALEcInp[4096];
    GetValueS(ifp,"SALEc.input",SALEcInp,"SALEc.inp");
    InputFile * sifp = OpenInputFile(SALEcInp);
    int npgx = GetValueI(sifp,"processor.npgx","2");
    int npgy = GetValueI(sifp,"processor.npgy","2");
    int npgz = GetValueI(sifp,"processor.npgz","2");
    double dx = GetValueD(sifp,"mesh.dx","-1.0");
    double dy = GetValueD(sifp,"mesh.dy","-1.0");
    double dz = GetValueD(sifp,"mesh.dz","-1.0");
    _ec->v0 = dx*dy*dz;
    _ec->nproc = npgx*npgy*npgz;
    CloseInputFile(sifp);

    GetValueS(ifp,"Ejecta.data",_ec->prefix,"bm");
    _ec->max_step = GetValueIk(ifp,"SALEc.step",1,"0");
    _ec->min_step = GetValueIk(ifp,"SALEc.step",0,"1");
    GetValueS(ifp,"Ejecta.output",_ec->output,"test");
    strcat(_ec->output,".vtp");
    _ec->len_allocated = 1024;
    _ec->data = malloc(sizeof(ejecta_t)*_ec->len_allocated);
    _ec->len = 0;
    return _ec->len_allocated;
}


void ejecta_collect_test_clean(ejecta_collect * _ec)
{
    if(_ec->len_allocated > 0)
    {
        free(_ec->data);
    }
}


double numerical_ejecta_orbit_moon_a(double * x,double * v, double dt)
{
    double Rm = 1.74e6;
    double gs = -1.622;

    double r0L = sqrt(x[0]*x[0] + x[1]*x[1] + (x[2]+Rm)*(x[2]+Rm));
    double gr0 = gs*(Rm/r0L)*(Rm/r0L);

    if(r0L < Rm) return (r0L - Rm);


    double xf[3] = {x[0],x[1],x[2]};
    double r0[3] = {x[0],x[1],x[2]+Rm};
    for(int k=0;k<3;++k)
    {
        x[k] = xf[k] + v[k]*dt + gr0*r0[k]/r0L*dt*dt*0.5;
        v[k] = v[k] + gr0*r0[k]/r0L*dt;
    }
    return sqrt(x[0]*x[0] + x[1]*x[1] + (x[2]+Rm)*(x[2]+Rm)) - Rm;
}


int numerical_ejecta_orbit_moon_b(double * x,double * v, double dt)
{
    double Rm = 1.74e6;
    double v0L = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    double local_dt0 = Rm/v0L/1024.0;
    int local_steps = (int)floor(dt/local_dt0) + 1;
    double local_dt = dt/local_steps;

    for(int k=0;k<local_steps;++k)
    {
        numerical_ejecta_orbit_moon_a(x,v,local_dt);
    }
    return local_steps;
}

void numerical_ejecta_orbit_moon(ejecta_collect * _ec, double dt)
{
    for(int k=0;k<_ec->len;++k)
    {
        ejecta_t * _cur = _ec->data + k;
        numerical_ejecta_orbit_moon_b(_cur->pos,_cur->vel,dt);
    }
}


int ejecta_collect_to_vtp(ejecta_collect * _ec, const char * vtp_name)
{
    if(_ec->len <= 0) return 0;

    const unsigned int tr_len = _ec->len;
    float * tr_pos = malloc(sizeof(float)*tr_len*3);
    float * tr_vel = malloc(sizeof(float)*tr_len*3);
    float * tr_mid  = malloc(sizeof(float)*tr_len);
    float * tr_gid  = malloc(sizeof(float)*tr_len);
    float * tr_mpre = malloc(sizeof(float)*tr_len);
    float * tr_mtem = malloc(sizeof(float)*tr_len);
    float * tr_epre = malloc(sizeof(float)*tr_len);
    float * tr_etem = malloc(sizeof(float)*tr_len);
    float * tr_dump = malloc(sizeof(float)*tr_len);
    float * tr_eden = malloc(sizeof(float)*tr_len);
    float * tr_t    = malloc(sizeof(float)*tr_len);
    float * tr_theta = malloc(sizeof(float)*tr_len);

    for(int k =0;k<_ec->len;++k)
    {
        ejecta_t * _cur = _ec->data + k;
        tr_pos[3*k + 0] = _cur->pos[0];
        tr_pos[3*k + 1] = _cur->pos[1];
        tr_pos[3*k + 2] = _cur->pos[2];

        tr_vel[3*k + 0] = _cur->vel[0];
        tr_vel[3*k + 1] = _cur->vel[1];
        tr_vel[3*k + 2] = _cur->vel[2];

        tr_mid[k] = 1.0*_cur->matid;
        tr_gid[k] = 1.0*_cur->id;
        tr_dump[k] = 1.0*_cur->rank;

        tr_mpre[k] = _cur->maxpre;
        tr_mtem[k] = _cur->maxtem;
        tr_epre[k] = _cur->pre;
        tr_etem[k] = _cur->tem;

        tr_t[k] = _cur->t;
        tr_theta[k] = _cur->theta;

    }


    const char * vtp_data_format = "binary";
    FILE * fp = fopen(vtp_name,"w");
    vtp_file_header(fp,tr_len);
    vtk_point_data_header(fp);
    vtk_dataarray_vec_f(fp,"matid",vtp_data_format,tr_mid,tr_len,1);
    vtk_dataarray_vec_f(fp,"id",vtp_data_format,tr_gid,tr_len,1);
    vtk_dataarray_vec_f(fp,"pre",vtp_data_format,tr_epre,tr_len,1);
    vtk_dataarray_vec_f(fp,"mpre",vtp_data_format,tr_mpre,tr_len,1);
    vtk_dataarray_vec_f(fp,"tem",vtp_data_format,tr_etem,tr_len,1);
    vtk_dataarray_vec_f(fp,"mtem",vtp_data_format,tr_mtem,tr_len,1);
    vtk_dataarray_vec_f(fp,"vel",vtp_data_format,tr_vel,tr_len,3);
    vtk_dataarray_vec_f(fp,"dump",vtp_data_format,tr_dump,tr_len,1);
    vtk_dataarray_vec_f(fp,"den",vtp_data_format,tr_eden,tr_len,1);
    vtk_dataarray_vec_f(fp,"t",vtp_data_format,tr_t,tr_len,1);
    vtk_dataarray_vec_f(fp,"theta",vtp_data_format,tr_theta,tr_len,1);
    vtk_point_data_trailer(fp);
    vtk_point_header(fp);
    vtk_dataarray_vec_f(fp,"coordinate",vtp_data_format,tr_pos,tr_len,3);
    vtk_point_trailer(fp);
    vtp_file_trailer(fp);
    fclose(fp);

    free(tr_pos);
    free(tr_vel);
    free(tr_mid);
    free(tr_gid);
    free(tr_mpre);
    free(tr_mtem);
    free(tr_epre);
    free(tr_etem);
    free(tr_dump);
    free(tr_eden);
    free(tr_t);
    free(tr_theta);

    return 1;
}

int ejecta_collect_to_vtm(ejecta_collect * _ec, const char * _vts_name, int nx)
{
    assert(nx >= 4);
    // create sphere grid

    // calculate sum thickness

    // write to vts/vtm files
    return 1;
}


double approximate_ejecta(double *x, double *v, double R, double g0)
{
    double u = g0*R*R;
    double h[3];
    VecCross(h, x, v, 3);

    double h1 = VecLen(h,3), r1 = VecLen(x, 3), v1 = VecLen(v,3);
    double a = -0.5*u/(-u/r1 + 0.5*v1*v1);
    if(a <= 0)
    {
        // unphysical solution
        return -2.0;
    }

    double e = sqrt(1.0 - h1*h1/(a*u));
    if(e >= 1.0 )
    {
        //  parabola or hyper bola
        return -1.0;
    }

    double p = a*(1.0 - e*e);
    double ct1 = (1.0 - p/R)/e, ct0 = (1.0 - p/r1)/e;
    double Q = acos(ct1) + acos(ct0);

    double xl[3] = {0}; // position of land
    VecNormalize(h, 3);
    VecNormalize(x, 3);
    double HxR[3];
    VecCross(HxR, h, x, 3);
    VecLinear(xl, x, cos(Q), HxR, sin(Q), 3);
    VecAdd(xl, h, (1.0 - cos(Q))* VecDot(h,x,3), 3);
    VecScale(xl, R, 3);

    double vl[3] = {0};// velocity of land
    double ne[3] = {0}; // eccentricity vector
    double HxV[3] = {0};
    VecCross(HxV, h, v, 3);
    VecLinear(ne, x, -1.0, HxV, -h1/u, 3);

    double HxE[3] = {0};
    VecCross(HxE, h, ne, 3);
    VecLinear(vl, HxE, u/h1, HxR, u/h1, 3);

    // copy result into x,v
    VecCopy(x, xl, 3);
    VecCopy(v, vl, 3);

    // normalize x
    VecNormalize(x, 3);
    VecScale(x, R, 3);
    return 1;
}

void analytical_ejecta_orbit_moon(ejecta_collect * _ec, double R, double g0)
{
    for(int k=0;k<_ec->len;++k)
    {
        ejecta_t * _cur = _ec->data + k;

        double _pos[3] = {_cur->pos[0], _cur->pos[1], _cur->pos[2] + R};
        double _vel[3] = {_cur->vel[0], _cur->vel[1], _cur->vel[2]};
        double _t = approximate_ejecta(_pos, _vel, R, g0);
        _cur->a = _t;
        if(_t > 0)
        {
            VecCopy(_cur->pos, _pos, 3);
            VecCopy(_cur->vel, _vel, 3);

            double XxV[3] = {0};
            VecCross(XxV,_pos, _vel, 3);
            if(XxV[2] < 0.0)
            {
                _cur->theta = atan2(-XxV[2], -XxV[0]) * R;
            }
            else
            {
                _cur->theta = atan2(XxV[2], XxV[0]) * R;
            }
            _cur->pos[2] -= R;
        }
        else
        {
            VecZero(_cur->pos, 3);
            VecZero(_cur->vel, 3);
            _cur->vel[2] = 1.0;
            _cur->theta = 0.;
        }
    }
}

void calculate_ejecta_thickness(citcoms_sphere * _cs, ejecta_collect * _ec, double R)
{
    const double v0 = _ec->v0;
    const double bandwidth = 20.0e3;
    const double s0 = bandwidth*bandwidth*M_PI;
    fprintf(stdout,"v0 is %f km3; %f km\n", v0*1e-9, v0/s0*1e-3);
    #pragma omp parallel for num_threads(12) shared(_cs,_ec,R) default(none)
    for(int j=0; j<_cs->nproc_surf; ++j)
    {
        citcoms_sphere_dump * _csd = _cs->cap + j;
        for(int i=0; i<_csd->nno; ++i)
        {
            int eid[4] = {0,0,0,0};
            citcoms_eid(i+1,eid,_csd->nox,_csd->noy,_csd->noz);
            for(int k=0;k<_ec->len;++k)
            {
                ejecta_t *_cur = _ec->data + k;
                double cpos[3] = {_cur->pos[0], _cur->pos[1], _cur->pos[2] + R};
                // if(_cur->a < 0)
                //     continue;
                double ipos[3] = {R * _csd->pos[3 * i + 0], R * _csd->pos[3 * i + 1], R * _csd->pos[3 * i + 2]};
                double distance_ki = VecDis(ipos, cpos, 3);
                if(distance_ki <= bandwidth)
                {
                    _csd->pdata[i * _csd->noc + 0] += (float) v0 / s0;
                    _csd->pdata[i * _csd->noc + 1] += 1.0;
                }
            }
        }
    }

    for(int j=0; j<_cs->nproc_surf; ++j)
    {
        citcoms_sphere_dump * _csd = _cs->cap + j;
        for(int i=0; i<_csd->nno; ++i)
        {
            if(_csd->pdata[i * _csd->noc + _csd->noc - 1] >= 1)
                _csd->pdata[i * _csd->noc + 0] = 0;
            /// adjust coordinates
            _csd->pos[3*i + 0] = R*_csd->pos[3*i + 0];
            _csd->pos[3*i + 1] = R*_csd->pos[3*i + 1];
            _csd->pos[3*i + 2] = R*_csd->pos[3*i + 2] - R;
        }
    }
}

