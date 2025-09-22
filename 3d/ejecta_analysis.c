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
    float * tr_dis = malloc(sizeof(float)*tr_len);
    float * tr_r   = malloc(sizeof(float)*tr_len);

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

        tr_r[k] = _cur->r;
        tr_dis[k] = _cur->dis;
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
    vtk_dataarray_vec_f(fp,"dis",vtp_data_format,tr_dis,tr_len,1);
    vtk_dataarray_vec_f(fp,"r",vtp_data_format,tr_r,tr_len,1);

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
    free(tr_dis);
    free(tr_r);
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

    // for(int j=0; j<_cs->nproc_surf; ++j)
    // {
    //     citcoms_sphere_dump * _csd = _cs->cap + j;
    //     for(int i=0; i<_csd->nno; ++i)
    //     {
    //         if(_csd->pdata[i * _csd->noc + _csd->noc - 1] >= 1)
    //             _csd->pdata[i * _csd->noc + 0] = 0;
    //         /// adjust coordinates
    //         _csd->pos[3*i + 0] = R*_csd->pos[3*i + 0];
    //         _csd->pos[3*i + 1] = R*_csd->pos[3*i + 1];
    //         _csd->pos[3*i + 2] = R*_csd->pos[3*i + 2] - R;
    //     }
    // }
}

void calculate_ejecta_thickness_pg(citcoms_sphere * _cs, ejecta_collect * _ec, double R)
{
    const double v0 = _ec->v0;
    const double bandwidth = 20.0e3;
    const double s0 = bandwidth*bandwidth*M_PI;
    fprintf(stdout,"v0 is %f km3; %f km\n", v0*1e-9, v0/s0*1e-3);

    #pragma omp parallel for num_threads(4) shared(_cs,_ec,R) default(none)
    for(int j=0; j<_cs->nproc_surf; ++j)
    {
        citcoms_sphere_dump * _csd = _cs->cap + j;
        int nox = _csd->nox, noy=_csd->noy, noz=_csd->noz, noc=_csd->noc;

        int nodeS[4] = {
                citcoms_offset(1,1,1,nox,noy,noz) - 1,
                citcoms_offset(nox,1,1,nox,noy,noz) - 1,
                citcoms_offset(1,noy,1,nox,noy,noz) - 1,
                citcoms_offset(nox,noy,1,nox,noy,noz) - 1
        };

        double nodePos[4][3] = {0};
        for(int s=0;s<4;++s)
        {
            for(int d=0;d<3;++d)
                nodePos[s][d] = _csd->pos[3*nodeS[s] + d];
        }

        for(int k=0;k<_ec->len;++k)
        {
            ejecta_t *_cur = _ec->data + k;
            double cpos[3] = {_cur->pos[0], _cur->pos[1], _cur->pos[2] + R};

            double cpos_len = VecLen(cpos, 3);
            double npos[3] = {cpos[0]/cpos_len, cpos[1]/cpos_len, cpos[2]/cpos_len};
            double l0 = VecLen(_cs->cap_info[j+1].R,3), l1=VecDot(npos, _cs->cap_info[j+1].R,3);
            if(l1 < 0.9*_cs->cap_info[j+1].dx)
                continue;

            double lxt = (nox-1)*solve_local(npos, _cs->cap_info[j+1].P[2], _cs->cap_info[j+1].P[1], _cs->cap_info[j+1].P[0]) + 1;
            double lyt = (noy-1)*solve_local(npos, _cs->cap_info[j+1].Q[2], _cs->cap_info[j+1].Q[1], _cs->cap_info[j+1].Q[0]) + 1;
            if(!(lxt <= nox && lxt >= 0 && lyt <= noy && lyt >= 0))
            {
                continue;
            }

            // int Ixt = round(lxt), Iyt = round(lyt);
            // int eid[4] = {0, Ixt+1, Iyt+1, 1};
            // eid[1] = (lxt - Ixt >= 0) ? Ixt : Ixt - 1;
            // eid[2] = (lyt - Iyt >= 0) ? Iyt : Iyt - 1;

            int eid[4] = {0, floor(lxt), floor(lyt), 1};
            if(eid[1]*(eid[1] - nox)*eid[2]*(eid[2] - noy) == 0)
                continue;
            const int n2 = citcoms_offset(eid[1],eid[2],1,nox-1,noy-1,1)-1;

            _csd->data[n2*noc + 0] += v0/(_csd->area[n2]*R*R);
        }

    }
}

void find_ejecta_contour(citcoms_sphere * _cs, double * _t0, int marker_id)
{
    int max_loc_i, max_loc_j, max_loc_n2, max_loc_p;
    double max_thickness = 0., max_pos[3] = {0.};
    for(int j=0; j<_cs->nproc_surf; ++j)
    {
        citcoms_sphere_dump * _csd = _cs->cap + j;
        int nox = _csd->nox, noy=_csd->noy, noz=_csd->noz, noc=_csd->noc;
        for(int ix=0;ix<nox-1;++ix)
            for(int jy=0;jy<noy-1;++jy)
            {
                int kz = 0;
                int eid[4] = {0, ix+1, jy+1, kz+1};
                const int n2 = citcoms_offset(ix+1,jy+1,1,nox-1,noy-1,1)-1;
                if(_csd->data[n2*noc + 0] > max_thickness)
                {
                    max_thickness = _csd->data[n2*noc + 0];
                    max_loc_i = ix;
                    max_loc_j = jy;
                    max_loc_n2 = n2;
                    max_loc_p = j;

                    int n2ien[4] = {citcoms_offset(eid[1]+1,eid[2]+1,eid[3]+0,nox, noy, 1) - 1,
                                    citcoms_offset(eid[1]+0,eid[2]+1,eid[3]+0,nox, noy, 1) - 1,
                                    citcoms_offset(eid[1]+0,eid[2]+0,eid[3]+0,nox, noy, 1) - 1,
                                    citcoms_offset(eid[1]+1,eid[2]+0,eid[3]+0,nox, noy, 1) - 1};

                    VecZero(max_pos,3);
                    for(int s=0;s<4;++s)
                    {
                        for(int d=0;d<3;++d)
                        {
                            max_pos[d] += _csd->pos[3*n2ien[s] + d] * 0.25;
                        }
                    }
                    VecNormalize(max_pos,3);
                }
            }
    }

    int npts = 128;
    double * con_pts = malloc(sizeof(double)*npts*3);
    double * con_t =  malloc(sizeof(double)*npts*3);

    double con_I[3] = {max_pos[2], 0, -max_pos[0]};
    double con_J[3] = {0};
    VecNormalize(con_I, 3);
    VecCross(con_J, con_I, max_pos, 3);
    VecNormalize(con_J, 3);
    double dx0 = sqrt(2.0*_cs->cap[max_loc_p].area[max_loc_n2]);
    for(int k=0;k<npts;++k)
    {
        double k_direction[3] = {0};
        VecLinear(k_direction, con_I, cos(k*M_PI*2/npts), con_J, sin(k*M_PI*2/npts), 3);
        VecNormalize(k_direction, 3);
        VecLinear(con_pts + 3*k, max_pos, 1.0, k_direction, dx0, 3);
    }


    free(con_pts);
    free(con_t);
}

void marker_in_cap(citcoms_sphere * _cs, double _t0, int m)
{
    int update_pts = 1;
    while(update_pts >= 1)
    {
        update_pts = 0;
        citcoms_sphere_dump *_csd = _cs->cap + m;
        int nox = _csd->nox, noy = _csd->noy, noc = _csd->noc;
        for(int ix = 0; ix < nox - 1; ++ix)
        {
            for(int jy = 0; jy < noy - 1; ++jy)
            {
                int eid[4] = {0, ix + 1, jy + 1, 1};
                int n2ien[4] = {citcoms_offset(eid[1] + 1, eid[2] + 1, eid[3] + 0, nox, noy, 1) - 1,
                                citcoms_offset(eid[1] + 0, eid[2] + 1, eid[3] + 0, nox, noy, 1) - 1,
                                citcoms_offset(eid[1] + 0, eid[2] + 0, eid[3] + 0, nox, noy, 1) - 1,
                                citcoms_offset(eid[1] + 1, eid[2] + 0, eid[3] + 0, nox, noy, 1) - 1};
                int n2 = citcoms_offset(ix + 1, jy + 1, 1, nox - 1, noy - 1, 1) - 1;
                if(_csd->data[n2 * noc + 0] < _t0)
                    continue;
                if(_csd->pmarker[n2ien[0]] + _csd->pmarker[n2ien[1]] + _csd->pmarker[n2ien[2]] +
                   _csd->pmarker[n2ien[3]] < 1)
                    continue;

                for(int s = 0; s < 4; ++s)
                {
                    if(_csd->pmarker[n2ien[s]] == 1)
                        continue;
                    _csd->pmarker[n2ien[s]] = 1;
                    update_pts++;
                }
            }
        }
    }

    int * update_ne = malloc(sizeof(int)*_cs->nproc_surf);
    for(int k=0;k<_cs->nproc_surf;++k)
    {
        update_ne[k] = 0;
    }

    citcoms_sphere_dump *_csd = _cs->cap + m;
    int nox = _csd->nox, noy = _csd->noy, noc = _csd->noc;
    for(int ix = 0; ix < nox; ++ix)
    {
        for(int jy = 0; jy < noy; ++jy)
        {
            if(ix*(ix-nox+1)*jy*(jy-noy+1) != 0)
                continue;
            int n2 = citcoms_offset(ix+1,jy+1,1,nox,noy,1) - 1;
            if(_csd->pmarker[n2] < 1.0f)
                continue;
            double pos[3] = {_csd->pos[n2*3+0], _csd->pos[n2*3+1], _csd->pos[n2*3+2]};
            VecNormalize(pos,3);
            for(int j=0;j<_cs->nproc_surf;++j)
            {
                if(j == m)
                    continue;
                double l1=VecDot(pos, _cs->cap_info[j+1].R,3);
                if(l1 < 0.9*_cs->cap_info[j+1].dx)
                    continue;
                double lxt = (nox-1)*solve_local(pos, _cs->cap_info[j+1].P[2], _cs->cap_info[j+1].P[1], _cs->cap_info[j+1].P[0]);
                double lyt = (noy-1)*solve_local(pos, _cs->cap_info[j+1].Q[2], _cs->cap_info[j+1].Q[1], _cs->cap_info[j+1].Q[0]);
                int Ixt = (int)round(lxt), Iyt= (int)round(lyt);
                if(Ixt>=0 && Ixt<= nox -1 && Iyt>=0 && Iyt<=nox-1)
                {
                    int n2v = citcoms_offset(Ixt+1,Iyt+1,1,nox,noy,1)-1;
                    if(_cs->cap[j].pmarker[n2v] < 1.0f)
                    {
                        _cs->cap[j].pmarker[n2v] = 1.0f;
                        // fprintf(stdout,"<%d,%d,%d> = <%d,%d,%d>\n", m, ix, jy, j, Ixt, Iyt);
                        update_ne[j] += 1;
                    }
                }
            }
        }
    }

    for(int j=0;j<_cs->nproc_surf;++j)
    {
        if(update_ne[j] == 0)
            continue;
        marker_in_cap(_cs, _t0, j);
    }

    free(update_ne);
}

void remove_center_ejecta_0(citcoms_sphere * _cs, double _t0)
{
    int max_loc_i, max_loc_j, max_loc_n2, max_loc_p;
    double max_thickness = 0., max_pos[3] = {0.};
    for(int j=0; j<_cs->nproc_surf; ++j)
    {
        citcoms_sphere_dump * _csd = _cs->cap + j;
        int nox = _csd->nox, noy=_csd->noy, noz=_csd->noz, noc=_csd->noc;
        for(int ix=0;ix<nox-1;++ix)
            for(int jy=0;jy<noy-1;++jy)
            {
                int kz = 0;
                int eid[4] = {0, ix+1, jy+1, kz+1};
                const int n2 = citcoms_offset(ix+1,jy+1,1,nox-1,noy-1,1)-1;
                if(_csd->data[n2*noc + 0] > max_thickness)
                {
                    max_thickness = _csd->data[n2*noc + 0];
                    max_loc_i = ix;
                    max_loc_j = jy;
                    max_loc_n2 = n2;
                    max_loc_p = j;

                    int n2ien[4] = {citcoms_offset(eid[1]+1,eid[2]+1,eid[3]+0,nox, noy, 1) - 1,
                                    citcoms_offset(eid[1]+0,eid[2]+1,eid[3]+0,nox, noy, 1) - 1,
                                    citcoms_offset(eid[1]+0,eid[2]+0,eid[3]+0,nox, noy, 1) - 1,
                                    citcoms_offset(eid[1]+1,eid[2]+0,eid[3]+0,nox, noy, 1) - 1};

                    VecZero(max_pos,3);
                    for(int s=0;s<4;++s)
                    {
                        for(int d=0;d<3;++d)
                        {
                            max_pos[d] += _csd->pos[3*n2ien[s] + d] * 0.25;
                        }
                    }
                    VecNormalize(max_pos,3);
                }
            }
    }

    {
        // marker the initial center
        int eid[4] = {0, max_loc_i+1, max_loc_j+1, 1};
        int nox = _cs->cap[max_loc_p].nox, noy = _cs->cap[max_loc_p].noy;
        int n2ien[4] = {citcoms_offset(eid[1]+1,eid[2]+1,eid[3]+0,nox, noy, 1) - 1,
                        citcoms_offset(eid[1]+0,eid[2]+1,eid[3]+0,nox, noy, 1) - 1,
                        citcoms_offset(eid[1]+0,eid[2]+0,eid[3]+0,nox, noy, 1) - 1,
                        citcoms_offset(eid[1]+1,eid[2]+0,eid[3]+0,nox, noy, 1) - 1};
        _cs->cap[max_loc_p].pmarker[n2ien[0]] = 1;
        _cs->cap[max_loc_p].pmarker[n2ien[1]] = 1;
        _cs->cap[max_loc_p].pmarker[n2ien[2]] = 1;
        _cs->cap[max_loc_p].pmarker[n2ien[3]] = 1;
    }
    marker_in_cap(_cs, _t0, max_loc_p);

    for(int j=0; j<_cs->nproc_surf; ++j)
    {
        citcoms_sphere_dump *_csd = _cs->cap + j;
        int nox = _csd->nox, noy = _csd->noy, noc = _csd->noc;
        for(int ix = 0; ix < nox - 1; ++ix)
        {
            for(int jy = 0; jy < noy - 1; ++jy)
            {
                int eid[4] = {0, ix + 1, jy + 1, 1};
                int n2ien[4] = {citcoms_offset(eid[1] + 1, eid[2] + 1, eid[3] + 0, nox, noy, 1) - 1,
                                citcoms_offset(eid[1] + 0, eid[2] + 1, eid[3] + 0, nox, noy, 1) - 1,
                                citcoms_offset(eid[1] + 0, eid[2] + 0, eid[3] + 0, nox, noy, 1) - 1,
                                citcoms_offset(eid[1] + 1, eid[2] + 0, eid[3] + 0, nox, noy, 1) - 1};
                int n2 = citcoms_offset(ix + 1, jy + 1, 1, nox - 1, noy - 1, 1) - 1;

                if(_csd->pmarker[n2ien[0]] + _csd->pmarker[n2ien[1]] + _csd->pmarker[n2ien[2]] +
                   _csd->pmarker[n2ien[3]] < 1)
                    continue;

                _csd->data[n2*noc + 0] = 0.f;
            }
        }
    }

}

void adjust_coordinate(citcoms_sphere * _cs, double R)
{
    for(int j=0; j<_cs->nproc_surf; ++j)
    {
        citcoms_sphere_dump * _csd = _cs->cap + j;
        for(int i=0; i<_csd->nno; ++i)
        {
            int noc=_csd->noc;
            if(_csd->pdata[i * noc + noc - 1] >= 1)
                _csd->pdata[i * noc + 0] = 0;
            /// adjust coordinates
            _csd->pos[3*i + 0] = R*_csd->pos[3*i + 0];
            _csd->pos[3*i + 1] = R*_csd->pos[3*i + 1];
            _csd->pos[3*i + 2] = R*_csd->pos[3*i + 2] - R;

            // _csd->data[i*noc + 0] /= _csd->area[i]*R*R;
        }
    }
}

int above_threshold(citcoms_sphere * _cs, double * pars)
{
    int m  = (int) pars[0];
    int n2 = (int) pars[1];
    double t0 = pars[2];
    int fn = (int) pars[3];
    citcoms_sphere_dump * _csd = _cs->cap + m;
    int noc = _csd->noc;

    if(fn == -1)
    {
        return _csd->data[n2*noc + 0] < t0;
    }
    else if(fn == 1)
    {
        return _csd->data[n2*noc + 0] > t0;
    }
    else
    {
        fprintf(stdout,"%s:undefined option fn=%d\n",__func__,fn);
        exit(0);
        return 0;
    }
}

void marker_in_fn(citcoms_sphere * _cs,int m, int (*fn)(citcoms_sphere *, double *), double *pars)
{
    citcoms_sphere_dump *_csd = _cs->cap + m;
    int nox = _csd->nox, noy = _csd->noy, noc = _csd->noc;
    int update_pts = 1;
    while(update_pts >= 1)
    {
        update_pts = 0;
        for(int ix = 0; ix < nox - 1; ++ix)
        {
            for(int jy = 0; jy < noy - 1; ++jy)
            {
                int eid[4] = {0, ix + 1, jy + 1, 1};
                int n2ien[4] = {citcoms_offset(eid[1] + 1, eid[2] + 1, eid[3] + 0, nox, noy, 1) - 1,
                                citcoms_offset(eid[1] + 0, eid[2] + 1, eid[3] + 0, nox, noy, 1) - 1,
                                citcoms_offset(eid[1] + 0, eid[2] + 0, eid[3] + 0, nox, noy, 1) - 1,
                                citcoms_offset(eid[1] + 1, eid[2] + 0, eid[3] + 0, nox, noy, 1) - 1};
                int n2 = citcoms_offset(ix + 1, jy + 1, 1, nox - 1, noy - 1, 1) - 1;
                if(0 == fn(_cs,(double []){m, n2, pars[2], pars[3]}))
                    continue;
                if(_csd->pmarker[n2ien[0]] + _csd->pmarker[n2ien[1]] + _csd->pmarker[n2ien[2]] +
                   _csd->pmarker[n2ien[3]] < 1)
                    continue;

                for(int s = 0; s < 4; ++s)
                {
                    if(_csd->pmarker[n2ien[s]] == 1)
                        continue;
                    _csd->pmarker[n2ien[s]] = 1;
                    update_pts++;
                }
            }
        }
    }

    int * update_ne = malloc(sizeof(int)*_cs->nproc_surf);
    for(int k=0;k<_cs->nproc_surf;++k)
    {
        update_ne[k] = 0;
    }


    for(int ix = 0; ix < nox; ++ix)
    {
        for(int jy = 0; jy < noy; ++jy)
        {
            if(ix*(ix-nox+1)*jy*(jy-noy+1) != 0)
                continue;
            int n2 = citcoms_offset(ix+1,jy+1,1,nox,noy,1) - 1;
            if(_csd->pmarker[n2] < 1.0f)
                continue;
            double pos[3] = {_csd->pos[n2*3+0], _csd->pos[n2*3+1], _csd->pos[n2*3+2]};
            VecNormalize(pos,3);
            for(int j=0;j<_cs->nproc_surf;++j)
            {
                if(j == m)
                    continue;
                double l1=VecDot(pos, _cs->cap_info[j+1].R,3);
                if(l1 < 0.9*_cs->cap_info[j+1].dx)
                    continue;
                double lxt = (nox-1)*solve_local(pos, _cs->cap_info[j+1].P[2], _cs->cap_info[j+1].P[1], _cs->cap_info[j+1].P[0]);
                double lyt = (noy-1)*solve_local(pos, _cs->cap_info[j+1].Q[2], _cs->cap_info[j+1].Q[1], _cs->cap_info[j+1].Q[0]);
                int Ixt = (int)round(lxt), Iyt= (int)round(lyt);
                if(Ixt>=0 && Ixt<= nox -1 && Iyt>=0 && Iyt<=nox-1)
                {
                    int n2v = citcoms_offset(Ixt+1,Iyt+1,1,nox,noy,1)-1;
                    if(_cs->cap[j].pmarker[n2v] < 1.0f)
                    {
                        _cs->cap[j].pmarker[n2v] = 1.0f;
                        update_ne[j] += 1;
                    }
                }
            }
        }
    }

    for(int j=0;j<_cs->nproc_surf;++j)
    {
        if(update_ne[j] == 0)
            continue;
        marker_in_fn(_cs, j, fn, pars);
    }

    free(update_ne);
}


int marker_from_fn(citcoms_sphere * _cs,int start_cap, int (*fn)(citcoms_sphere *, double *), double *pars)
{

    marker_in_fn(_cs,start_cap, fn, pars);

    for(int j=0; j<_cs->nproc_surf; ++j)
    {
        citcoms_sphere_dump *_csd = _cs->cap + j;
        int nox = _csd->nox, noy = _csd->noy, noc = _csd->noc;
        for(int ix = 0; ix < nox - 1; ++ix)
        {
            for(int jy = 0; jy < noy - 1; ++jy)
            {
                int eid[4] = {0, ix + 1, jy + 1, 1};
                int n2ien[4] = {citcoms_offset(eid[1] + 1, eid[2] + 1, eid[3] + 0, nox, noy, 1) - 1,
                                citcoms_offset(eid[1] + 0, eid[2] + 1, eid[3] + 0, nox, noy, 1) - 1,
                                citcoms_offset(eid[1] + 0, eid[2] + 0, eid[3] + 0, nox, noy, 1) - 1,
                                citcoms_offset(eid[1] + 1, eid[2] + 0, eid[3] + 0, nox, noy, 1) - 1};
                int n2 = citcoms_offset(ix + 1, jy + 1, 1, nox - 1, noy - 1, 1) - 1;

                if(_csd->pmarker[n2ien[0]] + _csd->pmarker[n2ien[1]] + _csd->pmarker[n2ien[2]] +
                   _csd->pmarker[n2ien[3]] < 1)
                    continue;
                _csd->data[n2*noc + 0] = 0.f;
            }
        }
    }
    return 0;
}

void remove_center_ejecta(citcoms_sphere * _cs, double _t0, double * crater)
{
    int max_loc_i, max_loc_j, max_loc_n2, max_loc_p;
    double max_thickness = 0., max_pos[3] = {0.};
    for(int j=0; j<_cs->nproc_surf; ++j)
    {
        citcoms_sphere_dump * _csd = _cs->cap + j;
        int nox = _csd->nox, noy=_csd->noy, noz=_csd->noz, noc=_csd->noc;
        for(int ix=0;ix<nox-1;++ix)
            for(int jy=0;jy<noy-1;++jy)
            {
                int kz = 0;
                int eid[4] = {0, ix+1, jy+1, kz+1};
                const int n2 = citcoms_offset(ix+1,jy+1,1,nox-1,noy-1,1)-1;
                if(_csd->data[n2*noc + 0] > max_thickness)
                {
                    max_thickness = _csd->data[n2*noc + 0];
                    max_loc_i = ix;
                    max_loc_j = jy;
                    max_loc_n2 = n2;
                    max_loc_p = j;

                    int n2ien[4] = {citcoms_offset(eid[1]+1,eid[2]+1,eid[3]+0,nox, noy, 1) - 1,
                                    citcoms_offset(eid[1]+0,eid[2]+1,eid[3]+0,nox, noy, 1) - 1,
                                    citcoms_offset(eid[1]+0,eid[2]+0,eid[3]+0,nox, noy, 1) - 1,
                                    citcoms_offset(eid[1]+1,eid[2]+0,eid[3]+0,nox, noy, 1) - 1};

                    VecZero(max_pos,3);
                    for(int s=0;s<4;++s)
                    {
                        for(int d=0;d<3;++d)
                        {
                            max_pos[d] += _csd->pos[3*n2ien[s] + d] * 0.25;
                        }
                    }
                    VecNormalize(max_pos,3);
                }
            }
    }

    {
        // marker the initial center
        int eid[4] = {0, max_loc_i+1, max_loc_j+1, 1};
        int nox = _cs->cap[max_loc_p].nox, noy = _cs->cap[max_loc_p].noy;
        int n2ien[4] = {citcoms_offset(eid[1]+1,eid[2]+1,eid[3]+0,nox, noy, 1) - 1,
                        citcoms_offset(eid[1]+0,eid[2]+1,eid[3]+0,nox, noy, 1) - 1,
                        citcoms_offset(eid[1]+0,eid[2]+0,eid[3]+0,nox, noy, 1) - 1,
                        citcoms_offset(eid[1]+1,eid[2]+0,eid[3]+0,nox, noy, 1) - 1};
        _cs->cap[max_loc_p].pmarker[n2ien[0]] = 1;
        _cs->cap[max_loc_p].pmarker[n2ien[1]] = 1;
        _cs->cap[max_loc_p].pmarker[n2ien[2]] = 1;
        _cs->cap[max_loc_p].pmarker[n2ien[3]] = 1;
    }

    marker_from_fn(_cs, max_loc_p, above_threshold, (double []){0,0, _t0,1});

    for(int j=0; j<_cs->nproc_surf; ++j)
    {
        citcoms_sphere_dump *_csd = _cs->cap + j;
        int nox = _csd->nox, noy = _csd->noy, noc = _csd->noc;
        for(int ix = 0; ix < nox - 1; ++ix)
        {
            for(int jy = 0; jy < noy - 1; ++jy)
            {
                int eid[4] = {0, ix + 1, jy + 1, 1};
                int n2ien[4] = {citcoms_offset(eid[1] + 1, eid[2] + 1, eid[3] + 0, nox, noy, 1) - 1,
                                citcoms_offset(eid[1] + 0, eid[2] + 1, eid[3] + 0, nox, noy, 1) - 1,
                                citcoms_offset(eid[1] + 0, eid[2] + 0, eid[3] + 0, nox, noy, 1) - 1,
                                citcoms_offset(eid[1] + 1, eid[2] + 0, eid[3] + 0, nox, noy, 1) - 1};
                int n2 = citcoms_offset(ix + 1, jy + 1, 1, nox - 1, noy - 1, 1) - 1;

                if(_csd->pmarker[n2ien[0]] + _csd->pmarker[n2ien[1]] + _csd->pmarker[n2ien[2]] +
                   _csd->pmarker[n2ien[3]] < 1)
                    continue;
                _csd->data[n2*noc + 0] = 0.f;
            }
        }
    }

    /// marker from center, calculate crater area/radius
    {
        for(int j=0; j<_cs->nproc_surf; ++j)
        {
            citcoms_sphere_dump *_csd = _cs->cap + j;
            int nno = _csd->nno, noc = _csd->noc;
            memset(_csd->pmarker, 0, sizeof(_csd->pmarker[0])*nno);
        }
        // marker the initial center
        int eid[4] = {0, max_loc_i+1, max_loc_j+1, 1};
        int nox = _cs->cap[max_loc_p].nox, noy = _cs->cap[max_loc_p].noy;
        int n2ien[4] = {citcoms_offset(eid[1]+1,eid[2]+1,eid[3]+0,nox, noy, 1) - 1,
                        citcoms_offset(eid[1]+0,eid[2]+1,eid[3]+0,nox, noy, 1) - 1,
                        citcoms_offset(eid[1]+0,eid[2]+0,eid[3]+0,nox, noy, 1) - 1,
                        citcoms_offset(eid[1]+1,eid[2]+0,eid[3]+0,nox, noy, 1) - 1};
        _cs->cap[max_loc_p].pmarker[n2ien[0]] = 1;
        _cs->cap[max_loc_p].pmarker[n2ien[1]] = 1;
        _cs->cap[max_loc_p].pmarker[n2ien[2]] = 1;
        _cs->cap[max_loc_p].pmarker[n2ien[3]] = 1;
    }

    marker_from_fn(_cs, max_loc_p, above_threshold, (double []){0,0, _t0,-1});

    /// using ejecta blank area to derive (cenetr, radius) of crater
    float center[3] = {0,0,0};
    int ncell = 0;
    double sum_area = 0;
    for(int j=0; j<_cs->nproc_surf; ++j)
    {
        citcoms_sphere_dump *_csd = _cs->cap + j;
        int nox = _csd->nox, noy = _csd->noy, noc = _csd->noc;
        for(int ix = 0; ix < nox - 1; ++ix)
        {
            for(int jy = 0; jy < noy - 1; ++jy)
            {
                int eid[4] = {0, ix + 1, jy + 1, 1};
                int n2ien[4] = {citcoms_offset(eid[1] + 1, eid[2] + 1, eid[3] + 0, nox, noy, 1) - 1,
                                citcoms_offset(eid[1] + 0, eid[2] + 1, eid[3] + 0, nox, noy, 1) - 1,
                                citcoms_offset(eid[1] + 0, eid[2] + 0, eid[3] + 0, nox, noy, 1) - 1,
                                citcoms_offset(eid[1] + 1, eid[2] + 0, eid[3] + 0, nox, noy, 1) - 1};
                int n2 = citcoms_offset(ix + 1, jy + 1, 1, nox - 1, noy - 1, 1) - 1;

                if(_csd->pmarker[n2ien[0]] + _csd->pmarker[n2ien[1]] + _csd->pmarker[n2ien[2]] +
                   _csd->pmarker[n2ien[3]] < 1)
                    continue;

                float ecenter[3] = {0};
                for(int s=0;s<4;++s)
                {
                    VecAddF(ecenter, _csd->pos + 3*n2ien[s],0.25f,3);
                }
                VecAddF(center, ecenter, 1.0f, 3);
                ncell++;
                sum_area += _csd->area[n2];
            }
        }
    }

    double radius = sqrt(sum_area/M_PI);
    double center2[3];
    VecF2D(center2, center, 3);
    // VecScale(center2, 1.74e6/ncell, 3);
    VecNormalize(center2,3);
    fprintf(stdout,"estimated radius %10.5f (%f,%f,%f)\n", radius*1.74e6, center2[0], center2[1], center2[2]);

    if(crater != NULL)
    {
        crater[0] = radius;
        crater[1] = sum_area;
        crater[2] = center2[0];
        crater[3] = center2[1];
        crater[4] = center2[2];
    }
}

void calculate_land_skew(ejecta_collect * _ec, double * _c, double R)
{
    // _c: results from remove center ejecta
    // radius[0], area[1], center[2,3,4]

    double cpos[3] = {_c[2]*R, _c[3]*R, _c[4]*R};
    double L0 = _c[4]*R/ VecLen(cpos,3);
    L0 = acos(L0) * R;

    for(int k=0;k<_ec->len;++k)
    {
        ejecta_t *_cur = _ec->data + k;
        if(VecLen(_cur->pos,3) < 1e-2)
        {
            _cur->dis = _cur->r = 0;
            continue;
        }
        double epos[3] = {_cur->pos[0], _cur->pos[1], _cur->pos[2] + R}; // change origin to moon center
        double evel[3] = {_cur->vel[0], _cur->vel[1], _cur->vel[2]};
        _cur->dis = R*acos(VecDot(cpos,epos,3)/(VecLen(cpos,3)* VecLen(epos,3)));
        double XxV[3];
        VecCross(XxV, epos, evel, 3);
        VecNormalize(XxV, 3);
        VecNormalize(cpos,3);
        VecNormalize(epos,3);

        double pos_direction[3];
        VecLinear(pos_direction, cpos, -1.0, epos, 1.0, 3);
        double tmp0 = VecDot(pos_direction, epos, 3);
        double tg_direction[3];
        VecLinear(tg_direction, pos_direction, 1.0, epos, -tmp0, 3);
        VecNormalize(tg_direction, 3);

        double tg_vel[3];
        double tmp1 = VecDot(evel, epos, 3);
        VecLinear(tg_vel, evel, 1.0, epos, -tmp1, 3);
        VecNormalize(tg_vel, 3);

        double CosDelta = fabs(VecDot(cpos,XxV,3));
        //_cur->r = (M_PI/2.0 - acos(CosDelta))*R;
        _cur->r = acos(VecDot(tg_vel, tg_direction,3))/M_PI*180.0;
        if(_cur->r > 90.0)
        {
            _cur->r = - 180.0 + _cur->r;
        }

        _cur->theta = L0 - _cur->theta;
    }
}
