//
// Created by huachengli on 10/17/24.
//

#include "spa_detect.h"


double fr0(double *pos, double *ctx)
{
    /// check weather ctx0 <= |pos| <= ctx1
    double r0 = VecLen(pos,3);
    double r1 = ctx[0]/ctx[2];
    double r2 = ctx[1]/ctx[2];

    if(r0 <= r2 && r0 >= r1)
        return 1.0;
    else
        return 0.0;
}

double fr1(double * pos, double  *ctx)
{
    /// if pos in [ctx0, ctx1] return |pos|
    double r0 = VecLen(pos,3);
    double r1 = ctx[0]/ctx[2];
    double r2 = ctx[1]/ctx[2];

    if(r0 <= r2 && r0 >= r1)
        return r0;
    else
        return 0.0;
}

int spa_filter(citcoms_sphere * _cs, double * ctx)
{
    /// calculate some value after integrate
    const double r1 = ctx[0]/ctx[2];
    const double r2 = ctx[1]/ctx[2];
    const int _noc = _cs->cap[0].noc;

    /*
     * the effective crust depth : data[2]
     * the average of crust fraction : data[0]
     * volume of integrate area : data[1]
     * the average depth of crust materials : data[3]
     */
    for(int k=0; k<_cs->nproc_surf;++k)
    {
        citcoms_sphere_dump * _csd = _cs->cap + k;
        for(int j=0; j<_csd->nel; ++j)
        {
            double avg_frac = _csd->data[_noc*j + 0]/_csd->data[_noc*j + 1];
            double avg_depth = _csd->data[_noc*j + 3]/_csd->data[_noc*j + 0];
            _csd->data[_noc*j + 0] = _csd->data[_noc*j + 0]/_csd->data[_noc*j + 1];
            _csd->data[_noc*j + 3] = _csd->data[_noc*j + 3]/_csd->data[_noc*j + 0];
            double eff_thick = pow(r2,3) - avg_frac * (pow(r2,3) - pow(r1,3));

            eff_thick = 1.0 - pow(eff_thick,1.0/3.0);
            avg_depth = 1.0 - avg_depth;

            _csd->data[_noc*j + 0] = avg_frac;
            _csd->data[_noc*j + 2] = ctx[2] * eff_thick;
            _csd->data[_noc*j + 3] = _csd->data[_noc*j + 2] > 200.0 ? ctx[2] * avg_depth : 0.0f;
        }
    }

    /*
     * calculate some integrate on surface
     *
     */
    double sum_area = 0.0;
    double spa_area = 0.0;
    double spa_center[3] = {0., 0., 0.};
    for(int k=0; k<_cs->nproc_surf;++k)
    {
        citcoms_sphere_dump * _csd = _cs->cap + k;
        const int nox = _csd->nox;
        const int noy = _csd->noy;
        const int noc = _csd->noc;

        for(int ix=0;ix<nox-1;++ix)
        {
            for(int jy=0;jy<noy-1;++jy)
            {
                int eid[4] = {0, ix+1, jy+1, 1};
                int n2ien[4] = {citcoms_offset(eid[1]+1,eid[2]+1,eid[3]+0,nox, noy, 1) - 1,
                               citcoms_offset(eid[1]+0,eid[2]+1,eid[3]+0,nox, noy, 1) - 1,
                               citcoms_offset(eid[1]+0,eid[2]+0,eid[3]+0,nox, noy, 1) - 1,
                               citcoms_offset(eid[1]+1,eid[2]+0,eid[3]+0,nox, noy, 1) - 1};

                double Xg[4][3];
                for(int j=0; j<4; ++j)
                {
                    for(int i=0; i<3; ++i)
                    {
                        Xg[j][i] = _csd->pos[3*n2ien[j] + i];
                    }
                }

                double res[4];
                DeriveArea(Xg, res);
                double element_area = VecLen(res,3);
                sum_area += element_area;

                const int n2 = citcoms_offset(ix+1,jy+1,1,nox-1,noy-1,1)-1;

                if(_csd->data[n2*noc + 2] <= 30.0e3)
                {
                    spa_area += element_area;
                    for(int j=0; j<4; ++j)
                    {
                        for(int i=0; i<3;++i)
                        {
                            // spa_center[i] += _csd->pos[3*n2ien[j] + i] * element_area * 0.25;
                            spa_center[i] += Xg[j][i] * 0.25* element_area;
                        }
                    }
                }
            }
        }
    }

    fprintf(stdout, "%s: sum area is %f*PI\n",__func__, sum_area/M_PI);
    fprintf(stdout,"%s: spa area is %f*PI approx %f m\n",__func__,spa_area/M_PI, sqrt(spa_area/M_PI)*ctx[2]);

    for(int i=0;i<3;++i)
        spa_center[i] /= spa_area;
    const double r_spa = sqrt(spa_area/M_PI);

    float x0 = spa_center[0], y0=spa_center[1], z0=spa_center[2];
    float r0 = sqrtf(x0*x0 + y0*y0 + z0*z0);
    x0 /= r0;
    y0 /= r0;
    z0 /= r0;

    for(int k=0; k<_cs->nproc_surf;++k)
    {
        citcoms_sphere_dump *_csd = _cs->cap + k;
        const int nox = _csd->nox;
        const int noy = _csd->noy;
        const int noc = _csd->noc;
        const int nno = _csd->nno;

        for(int j=0; j< nno; ++j)
        {
            float x = _csd->pos[j*3 + 0], y=_csd->pos[j*3 + 1], z=_csd->pos[j*3 + 2];
            float rx = sqrtf(x*x + y*y + z*z);
            /// cartesian distance
            double dist = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0));
            _csd->pdata[j*noc + 1] =  dist/r_spa;

            /// sphere distance
            double angle = acosf((x*x0 + y*y0 + z*z0)/rx);
            _csd->pdata[j*noc + 0] =  angle/r_spa;

            _csd->pdata[j*noc + 2] = angle > 2.5*r_spa ? 0.0 : 1.0; // mask used in temperature
        }
    }

    fprintf(stdout,"%s: center of spa (%e,%e,%e)\n",__func__,x0,y0,z0);
    return 0;
}

int apply_spa_filter(citcoms_dump * _cd, citcoms_sphere * _cs)
{
    for(int k=0; k<_cd->nproc;++k)
    {
        citcoms_temp_dump  * _ctd = _cd->temp + k;
        citcoms_sphere_dump *_csd = _cs->cap + k/_cd->nprocz;
        const int nox = _ctd->nox;
        const int noy = _ctd->noy;
        const int noz = _ctd->noz;

        assert(nox == _csd->nox);
        assert(noy == _csd->noy);
        assert(noz == _csd->noz);

        const int noc = _csd->noc;

        for(int ix=0;ix<nox;++ix)
        {
            for(int jy=0;jy<noy;++jy)
            {
                for(int kz=0;kz<noz;++kz)
                {
                    //const int n2 = jy + noy*ix;
                    const int n2 = citcoms_offset(ix+1,jy+1,1,nox,noy,1)-1;
                    const int n3 = citcoms_offset(ix+1,jy+1,kz+1,nox,noy,noz);
                    // _ctd->data[n3] = 100.0*_csd->pdata[n2*noc + 2];
                    // _ctd->data[n3] *= _csd->pdata[n2*noc + 2];
                }
            }
        }
    }
    return  0;
}

int VIntCitcomsTempDump(citcoms_temp_dump * _ctd, float * dump, int len_dump,double (*f[])(double *,double *),  double * ctx)
{
    const int nodes = _ctd->nox * _ctd->noy;
    const int elements = (_ctd->nox-1) * (_ctd->noy-1);
    /// dump should be allocated before called
    if(nodes == len_dump)
    {
        /// integrate along lines (z direction)
        return 1;
    }

    if(elements == len_dump)
    {
        /// calculate in elements
        for(int ix=0;ix<_ctd->nox-1;++ix)
            for(int jy=0;jy<_ctd->noy-1;++jy)
                for(int kz=0; kz<_ctd->noz-1; ++kz)
                {
                    int eid[4] = {0, ix+1, jy+1, kz+1};
                    int shape[3] = {_ctd->nox, _ctd->noy, _ctd->noz};
                    double res[3] = {0};
                    StructedGridIntf(_ctd->X, _ctd->crust, shape, eid, res, f[0], ctx);
                    const int noy = _ctd->noy, nox = _ctd->nox;
                    const int n2 = citcoms_offset(ix+1,jy+1,1,nox-1,noy-1,1)-1;
                    dump[n2*4 + 0] += (float) res[0];
                    dump[n2*4 + 1] += (float) res[1];
                    dump[n2*4 + 2] += (float) res[2];
                    StructedGridIntf(_ctd->X, _ctd->crust, shape, eid, res, f[1], ctx);
                    dump[n2*4 + 3] += (float) res[0];
                }
        return 2;
    }
    return 0;
}


int VIntCitcomsTempDump2(VtsInfo * _vsf, float * dump, int len_dump,double (*f[])(double *,double *),  double * ctx)
{
    const unsigned int nox = _vsf->Nxp[0];
    const unsigned int noy = _vsf->Nxp[1];
    const unsigned int noz = _vsf->Nxp[2];
    const unsigned int nodes =  nox*noy;
    const unsigned int elements = (nox-1)*(noy-1);
    /// dump should be allocated before called
    if(nodes == len_dump)
    {
        /// integrate along lines (z direction)
        return 1;
    }

    if(elements == len_dump)
    {
        int coord_fId = find_pointfield("coordinate", _vsf);
        int comp1_fId = find_pointfield("composition1",_vsf);
        /// calculate in elements
        for(int ix=0;ix<nox-1;++ix)
            for(int jy=0;jy<noy-1;++jy)
                for(int kz=0; kz<noz-1; ++kz)
                {
                    int eid[4] = {0, ix+1, jy+1, kz+1};
                    int shape[3] = {nox, noy, noz};
                    double res[3] = {0};
                    StructedGridIntf2(_vsf->PointField[coord_fId].Data, _vsf->PointField[comp1_fId].Data, shape, eid, res, f[0], ctx);
                    //const int n2 = jy + (noy - 1)*ix;
                    const int n2 = citcoms_offset(ix+1,jy+1,1,nox-1,noy-1,1)-1;
                    dump[n2*4 + 0] += (float) res[0];
                    dump[n2*4 + 1] += (float) res[1];
                    dump[n2*4 + 2] += (float) res[2];
                    StructedGridIntf2(_vsf->PointField[coord_fId].Data, _vsf->PointField[comp1_fId].Data, shape, eid, res, f[1], ctx);
                    dump[n2*4 + 3] += (float) res[0];
                }
        return 2;
    }
    return 0;
}


int SphereIntegrateCitcomsDump(citcoms_dump * _cd, const char * _prefix)
{
    /// (1) calculate integrate in r-interval,
    /// (2) write results to vtm file

    /// (1) calculate
    const int _noc = 4;
    citcoms_sphere * _cs = init_citcoms_sphere(_cd, _noc);

    /// calculate volume of crust between ctx[0], ctx[1]
    const double depth = 400.0e3;
    const double Rm = _cd->TransformR;
    double ctx[] = {Rm - depth, Rm, Rm};

    /// functions used for Integrate
    double (*fctx[])(double*, double*) = {fr0, fr1};

    // for(int k=0; k<_cd->nproc; ++k)
    // {
    //     citcoms_temp_dump * _ctd = _cd->temp + k;
    //     const int cap_id = k/ _cd->nprocz;
    //     VIntCitcomsTempDump(_ctd, _cs->cap[cap_id].data, _cs->cap[cap_id].nel,fctx, ctx);
    // }

    #pragma omp parallel for num_threads(LOADTHREADS) default(shared)
    for(int cap_id = 0; cap_id < _cs->nproc_surf; ++cap_id)
    {
        for(int k=0; k<_cd->nprocz; ++k)
        {
            citcoms_temp_dump * _ctd = _cd->temp +  _cd->nprocz*cap_id + k;
            VIntCitcomsTempDump(_ctd, _cs->cap[cap_id].data, _cs->cap[cap_id].nel,fctx, ctx);
        }
    }


    spa_filter(_cs, ctx);
    apply_spa_filter(_cd, _cs);

    /// (2) write to vtm & clean
    write_citcoms_sphere(_cs, _prefix);
    clean_citcoms_sphere(_cs);
    return 0;
}

int SphereIntegrateCitcomsDump2(CitcomsData * _cd, const char * _prefix)
{
    /// (1) calculate integrate in r-interval,
    /// (2) write results to vtm file

    /// (1) calculate
    const int _noc = 4;
    citcoms_sphere * _cs = init_citcoms_sphere2(_cd, _noc);

    /// calculate volume of crust between ctx[0], ctx[1]
    const double depth = 400.0e3;
    const double Rm = 1.74e6;
    double ctx[] = {Rm - depth, Rm, Rm};

    /// functions used for Integrate
    double (*fctx[])(double*, double*) = {fr0, fr1};

    // for(int k=0; k<_cd->nproc; ++k)
    // {
    //     citcoms_temp_dump * _ctd = _cd->temp + k;
    //     const int cap_id = k/ _cd->nprocz;
    //     VIntCitcomsTempDump(_ctd, _cs->cap[cap_id].data, _cs->cap[cap_id].nel,fctx, ctx);
    // }

    // #pragma omp parallel for num_threads(LOADTHREADS) default(shared)
    // for(int cap_id = 0; cap_id < _cs->nproc_surf; ++cap_id)
    // {
    //     for(int k=0; k<_cd->nprocz; ++k)
    //     {
    //         VtsInfo * _vsf = _cd->VSF +  _cd->nprocz*cap_id + k;
    //         VIntCitcomsTempDump2(_vsf, _cs->cap[cap_id].data, _cs->cap[cap_id].nel,fctx, ctx);
    //     }
    // }

    // spa_filter(_cs, ctx);

    //apply_spa_filter(_cd, _cs);
    set_ring_scope(_cs,"ring.txt");
    int calculate_effective_depth(citcoms_sphere * _cs, CitcomsData * _cd);
    calculate_effective_depth(_cs, _cd);
    
    /// (2) write to vtm & clean
    write_citcoms_sphere(_cs, _prefix);
    clean_citcoms_sphere(_cs);
    return 0;
}

int check_in_polygon(double * x, double * pl, int n)
{
    int in_polygon = 1;
    for(int k=0; k<n;++k)
    {
        double * A = pl + 3*k;
        double * B = pl + 3*((k+1)%n);

        double AxB[3];
        VecCross(AxB, A, B, 3);
        VecNormalize(AxB,3);

        double xdAxB = VecDot(x, AxB, 3);

        if(xdAxB < -1.0e-5)
        {
            in_polygon = 0;
            break;
        }
    }
    return in_polygon;
}

int set_polygon_marker(double * pl, int n, citcoms_sphere * _cs)
{
    for(int k=0;k<_cs->nproc_surf;++k)
    {
        citcoms_sphere_dump * _csd = _cs->cap + k;

        const int nno = _csd->nno;
        const int nox = _csd->nox;
        const int noy = _csd->noy;
        const int noc = _csd->noc;
        /// set marker on points
        for(int i=0;i<nno;++i)
        {
            double ipos[3] = {_csd->pos[3*i + 0], _csd->pos[3*i + 1], _csd->pos[3*i + 2]};
            _csd->pdata[i*noc + noc - 1] += check_in_polygon(ipos, pl, n);
        }

        /// set on integral points
        for(int ix=0;ix<nox-1;++ix)
        {
            for(int jy = 0; jy < noy - 1; ++jy)
            {
                int eid[4] = {0, ix + 1, jy + 1, 1};
                int n2ien[4] = {citcoms_offset(eid[1] + 1, eid[2] + 1, eid[3] + 0, nox, noy, 1) - 1,
                                citcoms_offset(eid[1] + 0, eid[2] + 1, eid[3] + 0, nox, noy, 1) - 1,
                                citcoms_offset(eid[1] + 0, eid[2] + 0, eid[3] + 0, nox, noy, 1) - 1,
                                citcoms_offset(eid[1] + 1, eid[2] + 0, eid[3] + 0, nox, noy, 1) - 1};

                double Xg[4][3];
                for(int j = 0; j < 4; ++j)
                {
                    for(int i = 0; i < 3; ++i)
                    {
                        Xg[j][i] = _csd->pos[3 * n2ien[j] + i];
                    }
                }

                double Xgi[4][3];
                for(int j=0;j<4;++j)
                {
                    XgIpB(Xgi[j],Xg,GIPS2d[j]);
                }

                const int n2 = citcoms_offset(ix+1,jy+1,1,nox-1,noy-1,1) -1;
                for(int j=0;j<4;++j)
                {
                    _csd->marker[n2*4 + j] += check_in_polygon(Xgi[j], pl, n);
                }
            }
        }
    }
    return 0;
}

int set_ring_scope(citcoms_sphere * _cs, const char * fname)
{
    FILE * fp = fopen(fname,"r");
    assert(fp!=NULL);
    int npos = 0;
    fscanf(fp,"%d", &npos);
    assert(npos>=1);

    /// points coordinates
    double * pts = malloc(sizeof(double)*npos*3);
    for(int k=0;k<npos;++k)
    {
        int ncoord = fscanf(fp,"%lf %lf %lf", pts+3*k+0, pts+3*k+1, pts+3*k+2);
        assert(ncoord==3);
        VecNormalize(pts+3*k, 3);
    }

    /// convex chord
    int nchords = 0;
    fscanf(fp,"%d",&nchords);
    assert(nchords >= 2);
    int * chords = malloc(sizeof(int)*nchords);
    if(nchords > 0)
    {
        for(int k=0;k<nchords;++k)
        {
            fscanf(fp,"%d",chords+k);
        }
    }

    /// build convex chords and check
    for(int k=0;k<nchords;++k)
    {
        int spt = chords[k];
        int ept = chords[(k + 1) % nchords];
        assert(spt != ept);
        int npt = 0;
        double *pt = NULL;
        if(spt < ept)
        {
            npt = ept - spt + 1;
        }
        else
        {
            npt = npos - (spt - ept - 1);
        }

        pt = malloc(sizeof(double) * npt * 3);
        for(int j = 0; j < npt; ++j)
        {
            int pindex = (spt + j)%npos;
            pt[j * 3 + 0] = pts[pindex * 3 + 0];
            pt[j * 3 + 1] = pts[pindex * 3 + 1];
            pt[j * 3 + 2] = pts[pindex * 3 + 2];
        }

        /// check in polygon
        set_polygon_marker(pt, npt, _cs);

        if(NULL!=pt)
            free(pt);
    }

    if(nchords >= 3)
    {
        int spt = chords[0];
        int ept = chords[nchords - 1];
        assert(spt != ept);
        int npt = nchords;
        double *pt = malloc(sizeof(double) * npt * 3);

        for(int j=0;j<nchords;++j)
        {
            pt[j*3 + 0] = pts[chords[j]*3 + 0];
            pt[j*3 + 1] = pts[chords[j]*3 + 1];
            pt[j*3 + 2] = pts[chords[j]*3 + 2];
        }

        /// check in polygon
        set_polygon_marker(pt, npt, _cs);

        if(NULL!=pt)
            free(pt);
    }

    free(pts);
    free(chords);
    fclose(fp);
    return 0;
}


int vts_calculate_effective_depth(citcoms_sphere_dump * _csd, VtsInfo * _vsf)
{
    const int nox = _vsf->Nxp[0];
    const int noy = _vsf->Nxp[1];
    const int noz = _vsf->Nxp[2];
    const int nodes =  nox*noy;
    const int elements = (nox-1)*(noy-1);

    int coord_fId = find_pointfield("coordinate", _vsf);
    int comp1_fId = find_pointfield("composition1",_vsf);
    int melting_fId = find_pointfield("melting",_vsf);

    float * points = _vsf->PointField[coord_fId].Data;
    float * crust  = _vsf->PointField[comp1_fId].Data;
    float * melting = _vsf->PointField[melting_fId].Data;

    /// calculate in elements
    for(int ix=0;ix<nox-1;++ix)
    {
        for(int jy=0;jy<noy-1;++jy)
            for(int kz=0; kz<noz-1; ++kz)
            {
                int eid[4] = {0, ix+1, jy+1, kz+1};
                int shape[3] = {nox, noy, noz};

                int ien[8] = {
                        citcoms_offset(eid[1]+1,eid[2]+1,eid[3]+0,shape[0],shape[1],shape[2]) - 1,
                        citcoms_offset(eid[1]+0,eid[2]+1,eid[3]+0,shape[0],shape[1],shape[2]) - 1,
                        citcoms_offset(eid[1]+0,eid[2]+0,eid[3]+0,shape[0],shape[1],shape[2]) - 1,
                        citcoms_offset(eid[1]+1,eid[2]+0,eid[3]+0,shape[0],shape[1],shape[2]) - 1,
                        citcoms_offset(eid[1]+1,eid[2]+1,eid[3]+1,shape[0],shape[1],shape[2]) - 1,
                        citcoms_offset(eid[1]+0,eid[2]+1,eid[3]+1,shape[0],shape[1],shape[2]) - 1,
                        citcoms_offset(eid[1]+0,eid[2]+0,eid[3]+1,shape[0],shape[1],shape[2]) - 1,
                        citcoms_offset(eid[1]+1,eid[2]+0,eid[3]+1,shape[0],shape[1],shape[2]) - 1,
                };

                int n2ien[4] = {citcoms_offset(eid[1]+1,eid[2]+1,1,nox, noy, 1) - 1,
                                citcoms_offset(eid[1]+0,eid[2]+1,1,nox, noy, 1) - 1,
                                citcoms_offset(eid[1]+0,eid[2]+0,1,nox, noy, 1) - 1,
                                citcoms_offset(eid[1]+1,eid[2]+0,1,nox, noy, 1) - 1};

                const int n2 = citcoms_offset(ix+1,jy+1,1,nox-1,noy-1,1)-1;
                double Xg[8][3];

                double Vg[8], Yg[8];
                double Hg[8], Mg[8], Eg[8], C1g[8];

                // Hg : the depth of points
                // Mg : marker of internal of spa
                // C1g : composition 1
                // Eg : melting frac

                double HC1g[8];

                for(int k=0;k<8;++k)
                {
                    Xg[k][0] = points[ien[k]*3 + 0];
                    Xg[k][1] = points[ien[k]*3 + 1];
                    Xg[k][2] = points[ien[k]*3 + 2];

                    Yg[k] = 1.0;
                    Vg[k] = 1.0;
                    Mg[k] = _csd->pdata[_csd->noc*n2ien[k%4] + _csd->noc - 1];
                    Hg[k] = 1.0 - VecLen(Xg[k], 3);

                    C1g[k] = crust[ien[k]];
                    Eg[k] = melting[ien[k]];

                    HC1g[k] = Hg[k] * C1g[k];
                }

                double Ires[4] = {0., 0., 0., 0.};

                for(int k=0;k<NIpV;k++)
                {
                    // if(_csd->marker[n2*4 + k%4] < 1)
                    //     continue;

                    // double tmp[3];
                    // XgIpB(tmp, Xg, GIPS3d[k]);
                    // if(VecLen(tmp,3) < 0.3)
                    //     continue;

                    Ires[1] += GIWS3d[k] * detJV(Xg,GIPS3d[k]); // integral volume
                    double Egk = DataIpV(Eg,GIPS3d[k]);
                    if(Egk > 0.01)
                    {
                        Ires[0] += GIWS3d[k] * detJV(Xg,GIPS3d[k]) * DataIpV(C1g,GIPS3d[k]); // fraction sum
                        Ires[2] += GIWS3d[k] * detJV(Xg,GIPS3d[k]) * DataIpV(HC1g,GIPS3d[k]); // effective depth
                    }
                }
                _csd->data[n2*_csd->noc + 0] += Ires[0];
                _csd->data[n2*_csd->noc + 1] += Ires[1];
                _csd->data[n2*_csd->noc + 2] += Ires[2];
            }
    }
    return 0;
}

int calculate_effective_depth(citcoms_sphere * _cs, CitcomsData * _cd)
{
    #pragma omp parallel for num_threads(LOADTHREADS) default(shared)
    for(int cap_id = 0; cap_id < _cs->nproc_surf; ++cap_id)
    {
        for(int k=0; k<_cd->nprocz; ++k)
        {
            vts_calculate_effective_depth(_cs->cap + cap_id, _cd->VSF +  _cd->nprocz*cap_id + k);
        }
    }

    double sum_area = 0.0;
    for(int k=0; k<_cs->nproc_surf;++k)
    {
        citcoms_sphere_dump *_csd = _cs->cap + k;
        for(int j=0;j<_csd->nel;++j) sum_area += _csd->area[j];
    }
    fprintf(stdout,"sum area:%f*PI\n",sum_area/M_PI);

    /// divide volume
    for(int k=0; k<_cs->nproc_surf;++k)
    {
        citcoms_sphere_dump * _csd = _cs->cap + k;
        const int _noc = _csd->noc;
        for(int j=0; j<_csd->nel; ++j)
        {
            if(_csd->data[_noc*j + 1] <= 0.)
                continue;

            double r2 = 1.0;
            double r1 = 0.1954;
            double avg_frac = _csd->data[_noc*j + 0]/_csd->data[_noc*j + 1];
            double avg_depth = _csd->data[_noc*j + 2]/_csd->data[_noc*j + 0];
            double eff_thick = pow(r2,3) - avg_frac * (pow(r2,3) - pow(r1,3));

            eff_thick = 1.0 - pow(eff_thick,1.0/3.0);

            eff_thick *= 1.74e6;
            avg_depth *= 1.74e6;

            double eff_thick2 = _csd->data[_noc*j + 0]/_csd->area[j]*1.74e6;

            _csd->data[_noc*j + 0] = pow(1.0-3.0*_csd->data[_noc*j + 1]/_csd->area[j],1.0/3.0);
            _csd->data[_noc*j + 1] = eff_thick2;
            _csd->data[_noc*j + 2] = eff_thick;
            _csd->data[_noc*j + 3] = eff_thick > 500.0 ? avg_depth : 0.0f;
        }
    }

    return 0;
}