//
// Created by huachengli on 8/26/24.
//

#include "citcoms_related.h"
#include <omp.h>
#include "lmath.h"
#include "lmath3d8.h"
#include <sys/time.h>

int load_citcoms_temp_dump(citcoms_temp_dump * _ctd, const char * fname)
{
    FILE * fp = fopen(fname,"rb");
    int header[6];
    fread(header, sizeof(int), 6, fp);
    if(header[0] != 0x789)
    {
        fprintf(stdout,"%s:%s unmatched header %d/%d\n",__func__,fname,header[0],0x789);
        return 0;
    }

    _ctd->nox = header[1];
    _ctd->noy = header[2];
    _ctd->noz = header[3];
    _ctd->nno = header[4];
    _ctd->ncaps = header[5];

    assert(_ctd->ncaps > 0);

    _ctd->x = malloc(sizeof(double) * _ctd->ncaps * (_ctd->nno + 1));
    _ctd->y = malloc(sizeof(double) * _ctd->ncaps * (_ctd->nno + 1));
    _ctd->z = malloc(sizeof(double) * _ctd->ncaps * (_ctd->nno + 1));
    _ctd->data = malloc(sizeof(double) * _ctd->ncaps * (_ctd->nno + 1));
    _ctd->crust  = malloc(sizeof(double) * _ctd->ncaps * (_ctd->nno + 1));
    _ctd->dump  = malloc(sizeof(double) * _ctd->ncaps * (_ctd->nno + 1));

    for(int j=0;j<_ctd->ncaps;++j)
    {
        int dump_offset = j*_ctd->nno;
        fread(_ctd->x + dump_offset + 1, sizeof(double), _ctd->nno, fp);
        fread(_ctd->y + dump_offset + 1, sizeof(double), _ctd->nno, fp);
        fread(_ctd->z + dump_offset + 1, sizeof(double), _ctd->nno, fp);
        fread(_ctd->data + dump_offset, sizeof(double), _ctd->nno + 1, fp);
    }

    _ctd->X[0] = _ctd->x;
    _ctd->X[1] = _ctd->y;
    _ctd->X[2] = _ctd->z;

    fclose(fp);
    return _ctd->nno;
}

int clean_citcoms_temp_dump(citcoms_temp_dump * _ctd)
{
    free(_ctd->x);
    free(_ctd->y);
    free(_ctd->z);
    free(_ctd->data);
    free(_ctd->crust);
    free(_ctd->dump);
    return 0;
}

int load_citcoms_tracer_dump(citcoms_tracer_dump * _ctd, const char * fname)
{
    FILE * fp = fopen(fname,"rb");
    assert(NULL != fp);
    int header[7];
    fread(header, sizeof(int), 7, fp);
    if(header[0] !=  101 || header[1] != 120)
    {
        fprintf(stdout,"%s:%s unmatched header %d,%d/%d,%d\n",__func__,fname,header[0],header[1],101,120);
        return 0;
    }

    _ctd->ncaps = header[2];
    _ctd->num_basic_q = header[3];
    _ctd->num_extra_q = header[4];
    _ctd->nflavors = header[5];
    _ctd->itc = header[6];

    _ctd->ntracers = malloc(sizeof(int)*_ctd->ncaps);
    _ctd->basicq = malloc(sizeof(double *)*_ctd->ncaps);
    _ctd->extraq = malloc(sizeof(double *)*_ctd->ncaps);
    _ctd->ielement = malloc(sizeof(int *)*_ctd->ncaps);

    const size_t fix_num_basicq = 6;
    for(int j=0;j<_ctd->ncaps;++j)
    {
        fread(_ctd->ntracers + j, sizeof(int), 1, fp);
        _ctd->basicq[j] = malloc(sizeof(double)*fix_num_basicq*(_ctd->ntracers[j]+1));
        _ctd->extraq[j] = malloc(sizeof(double)*_ctd->num_extra_q*(_ctd->ntracers[j]+1));
        _ctd->ielement[j] = malloc(sizeof(int)*(_ctd->ntracers[j]+1));
    }

    size_t t1=0,t2=0,t3=0;
    for(int j=0;j<_ctd->ncaps;++j)
    {
        t1 = fread(_ctd->basicq[j], sizeof(double), fix_num_basicq*(_ctd->ntracers[j]+1), fp);
        t2 = fread(_ctd->extraq[j], sizeof(double), _ctd->num_extra_q*(_ctd->ntracers[j]+1), fp);
        t3 = fread(_ctd->ielement[j], sizeof(int), _ctd->ntracers[j]+1, fp);
        assert(t1 + t2 + t3 == (_ctd->num_extra_q + fix_num_basicq + 1)*(_ctd->ntracers[j]+1));
    }

    fclose(fp);
    return _ctd->nflavors;
}

int clean_citcoms_tracer_dump(citcoms_tracer_dump * _ctd)
{
    for(int j=0;j<_ctd->ncaps;++j)
    {
        free(_ctd->basicq[j]);
        free(_ctd->extraq[j]);
        free(_ctd->ielement[j]);
    }

    free(_ctd->ntracers);
    free(_ctd->basicq);
    free(_ctd->extraq);
    free(_ctd->ielement);

    return 0;
}

int load_citcoms_dump(citcoms_dump * _cd, InputFile * ifp)
{
    _cd->nproc = GetValueI(ifp,"citcoms.nproc","12");
    _cd->nproc_surf = GetValueI(ifp,"citcoms.nproc_surf","12");
    _cd->nprocx = GetValueI(ifp,"citcoms.nprocx","2");
    _cd->nprocy = GetValueI(ifp,"citcoms.nprocy","2");
    _cd->nprocz = GetValueI(ifp,"citcoms.nprocz","2");

    assert(_cd->nproc == _cd->nproc_surf*_cd->nprocx*_cd->nprocy*_cd->nprocz);
    _cd->nproc_surf *= _cd->nprocx*_cd->nprocy;

    _cd->TransformR = GetValueD(ifp,"citcoms.TransformR","1.0");
    _cd->pad_step = _cd->TransformR/2500.0;
    GetValueS(ifp,"citcoms.temperature_dump",_cd->temp_prefix,"NONE");
    GetValueS(ifp,"citcoms.tracer_dump",_cd->tracer_prefix,"NONE");
    assert(_cd->nproc > 1);
    _cd->temp = malloc(sizeof(citcoms_temp_dump)* _cd->nproc);
    _cd->tracer = malloc(sizeof(citcoms_tracer_dump)* _cd->nproc);
    assert(_cd->temp != NULL && _cd->tracer != NULL);
    fprintf(stdout,"%s:load citcoms temperature dump from %s.\n",__func__,_cd->temp_prefix);
    fprintf(stdout,"%s:load citcoms tracer dump from %s.\n",__func__,_cd->tracer_prefix);

    for(int j=0;j<_cd->nproc;++j)
    {
        char temp_dump_name[4097];
        char tracer_dump_name[4097];
        snprintf(temp_dump_name,4096,"%s.%d.dump",_cd->temp_prefix,j);
        snprintf(tracer_dump_name,4096,"%s.%d.dump",_cd->tracer_prefix,j);
        load_citcoms_temp_dump(_cd->temp+j,temp_dump_name);
        load_citcoms_tracer_dump(_cd->tracer+j,tracer_dump_name);
    }
    return _cd->nproc;
}

int clean_citcoms_dump(citcoms_dump * x)
{
    assert(x!=NULL && x->nproc > 0);
    for(int j=0;j<x->nproc;++j)
    {
        clean_citcoms_temp_dump(x->temp + j);
        clean_citcoms_tracer_dump(x->tracer + j);
    }
    free(x->temp);
    free(x->tracer);
    return 0;
}

citcoms_sphere * init_citcoms_sphere(citcoms_dump * _cd, int noc)
{
    citcoms_sphere * _cs = (citcoms_sphere *) malloc(sizeof(citcoms_sphere));
    assert(NULL != _cs);
    /// copy metadata from _cd
    _cs->nproc_surf = _cd->nproc_surf;
    _cs->nproc = _cd->nproc;
    _cs->nprocx = _cd->nprocx;
    _cs->nprocy = _cd->nprocy;
    _cs->nprocz = _cd->nprocz;
    /// allocate memory for sphere
    assert(noc >= 1);
    assert(_cs->nproc_surf >= 1);
    _cs->cap = (citcoms_sphere_dump *)malloc(sizeof(citcoms_sphere_dump)*_cs->nproc_surf);

    /// allocate sphere
    for(int k=0; k< _cs->nproc_surf; ++k)
    {
        const int iproc = k*_cd->nprocz + _cd->nprocz - 1;
        const int nox = _cd->temp[iproc].nox;
        const int noy = _cd->temp[iproc].noy;
        const int noz = _cd->temp[iproc].noz;
        const int nel = (nox - 1) * (noy - 1);
        const int nno = nox * noy;

        citcoms_temp_dump * _ctd = _cd->temp + iproc;

        _cs->cap[k].pos = (float *) malloc(sizeof(float) * nno * 3);
        _cs->cap[k].data = (float *) malloc(sizeof(float) * nel * noc);
        _cs->cap[k].pdata = (float *) malloc(sizeof(float) * nno * noc);

        _cs->cap[k].nno = nno;
        _cs->cap[k].nel = nel;
        _cs->cap[k].noc = noc;
        _cs->cap[k].nox = nox;
        _cs->cap[k].noy = noy;
        _cs->cap[k].noz = noz;

        /// copy coordinate to pos
        const int kz = noz - 1;
        for(int ix=0;ix<nox;++ix)
        {
            for(int jy=0;jy<noy;++jy)
            {
                const int n3 =citcoms_offset(ix+1,jy+1,kz+1,nox,noy,noz);
                const int n2 = citcoms_offset(ix+1,jy+1,1,nox,noy,1)-1;
                _cs->cap[k].pos[3*n2 + 0] = (float) _ctd->x[n3]*1.0001f;
                _cs->cap[k].pos[3*n2 + 1] = (float) _ctd->y[n3]*1.0001f;
                _cs->cap[k].pos[3*n2 + 2] = (float) _ctd->z[n3]*1.0001f;
            }
        }
        /// set initial value of data to 0.
        for(int j=0;j<noc*nel;++j)
            _cs->cap[k].data[j] = 0.0f;
        for(int j=0;j<noc*nno;++j)
            _cs->cap[k].pdata[j] = 8.0f;
        citcoms_sphere_dump * _csd = _cs->cap + k;
        _csd->area = _csd->marker = NULL;
    }
    return _cs;
}

citcoms_sphere * init_citcoms_sphere2(CitcomsData * _cd, int noc)
{
    citcoms_sphere * _cs = (citcoms_sphere *) malloc(sizeof(citcoms_sphere));
    assert(NULL != _cs);
    /// copy metadata from _cd
    _cs->nproc_surf = _cd->nproc_surf;
    _cs->nproc = _cd->nproc;
    _cs->nprocx = _cd->nprocx;
    _cs->nprocy = _cd->nprocy;
    _cs->nprocz = _cd->nprocz;
    /// allocate memory for sphere
    assert(noc >= 1);
    assert(_cs->nproc_surf >= 1);
    _cs->cap = (citcoms_sphere_dump *)malloc(sizeof(citcoms_sphere_dump)*_cs->nproc_surf);

    /// allocate sphere
    for(int k=0; k< _cs->nproc_surf; ++k)
    {
        const int iproc = k*_cd->nprocz + _cd->nprocz - 1;
        VtsInfo * _vsf = _cd->VSF + iproc;
        const int nox = _cd->nox;
        const int noy = _cd->noy;
        const int noz = _cd->noz;
        const int nel = (nox - 1) * (noy - 1);
        const int nno = nox * noy;

        _cs->cap[k].pos = (float *) malloc(sizeof(float) * nno * 3);
        _cs->cap[k].data = (float *) malloc(sizeof(float) * nel * noc);
        _cs->cap[k].pdata = (float *) malloc(sizeof(float) * nno * noc);

        _cs->cap[k].nno = nno;
        _cs->cap[k].nel = nel;
        _cs->cap[k].noc = noc;
        _cs->cap[k].nox = nox;
        _cs->cap[k].noy = noy;
        _cs->cap[k].noz = noz;

        /// copy coordinate to pos
        const unsigned int coord_fId = find_pointfield("coordinate",_vsf);
        assert(coord_fId < 100);
        /*for(int j=0;j<_vsf->PointNoF;++j)
        {
            fprintf(stdout,"%s(%s):%d\n",_vsf->PointField[j].Name,_vsf->PointField[j].Format,_vsf->PointField[j].NoC);
        }*/
        float * pts_coord = _vsf->PointField[coord_fId].Data;

        const int kz = noz - 1;
        for(int ix=0;ix<nox;++ix)
        {
            for(int jy=0;jy<noy;++jy)
            {
                const int n3 =citcoms_offset(ix+1,jy+1,kz+1,nox,noy,noz)-1;
                //const int n2 = jy + noy*ix;
                const int n2 = citcoms_offset(ix+1,jy+1,1,nox,noy,1)-1;
                _cs->cap[k].pos[3*n2 + 0] = (float) pts_coord[3*n3 + 0] * 1.0001f;
                _cs->cap[k].pos[3*n2 + 1] = (float) pts_coord[3*n3 + 1] * 1.0001f;
                _cs->cap[k].pos[3*n2 + 2] = (float) pts_coord[3*n3 + 2] * 1.0001f;
            }
        }

        /// set initial value of data to 0.
        for(int j=0;j<noc*nel;++j)
            _cs->cap[k].data[j] = 0.0f;
        for(int j=0;j<noc*nno;++j)
            _cs->cap[k].pdata[j] = 0.0f;

        _cs->cap[k].marker = (float *) malloc(sizeof(float) * nel * 4);
        for(int j=0;j<4*nel;++j)
            _cs->cap[k].marker[j] = 0.0f;
        /// record surface area
        citcoms_sphere_dump * _csd = _cs->cap + k;
        _cs->cap[k].area = (float *) malloc(sizeof(float) * nel);
        for(int ix=0;ix<nox-1;++ix)
        {
            for(int jy=0;jy<noy-1;++jy)
            {
                int eid[4] = {0, ix+1, jy+1, 1};
                int n2ien[4] = {citcoms_offset(eid[1]+1,eid[2]+1,eid[3]+0,nox, noy, 1) - 1,
                                citcoms_offset(eid[1]+0,eid[2]+1,eid[3]+0,nox, noy, 1) - 1,
                                citcoms_offset(eid[1]+0,eid[2]+0,eid[3]+0,nox, noy, 1) - 1,
                                citcoms_offset(eid[1]+1,eid[2]+0,eid[3]+0,nox, noy, 1) - 1};
                const int n2 = citcoms_offset(ix+1,jy+1,1,nox-1,noy-1,1)-1;

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
                _csd->area[n2] = VecLen(res,3);
            }
        }
        _csd->vstat = (float *) malloc(sizeof(float) * noz * noc * _cs->nprocz);
        for(int j=0;j<noc*noz*_cs->nprocz;++j)
        {
            _csd->vstat[j] = 0.;
        }
    }
    return _cs;
}

int clean_citcoms_sphere(citcoms_sphere * _cs)
{
    for(int k=0;k<_cs->nproc_surf;++k)
    {
        free(_cs->cap[k].pos);
        free(_cs->cap[k].data);
        free(_cs->cap[k].pdata);

        if(NULL!=_cs->cap[k].area)
            free(_cs->cap[k].area);
        if(NULL!=_cs->cap[k].marker)
            free(_cs->cap[k].marker);
        if(NULL!=_cs->cap[k].vstat)
            free(_cs->cap[k].vstat);
    }
    free(_cs->cap);
    free(_cs);
    return 0;
}

int write_citcoms_sphere(citcoms_sphere * _cs, CitcomsData * _cd,const char * _name)
{
    try_make_dir("vts");
    char vtm_name[4097];
    snprintf(vtm_name, 4096, "%s.vtm", _name);
    FILE * vtm_fp = fopen(vtm_name,"w");
    assert(vtm_fp != NULL);

    if(_cd!=NULL && _cd->len_attach>0)
    {
        const char header[] =
                "<?xml version=\"1.0\"?>\n"
                "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" compressor=\"vtkZLibDataCompressor\" byte_order=\"LittleEndian\">\n";
        fputs(header, vtm_fp);
        for(int k=0;k<_cd->len_attach;++k)
        {
            fputs(_cd->attach[k], vtm_fp);
        }
        fputs("  <vtkMultiBlockDataSet>\n", vtm_fp);
    }
    else
    {
        const char header[] =
                "<?xml version=\"1.0\"?>\n"
                "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" compressor=\"vtkZLibDataCompressor\" byte_order=\"LittleEndian\">\n"
                "  <vtkMultiBlockDataSet>\n";
        fputs(header, vtm_fp);
    }

    for(int k=0; k<_cs->nproc_surf; ++k)
    {
        citcoms_sphere_dump * _csd = _cs->cap + k;
        char vts_name[4097];
        snprintf(vts_name,4096,"vts/%s.%04d.vts",_name,k);
        fprintf(vtm_fp, "    <DataSet index=\"%d\" file=\"%s\"/>\n",k,vts_name);
        FILE * vts_fp = fopen(vts_name,"w");
        assert(vts_fp != NULL);
        char whole_extent[4096], piece_extent[4096];
        snprintf(whole_extent, 4096, "%d %d %d %d 0 0", 1, _csd->nox, 1, _csd->noy);
        snprintf(piece_extent, 4096, "%d %d %d %d 0 0", 1, _csd->nox, 1, _csd->noy);
        vts_file_header(vts_fp, piece_extent, whole_extent);
        vtk_point_data_header(vts_fp);
        vtk_dataarray_vec_f(vts_fp, "dist", "binary", _csd->pdata, _csd->nno, _csd->noc);
        vtk_point_data_trailer(vts_fp);
        vtk_cell_data_header(vts_fp);
        vtk_dataarray_vec_f(vts_fp, "dump", "binary", _csd->data, _csd->nel, _csd->noc);
        if(_csd->marker != NULL)
        {
            vtk_dataarray_vec_f(vts_fp, "marker", "binary", _csd->marker, _csd->nel, 4);
        }
        if(_csd->area != NULL)
        {
            vtk_dataarray_vec_f(vts_fp, "area", "binary", _csd->area, _csd->nel, 1);

        }
        vtk_cell_data_trailer(vts_fp);
        vtk_point_header(vts_fp);
        vtk_dataarray_vec_f(vts_fp, "coordinate", "binary", _csd->pos, _csd->nno, 3);
        vtk_point_trailer(vts_fp);
        vts_file_trailer(vts_fp);
        fclose(vts_fp);
    }
    fputs("  </vtkMultiBlockDataSet>\n", vtm_fp);
    fputs("</VTKFile>", vtm_fp);
    fclose(vtm_fp);
    return 0;
}

citcoms_dump * InitCitcomsDump(InputFile * ifp)
{
    citcoms_dump * x = malloc(sizeof(citcoms_dump));
    assert(x != NULL);
    load_citcoms_dump(x,ifp);
    return x;
}

int CloseCitcomsDump(citcoms_dump * x)
{
    clean_citcoms_dump(x);
    free(x);
    return 0;
}

SALEcData * CrInitSALEcData(InputFile * ifp)
 {
     char _salc_inp_name[4096];
     GetValueS(ifp,"SALEc.input",_salc_inp_name,"SALEc.inp");
     SALEcData * _sdata = malloc(sizeof(SALEcData));
     assert(_sdata != NULL);
     LoadInpInfo(_sdata,_salc_inp_name);
     char _salec_data_path[4097];
     char _salec_data_prefix[4096];
     GetValueS(ifp,"Citcoms.data",_salec_data_prefix,"ParaTest");
     int _salec_data_step = GetValueI( ifp,"Citcoms.step", "0");
     snprintf(_salec_data_path,4096,"%s.proc%%d.%d.vts",_salec_data_prefix,_salec_data_step);
     fprintf(stdout,"%s:load initial temperature from %s.\n", __func__, _salec_data_path);
     LoadVtsData(_sdata,_salec_data_path);
     fputs("\n",stdout);
     return _sdata;
 }

SALEcData * CrInitSALEcData_ref(InputFile * ifp)
{
    assert(NULL != ifp);
    int _salec_data_step = GetValueI( ifp,"Citcoms.ref", "-1");
    if(-1 == _salec_data_step) return NULL;

    char _salc_inp_name[4096];
    GetValueS(ifp,"SALEc.input",_salc_inp_name,"SALEc.inp");
    SALEcData * _sdata = malloc(sizeof(SALEcData));
    assert(_sdata != NULL);
    LoadInpInfo(_sdata,_salc_inp_name);
    char _salec_data_path[4097];
    char _salec_data_prefix[4096];
    GetValueS(ifp,"Citcoms.data",_salec_data_prefix,"ParaTest");
    snprintf(_salec_data_path,4096,"%s.proc%%d.%d.vts",_salec_data_prefix,_salec_data_step);
    fprintf(stdout,"%s:load reference temperature from %s.\n", __func__, _salec_data_path);
    LoadVtsData(_sdata,_salec_data_path);
    fputs("\n",stdout);
    return _sdata;
}

void CrCloseSALEcData(SALEcData * _sdata)
{
    if(NULL != _sdata)
    {
        CleanSALEcData(_sdata);
        free(_sdata);
    }
}

int UpdateCitcomsTempDump(citcoms_dump * _cd, SALEcData * _sdata, SALEcData * _rdata)
{
    // locate tem field id
    VtsInfo * _vsf = _sdata->VSF;
    unsigned long fId= find_cellfield("Temperature",_vsf);
    unsigned long fId_vaccum = find_cellfield("VOF-0",_vsf);
    unsigned long fId_mantle = find_cellfield("VOF-1",_vsf);
    const unsigned long NoC = _vsf->CellField[fId].NoC;
    assert(NoC == 1);

    int IdList[3] = {fId_vaccum, fId, fId_mantle};
    fprintf(stdout,"%s:update citcoms temperature dump.\n",__func__);

    int invalid_pts = 0;
    #pragma omp parallel for num_threads(32) reduction(+:invalid_pts) shared(_cd,_sdata,_rdata,IdList) default(none)
    for(int p=0;p<_cd->nproc;++p)
    {
        citcoms_temp_dump * _ctd = _cd->temp + p;
        const int IdListLen = sizeof(IdList)/ sizeof(IdList[0]);
        assert(_ctd != NULL);
        for(int j=0;j<_ctd->ncaps;++j)
        {
            int dump_offset = j*_ctd->nno;
            for(int i=1;i<=_ctd->nno;++i)
            {

                VTSDATAFLOAT _pos[3] = {
                        (float ) (_ctd->x[dump_offset + i]*_cd->TransformR),
                        (float ) (_ctd->y[dump_offset + i]*_cd->TransformR),
                        (float ) ((_ctd->z[dump_offset + i] - 1.0)*_cd->TransformR)
                };

                VTSDATAFLOAT * _data = malloc(sizeof(VTSDATAFLOAT)*(IdListLen + 1));
                VTSDATAFLOAT * _data_r = malloc(sizeof(VTSDATAFLOAT)*(IdListLen+1)); // value for reference data
                SALEcGetCDataN(_sdata,IdList,IdListLen,_pos,_data,NULL);

                if(_rdata != NULL)
                {
                    SALEcGetCDataN(_rdata,IdList,IdListLen,_pos,_data_r,NULL);
                }
                else
                {
                    _data_r[0] = 0.f;
                    _data_r[1] = 0.f;
                }

                if(_data[0] > 0.05 || _data_r[0] > 0.05)
                {
                    // padding data from internal
                    double d0 = sqrt(_pos[0]*_pos[0] + _pos[1]*_pos[1] + (_pos[2]+_cd->TransformR)*(_pos[2]+_cd->TransformR))/_cd->TransformR;
                    VTSDATAFLOAT pad_vec[3] = {_pos[0]/d0, _pos[1]/d0, (_pos[2]+_cd->TransformR)/d0};
                    VTSDATAFLOAT pad_step = _cd->TransformR / 2500.0;
                    for(int k=0;k<1500;++k)
                    {
                        for(int q=0;q<3;++q) _pos[q] = _pos[q] - pad_step*pad_vec[q];

                        if(_data[0] > 0.05)
                            SALEcGetCDataN(_sdata, IdList, IdListLen, _pos, _data,NULL);

                        if(_data_r[0] > 0.05)
                            SALEcGetCDataN(_rdata,IdList,IdListLen,_pos,_data_r,NULL);

                        if(_data[0] <= 0.05 && _data_r[0] <= 0.05){
                            break;
                        }
                    }
                    invalid_pts++;
                }

                double * _data_d = _ctd->data + dump_offset + i*NoC;
                for(int kc=0;kc<NoC;++kc){
                    _data_d[kc] = _data[kc + 1] - _data_r[kc + 1];
                }

                double * _crust = _ctd->crust + dump_offset + i; // assuming NoC of vof-2 is 0
                _crust[0] = _data[2];
                double * _crust_r = _ctd->dump + dump_offset + i;
                _crust_r[0] = _data_r[2];
                free(_data);
                free(_data_r);
            }
        }
    }
    return invalid_pts;
}


int StructedGridIntf(double ** X, double * _data, int *shape, int *eid, double * res, double (*f)(double *, double *), double * ctx)
{
    int ien[8] = {
            citcoms_offset(eid[1]+1,eid[2]+1,eid[3]+0,shape[0],shape[1],shape[2]),
            citcoms_offset(eid[1]+0,eid[2]+1,eid[3]+0,shape[0],shape[1],shape[2]),
            citcoms_offset(eid[1]+0,eid[2]+0,eid[3]+0,shape[0],shape[1],shape[2]),
            citcoms_offset(eid[1]+1,eid[2]+0,eid[3]+0,shape[0],shape[1],shape[2]),
            citcoms_offset(eid[1]+1,eid[2]+1,eid[3]+1,shape[0],shape[1],shape[2]),
            citcoms_offset(eid[1]+0,eid[2]+1,eid[3]+1,shape[0],shape[1],shape[2]),
            citcoms_offset(eid[1]+0,eid[2]+0,eid[3]+1,shape[0],shape[1],shape[2]),
            citcoms_offset(eid[1]+1,eid[2]+0,eid[3]+1,shape[0],shape[1],shape[2]),
    };

    double Xg[8][3];
    double Yg[8];
    double Vg[9];

    for(int k=0;k<8;++k)
    {
        Xg[k][0] = X[0][ien[k]];
        Xg[k][1] = X[1][ien[k]];
        Xg[k][2] = X[2][ien[k]];

        Yg[k] = 1.0;
        Vg[k] = 1.0;

        if(NULL != f)
        {
            Yg[k] *= f(Xg[k], ctx);
            Vg[k] *= f(Xg[k], ctx);
        }

        if(NULL != _data)
        {
            Yg[k] *= _data[ien[k]];
        }
    }

    // res[0] = fabs(DeriveVolumeAverage(Yg,Xg));
    // res[1] = fabs(DeriveVolume(Xg));
    res[0] = DeriveVolumeAverage(Yg,Xg);
    res[1] = DeriveVolumeAverage(Vg,Xg);
    // res[1] = DeriveVolume(Xg);
    res[2] = 0.5*(Yg[0] + Yg[6]);
    return 2;
}

int StructedGridIntf2(float * X, float * _data, int *shape, int *eid, double * res, double (*f)(double *, double *), double * ctx)
{
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

    double Xg[8][3];
    double Yg[8];
    double Vg[9];

    for(int k=0;k<8;++k)
    {
        Xg[k][0] = X[ien[k]*3 + 0];
        Xg[k][1] = X[ien[k]*3 + 1];
        Xg[k][2] = X[ien[k]*3 + 2];

        Yg[k] = 1.0;
        Vg[k] = 1.0;

        if(NULL != f)
        {
            Yg[k] *= f(Xg[k], ctx);
            Vg[k] *= f(Xg[k], ctx);
        }

        if(NULL != _data)
        {
            Yg[k] *= _data[ien[k]];
        }
    }

    res[0] = DeriveVolumeAverage(Yg,Xg);
    res[1] = DeriveVolumeAverage(Vg,Xg);
    res[2] = 0.5*(Yg[0] + Yg[6]);
    return 2;
}



int CheckCitcomsTracerDump(citcoms_dump * _cd)
{
    /*
     * check the id of tracer flavor
     * map the matid from SALEc to CITCOMS
     */
    int check_num = 0;
    for(int p=0;p<_cd->nproc;++p)
    {
        citcoms_tracer_dump * _ctd = _cd->tracer + p;
        for(int j=0;j<_ctd->ncaps;++j)
        {
            for(int i=1;i<=_ctd->ntracers[j];++i)
            {
                // only for the crust + mantle model.
                if(_ctd->extraq[j][i] > 2.5 || _ctd->extraq[j][i] < 0.5)
                {
                    check_num++;
                    fprintf(stdout,"%s: not fill flavor in (%d,%d,%d)\n",__func__,p,j,i);
                }
                // convert to citcoms tracer flavor id
                _ctd->extraq[j][i] = 2.0 - _ctd->extraq[j][i];
            }
        }
    }
    return check_num;
}

int UpdateCitcomsTracerDump(citcoms_dump * _cd, SALEcData * _sdata)
{
    /*
     * fill tracers flavors in tracer dump data
     */
    VtsInfo * _vsf = _sdata->VSF;
    unsigned long Id0 = find_cellfield("VOF-0",_vsf);
    unsigned long Id1 = find_cellfield("VOF-1",_vsf);
    unsigned long Id2 = find_cellfield("VOF-2",_vsf);
    unsigned long Id3 = find_cellfield("VOF-3",_vsf);
    int IdList[4] = {Id0,Id1,Id2,Id3};
    assert(Id0 < 100 && Id1 < 100 && Id2 < 100 && Id3 < 100);

    int undetermined_tracer = 0, all_tracers=0;
    citcoms_tracer_mixed CTM;
    citcoms_tracer_mixed_init(&CTM);

    fprintf(stdout,"%s:update citcoms tracer dump.\n",__func__);

    #pragma omp parallel for num_threads(32) reduction(+:undetermined_tracer,all_tracers) shared(_cd,CTM,IdList,_sdata,stdout) default(none)
    for(int p=0;p<_cd->nproc;++p)
    {
        citcoms_tracer_dump * _ctd = _cd->tracer + p;
        for(int j=0;j<_ctd->ncaps;++j)
        {
            for(int i=1;i<=_ctd->ntracers[j];++i)
            {
                VTSDATAFLOAT _pos[3] = {
                        (VTSDATAFLOAT)(_ctd->basicq[j][i + 3*(_ctd->ntracers[j]+1)]*_cd->TransformR),
                        (VTSDATAFLOAT)(_ctd->basicq[j][i + 4*(_ctd->ntracers[j]+1)]*_cd->TransformR),
                        (VTSDATAFLOAT)((_ctd->basicq[j][i + 5*(_ctd->ntracers[j]+1)]-1.0)*_cd->TransformR)
                };

                VTSDATAFLOAT vof2[4];
                int mask[4];
                SALEcGetCDataN(_sdata, IdList, 4, _pos, vof2, mask);
                VTSDATAFLOAT d0 = sqrt(_pos[0]*_pos[0] + _pos[1]*_pos[1] + (_pos[2]+_cd->TransformR)*(_pos[2]+_cd->TransformR));

                if(vof2[0] > 0.05)
                {
                    // padding data from internal
                    VTSDATAFLOAT pad_vec[3] = {_pos[0]/d0, _pos[1]/d0, (_pos[2]+_cd->TransformR)/d0};
                    VTSDATAFLOAT pad_step = _cd->TransformR / 1750.0;
                    for(int k=0;k<1500;++k)
                    {
                        for(int q=0;q<3;++q) _pos[q] = _pos[q] - pad_step*pad_vec[q];
                        SALEcGetCDataN(_sdata, IdList, 4, _pos, vof2,mask);
                        if(vof2[0] <= 0.05){
                            break;
                        }
                    }
                }

                int Mid = VecMaxArgF(vof2, 4);
                if( 0 == Mid)
                {
                    fprintf(stdout,"%s: can not set vof at tracer[%d,%d,%d]\n",__func__,p,j,i);
                    _ctd->extraq[j][i] = -1.0f;
                    continue;
                }

                if(vof2[Mid] >= 0.99)
                {
                    _ctd->extraq[j][i] = 1.0f*Mid;
                }
                else
                {
                    #pragma omp critical
                    {
                        tracer_mixed tmp;
                        tmp.id[0] = p;
                        tmp.id[1] = j;
                        tmp.id[2] = i;
                        memcpy(tmp.mask,mask, sizeof(int)*4);
                        memcpy(tmp.vof,vof2, sizeof(VTSDATAFLOAT)*4);
                        tmp.d0 = d0;
                        citcoms_tracer_mixed_push(&CTM,&tmp);
                    }
                    undetermined_tracer++;
                    _ctd->extraq[j][i] = -2.0f;
                }

                all_tracers++;
            }
        }
    }

    fprintf(stdout,"%s: %d complex tracer in dump (%d/%d), %d tracer in dump\n",__func__,undetermined_tracer,CTM.len,CTM.len_alloc,
            all_tracers);

    Clock(0);
    qsort(CTM.data,CTM.len, sizeof(tracer_mixed), tracer_mixed_cmp);
    fprintf(stdout,"%s:use %f sec to sort complex tracers.\n",__func__, Clock(1));
    // align tracers flavor
    int undetermined_cells = 0, invalid_cells = 0;
    int cur_cell = 0;
    VTSDATAFLOAT sum_vof[4] = {0.0f};
    VTSDATAFLOAT avg_d0 = 0.0f;
    int pts_count = 0;

    // export sid
    double depth_th = (1.0 - 80.0/1750.0)*_cd->TransformR;
    int num_sid = 0;

    for(int k=0;k<=CTM.len;++k)
    {
        if(k == CTM.len ||tracer_mixed_cmp(CTM.data+k, CTM.data + cur_cell) != 0)
        {
            if(avg_d0 < depth_th*pts_count)
            {
                for(int j=1;j<=pts_count;++j){
                    CTM.data[k - j].sid = 1;
                    num_sid++;
                }
            }

            int Mid = VecMaxArgF(sum_vof, 4);
            if(Mid == 0)
            {
                invalid_cells ++;
                // may be void ?
            }

            // move volume to mat1,2
            float frac = 1.0f*pts_count/(1.0*pts_count - (sum_vof[0] + sum_vof[3]));
            sum_vof[1] *= frac;
            sum_vof[2] *= frac;

            // seek number of flavor 1,2
            int num_f1 = (int) roundf(sum_vof[1]);
            int num_f2 = (int) roundf(sum_vof[2]);

            qsort(CTM.data+k-pts_count,pts_count, sizeof(tracer_mixed), tracer_mixed_vofcmp);
            for(int j=1; j<=pts_count;++j)
            {
                if(j <= num_f1)
                    CTM.data[k - j].flavor = 1;
                else
                    CTM.data[k - j].flavor = 2;
                 // fprintf(stdout,"(%d) %f,%f,%.1f,%d,%d\n",j,CTM.data[k - j].vof[1],CTM.data[k - j].vof[2],CTM.data[k - j].flavor,num_f1,num_f2);

                // copy flavor to extraq
                int _p = CTM.data[k-j].id[0];
                int _j = CTM.data[k-j].id[1];
                int _i = CTM.data[k-j].id[2];
                assert(_cd->tracer[_p].extraq[_j][_i] < 0.);
                _cd->tracer[_p].extraq[_j][_i] = CTM.data[k - j].flavor;
            }

            // if(num_f2 >= 2 && num_f1 >= 2) exit(8);

            // reset counter vals

            if(k < CTM.len)
            {
                CTM.data[k].sid = 0;
                CTM.data[k].flavor = -1.0f;

                cur_cell = k;
                undetermined_cells ++;
                VecZeroF(sum_vof,4);
                pts_count = 0;
                avg_d0 = 0.f;
            }
        }

        if(k < CTM.len)
        {
            VecAddF(sum_vof, CTM.data[k].vof, 1.0f,4);
            avg_d0 += CTM.data[k].d0;
            pts_count++;
        }
    }

    fprintf(stdout,"%s: export %d complex tracers to vtp\n", __func__ ,num_sid);
    citcoms_tracer_mixed_export(&CTM,_cd,1,num_sid,"sid.vtp");
    fprintf(stdout,"%s: %d,complex tracers are divided into %d cells, %f per cells, find %d invalid cells\n",__func__,
            CTM.len,undetermined_cells, 1.0*undetermined_tracer/undetermined_cells, invalid_cells);
    citcoms_tracer_mixed_clean(&CTM);
    return undetermined_tracer;
}

int UpdateCitcomsDump(citcoms_dump * _cdp, SALEcData * _sdata, SALEcData * _rdata)
{
    UpdateCitcomsTempDump(_cdp,_sdata,_rdata);
    SphereIntegrateCitcomsDump(_cdp, "dump2");
    UpdateCitcomsTracerDump(_cdp,_sdata);
    CheckCitcomsTracerDump(_cdp);
    // citcoms_tracer_dump_pvtp(_cdp,"tracer_dump");
    return 0;
}

int SALEcGetCData(SALEcData * _sdata, int fId, VTSDATAFLOAT * _pos, VTSDATAFLOAT * _data)
{
    /*
     * get cell data at specified position
     * for fId
     */
    assert(fId >= 0 && fId < 99);
    VtsInfo * _vsf = _sdata->VSF;
    int NoC = _vsf->CellField[fId].NoC;

    VTSDATAFLOAT weight[3];
    int BlockId = BlockSearch(_sdata,_pos);
    int offset[3] = {0};
     if(BlockId < 0) return 1;

    OffsetSerchC(_sdata->VSF + BlockId, _pos, offset,weight);

    for(int kd=0;kd<VTSDIM;kd++)
    {
        if(offset[kd] < _sdata->Noffset)
        {
            weight[kd] = 1.0f;
        } else if(offset[kd] >= _sdata->VSF[BlockId].Nxp[kd] - _sdata->Noffset-2)
        {
            weight[kd] = 0.0f;
        }
    }

    VTSDATAFLOAT * tfield = _vsf[BlockId].CellField[fId].Data;
    unsigned long tnx = _vsf[BlockId].Nxp[0]-1;
    unsigned long tny = _vsf[BlockId].Nxp[1]-1;
    unsigned long tnz = _vsf[BlockId].Nxp[2]-1;
    unsigned long tx  = offset[0];
    unsigned long ty  = offset[1];
    unsigned long tz  = offset[2];

    unsigned long taId[2][2][2];
    for(int ia=0;ia<2;ia++)
    {
        for(int ja = 0; ja < 2; ++ja)
        {
            for(int ka=0;ka<2;ka++)
                taId[ia][ja][ka] = _lId3(tx+ia,ty+ja,tz+ka,tnx,tny,tnz);
        }
    }

    for(int kc=0;kc<NoC;kc++)
    {
        VTSDATAFLOAT ta[2][2];
        VTSDATAFLOAT tb[2];
        for(int ia=0;ia<2;ia++)
        {
            for(int ja = 0; ja < 2; ++ja)
            {
                ta[ia][ja] = lerp(tfield[_lId2(kc,taId[ia][ja][0],NoC,1)],
                                  tfield[_lId2(kc,taId[ia][ja][1],NoC,1)],
                                  weight[2]);
            }
            tb[ia] = lerp(ta[ia][0],ta[ia][1],weight[1]);
        }
        _data[kc] = lerp(tb[0],tb[1],weight[0]);
    }
    return 0;
}

int SALEcGetCDataN(SALEcData * _sdata, int *fId, int length, VTSDATAFLOAT * _pos, VTSDATAFLOAT * _data, int * mask)
{
    /*
     * get cell data at specified position
     * for fId
     */
    assert(NULL != _sdata && length >= 1);
    for(int k=0;k<length;++k)
    {
        assert(fId[k] >= 0 && fId[k] < 100);
    }
    VtsInfo * _vsf = _sdata->VSF;

    VTSDATAFLOAT weight[3];
    int BlockId = BlockSearch(_sdata,_pos);
    int offset[3] = {0};
    if(BlockId < 0) return -1;

    OffsetSerchC(_sdata->VSF + BlockId, _pos, offset,weight);

    if(NULL != mask)
    {
        mask[0] = BlockId;
        mask[1] = offset[0];
        mask[2] = offset[1];
        mask[3] = offset[2];
    }


    for(int kd=0;kd<VTSDIM;kd++)
    {
        if(offset[kd] < _sdata->Noffset)
        {
            weight[kd] = 1.0f;
        } else if(offset[kd] >= _sdata->VSF[BlockId].Nxp[kd] - _sdata->Noffset-2)
        {
            weight[kd] = 0.0f;
        }
    }

    unsigned long tnx = _vsf[BlockId].Nxp[0]-1;
    unsigned long tny = _vsf[BlockId].Nxp[1]-1;
    unsigned long tnz = _vsf[BlockId].Nxp[2]-1;
    unsigned long tx  = offset[0];
    unsigned long ty  = offset[1];
    unsigned long tz  = offset[2];

    unsigned long taId[2][2][2];
    for(int ia=0;ia<2;ia++)
    {
        for(int ja = 0; ja < 2; ++ja)
        {
            for(int ka=0;ka<2;ka++)
                taId[ia][ja][ka] = _lId3(tx+ia,ty+ja,tz+ka,tnx,tny,tnz);
        }
    }

    int di = 0;
    for(int fi=0;fi<length;++fi)
    {
        VTSDATAFLOAT * tfield = _vsf[BlockId].CellField[fId[fi]].Data;
        unsigned int NoC = _vsf[BlockId].CellField[fId[fi]].NoC;

        for(int kc=0;kc<NoC;kc++)
        {
            VTSDATAFLOAT ta[2][2];
            VTSDATAFLOAT tb[2];
            for(int ia=0;ia<2;ia++)
            {
                for(int ja = 0; ja < 2; ++ja)
                {
                    ta[ia][ja] = lerp(tfield[_lId2(kc,taId[ia][ja][0],NoC,1)],
                                      tfield[_lId2(kc,taId[ia][ja][1],NoC,1)],
                                      weight[2]);
                }
                tb[ia] = lerp(ta[ia][0],ta[ia][1],weight[1]);
            }
            _data[di] = lerp(tb[0],tb[1],weight[0]);
            di++;
        }
    }
    return di;
}

int SALEcGetCDataMask(SALEcData * _sdata, VTSDATAFLOAT * _pos, int * _id)
{
    /*
     * get cell data at specified position
     * for fId
     */

    VTSDATAFLOAT weight[3];
    int BlockId = BlockSearch(_sdata,_pos);
    int offset[3] = {0};
    if(BlockId < 0) return -1;
    OffsetSerchC(_sdata->VSF + BlockId, _pos, offset,weight);
    _id[0] = BlockId;
    _id[1] = offset[0];
    _id[2] = offset[1];
    _id[3] = offset[2];
    return 0;
}

void citcoms_tracer_dump_vtp(citcoms_tracer_dump * _ctd, const char * name)
{
    int num_tracers = 0;
    for(int j=0;j<_ctd->ncaps;++j)
    {
        num_tracers += _ctd->ntracers[j];
    }

    const int q_len = 5;
    float * tr_pos = malloc(sizeof(float)* (num_tracers+1) * 3);
    float * tr_mat = malloc(sizeof(float)*(num_tracers+1));
    float * tr_exq = malloc(sizeof(float)*(num_tracers+1) * q_len);
    assert(tr_pos != NULL);

    int tr_index = 0;
    int padded_tracer = 0;
    for(int j=0;j<_ctd->ncaps;++j)
    {
        for(int i=1;i<=_ctd->ntracers[j];++i)
        {
            tr_pos[tr_index*3 + 0] = (float) _ctd->basicq[j][i + 3*(_ctd->ntracers[j]+1)];
            tr_pos[tr_index*3 + 1] = (float) _ctd->basicq[j][i + 4*(_ctd->ntracers[j]+1)];
            tr_pos[tr_index*3 + 2] = (float) _ctd->basicq[j][i + 5*(_ctd->ntracers[j]+1)];
            tr_mat[tr_index] = (float) (_ctd->extraq[j][i]);
            tr_exq[tr_index*q_len + 0] = (float) _ctd->basicq[j][i + 0*(_ctd->ntracers[j]+1)];
            tr_exq[tr_index*q_len + 1] = (float) _ctd->basicq[j][i + 1*(_ctd->ntracers[j]+1)];
            tr_exq[tr_index*q_len + 2] = (float) _ctd->basicq[j][i + 2*(_ctd->ntracers[j]+1)];
            tr_index++;

            double dr0 = 0.0;
            for(int p=0;p<3;++p)
            {
                dr0 = tr_pos[tr_index*3 + p]*tr_pos[tr_index*3 + p];
            }
            if(dr0 > 1.0 + 1.0e-5)
            {
                fprintf(stdout, "Error[%d]:%f\n", padded_tracer++, dr0 - 1.0);
            }
        }
    }

    FILE * fp = fopen(name,"w");
    const char * vtp_data_format = "binary";
    vtp_file_header(fp,num_tracers);
    vtk_point_data_header(fp);
    vtk_dataarray_vec_f(fp,"id",vtp_data_format,tr_mat,num_tracers,1);
    vtk_dataarray_vec_f(fp,"q",vtp_data_format,tr_exq,num_tracers,5);
    vtk_point_data_trailer(fp);
    vtk_point_header(fp);
    vtk_dataarray_vec_f(fp,"coordinate",vtp_data_format,tr_pos,num_tracers,3);
    vtk_point_trailer(fp);
    vtp_file_trailer(fp);

    fclose(fp);
    free(tr_pos);
    free(tr_mat);
    free(tr_exq);

    fprintf(stdout,"%s:%s have %d tracers.\n",__func__,name,num_tracers);
}

void citcoms_tracer_dump_pvtp(citcoms_dump * _cdp, const char * name)
{
    const char header[] =
            "<?xml version=\"1.0\"?>\n"
            "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" compressor=\"vtkZLibDataCompressor\" byte_order=\"LittleEndian\">\n"
            "  <vtkMultiBlockDataSet>\n";
    char fname[4097];
    snprintf(fname,4096,"%s.full.vtm",name);
    FILE * fp = fopen(fname,"w");
    fputs(header, fp);
    for(int k=0;k<_cdp->nproc;++k)
    {
        char _tmp_name[4096];
        snprintf(_tmp_name,4095,"%s.%04d.vtp",name,k);
        citcoms_tracer_dump_vtp(_cdp->tracer + k,_tmp_name);
        fprintf(fp, "    <DataSet index=\"%d\" file=\"%s\"/>\n",k,_tmp_name);
    }
    fputs("  </vtkMultiBlockDataSet>\n",fp);
    fputs("</VTKFile>",fp);
    fclose(fp);
}

int WriteCitcomsDump(citcoms_dump * _cdp)
{
    assert(_cdp->nproc > 0);
    fprintf(stdout,"%s:write temperature dump in to %s\n",__func__,_cdp->temp_prefix);
    fprintf(stdout,"%s:write tracer dump in to %s\n",__func__,_cdp->tracer_prefix);

    Clock(0);
    // #pragma omp parallel for num_threads(32) shared(_cdp) default(none)
    for(int j=0;j<_cdp->nproc;++j)
    {
        char temp_dump_name[4097];
        char tracer_dump_name[4097];
        snprintf(temp_dump_name,4096,"%s.%d",_cdp->temp_prefix,j);
        snprintf(tracer_dump_name,4096,"%s.%d",_cdp->tracer_prefix,j);
        write_citcoms_temp_dump(_cdp->temp + j, temp_dump_name);
        write_citcoms_tracer_dump(_cdp->tracer + j,tracer_dump_name);
    }
    fprintf(stdout,"%s: use %f sec to write new dump file\n",__func__, Clock(1));
    return _cdp->nproc;
}

int write_citcoms_temp_dump(citcoms_temp_dump * _ctd, const char * fname)
{
    FILE * fp = fopen(fname,"wb");
    assert(fp != NULL);
    int header[6] = {0x789, _ctd->nox,_ctd->noy,_ctd->noz,_ctd->nno,_ctd->ncaps};
    fwrite(header, sizeof(int), 6, fp);
    assert(_ctd->ncaps > 0);
    for(int j=0;j<_ctd->ncaps;++j)
    {
        int dump_offset = j*_ctd->nno;
        fwrite(_ctd->x + dump_offset + 1, sizeof(double), _ctd->nno, fp);
        fwrite(_ctd->y + dump_offset + 1, sizeof(double), _ctd->nno, fp);
        fwrite(_ctd->z + dump_offset + 1, sizeof(double), _ctd->nno, fp);
        fwrite(_ctd->data + dump_offset, sizeof(double), _ctd->nno + 1, fp);
    }
    fclose(fp);
    return _ctd->nno;
}

int write_citcoms_tracer_dump(citcoms_tracer_dump * _ctd, const char * fname)
{
    FILE * fp = fopen(fname,"wb");
    assert(fp != NULL);

    int header[7] = {101,120,_ctd->ncaps,_ctd->num_basic_q,_ctd->num_extra_q,_ctd->nflavors,_ctd->itc};
    fwrite(header, sizeof(int), 7, fp);

    for(int j=0;j<_ctd->ncaps;++j)
    {
        fwrite(_ctd->ntracers + j, sizeof(int), 1, fp);
    }

    const size_t fix_num_basicq = 6;
    for(int j=0;j<_ctd->ncaps;++j)
    {
        fwrite(_ctd->basicq[j], sizeof(double), fix_num_basicq*(_ctd->ntracers[j]+1), fp);
        fwrite(_ctd->extraq[j], sizeof(double), _ctd->num_extra_q*(_ctd->ntracers[j]+1), fp);
        fwrite(_ctd->ielement[j], sizeof(int), _ctd->ntracers[j]+1, fp);
    }
    fclose(fp);
    return _ctd->ncaps;
}

int citcoms_tracer_mixed_init(citcoms_tracer_mixed * _ctm)
{
    _ctm->len_alloc = 4096;
    _ctm->len = 0;
    _ctm->data = malloc(sizeof(tracer_mixed)*_ctm->len_alloc);
    return _ctm->len_alloc;
}

int citcoms_tracer_mixed_push(citcoms_tracer_mixed * _ctm, tracer_mixed * x)
{
    if(_ctm->len + 2 > _ctm->len_alloc)
    {
        _ctm->data = realloc(_ctm->data,2*_ctm->len_alloc* sizeof(tracer_mixed));
        _ctm->len_alloc *= 2;
    }
    memcpy(_ctm->data + _ctm->len, x, sizeof(tracer_mixed));
    _ctm->len++;
    return _ctm->len;
}

void citcoms_tracer_mixed_clean(citcoms_tracer_mixed *_ctm)
{
    if(_ctm == NULL)
        return;
    free(_ctm->data);
    _ctm->len_alloc = 0;
    _ctm->len = 0;
}

int tracer_mixed_cmp(const void * _a, const void * _b)
{
    tracer_mixed * a = (tracer_mixed *) _a;
    tracer_mixed * b = (tracer_mixed *) _b;

    int mask_diff = 0;

    for(int k=0; k < 4; ++k)
    {
        if(a->mask[k] == b->mask[k])
            continue;

        mask_diff = a->mask[k] - b->mask[k];
        break;
    }
    return mask_diff;
}

int tracer_mixed_vofcmp(const void * _a, const void * _b)
{
    tracer_mixed * a = (tracer_mixed *) _a;
    tracer_mixed * b = (tracer_mixed *) _b;

    int mask_diff = 0;

    for(int k=1; k < 4; ++k)
    {
        if(fabsf(a->vof[k] - b->vof[k]) < 1.0e-6)
            continue;
        if(a->vof[k] < b->vof[k])
        {
            mask_diff = -1;
        }
        else
        {
            mask_diff = 1;
        }
        break;
    }
    return mask_diff;
}

void citcoms_tracer_mixed_export(citcoms_tracer_mixed * _ctm, citcoms_dump * _cd, int sid, int len,const char * name)
{
    assert(len >= 1);
    float * ex_pos = malloc(sizeof(float)*3*len);
    float * ex_vof = malloc(sizeof(float)*4*len);
    float * ex_d0_ = malloc(sizeof(float)*2*len);
    int ex_ind = 0;
    for(int k=0;k<_ctm->len;++k)
    {
        if(_ctm->data[k].sid != sid)
            continue;
        double * _pb =  ((_cd->tracer + _ctm->data[k].id[0])->basicq)[_ctm->data[k].id[1]];
        int nt = (_cd->tracer + _ctm->data[k].id[0])->ntracers[_ctm->data[k].id[1]];
        VTSDATAFLOAT _pos[3] = {
                (VTSDATAFLOAT)(_pb[_ctm->data[k].id[2] + 3*(nt + 1)]*_cd->TransformR),
                (VTSDATAFLOAT)(_pb[_ctm->data[k].id[2] + 4*(nt + 1)]*_cd->TransformR),
                (VTSDATAFLOAT)((_pb[_ctm->data[k].id[2] + 5*(nt + 1)] - 1.0)*_cd->TransformR)
        };

        memcpy(ex_pos + ex_ind*3, _pos, sizeof(float)*3);
        memcpy(ex_vof + ex_ind*4, _ctm->data[k].vof, sizeof(float)*3);
        ex_d0_[ex_ind*2 + 0] =  _cd->TransformR - _ctm->data[k].d0;
        ex_d0_[ex_ind*2 + 1] = _ctm->data[k].flavor;
        ex_ind++;
    }

    FILE * fp = fopen(name,"w");
    assert(fp != NULL);
    const char * vtp_data_format = "binary";
    vtp_file_header(fp,len);
    vtk_point_data_header(fp);
    vtk_dataarray_vec_f(fp,"vof",vtp_data_format,ex_vof,len,4);
    vtk_dataarray_vec_f(fp,"d",vtp_data_format,ex_d0_,len,2);
    vtk_point_data_trailer(fp);
    vtk_point_header(fp);
    vtk_dataarray_vec_f(fp,"coordinate",vtp_data_format,ex_pos,len,3);
    vtk_point_trailer(fp);
    vtp_file_trailer(fp);
    fclose(fp);

    free(ex_pos);
    free(ex_vof);
    free(ex_d0_);
}

int citcoms_offset(int i, int j, int k, int nx, int ny, int nz)
{
    assert(1 <= i && i <= nx);
    assert(1 <= j && j <= ny);
    assert(1 <= k && k <= nz);
    return k + nz*(i-1) + nz*nx*(j-1);
}


void tracer_finder_init(tracer_finder * _tf,citcoms_dump * _cd,int p[4])
{
    fprintf(stdout,"%s: checking %d->A,%d->B\n",__func__, p[0],p[1]);

    citcoms_temp_dump * _ctd = _cd->temp;
    _tf->nox = _ctd->nox;
    _tf->noy = _ctd->noy;
    _tf->noz = _ctd->noz;

    _tf->elx = _tf->nox - 1;
    _tf->ely = _tf->noy - 1;
    _tf->elz = _tf->noz - 1;

    double *A, *B, *C, *D;
    double S[12];

    int nodeA = citcoms_offset(1,1,1,_tf->nox,_tf->noy,_tf->noz);
    int nodeB = citcoms_offset(_tf->nox,1,1,_tf->nox,_tf->noy,_tf->noz);
    int nodeD = citcoms_offset(1,_tf->noy,1,_tf->nox,_tf->noy,_tf->noz);
    int nodeC = citcoms_offset(_tf->nox,_tf->noy,1,_tf->nox,_tf->noy,_tf->noz);

    int Snode[4] = {nodeA, nodeB, nodeC, nodeD};

    for(int j=0;j<4;++j)
    for(int k=0;k<3;++k)
    {
        S[j*3 + k] = (_ctd + p[j])->X[k][Snode[j]];
    }

    A = S + 0;
    B = S + 3;
    C = S + 6;
    D = S + 9;


    set_projection_axis(_tf->P2, _tf->P1, _tf->P0, A, B, C, D);
    set_projection_axis(_tf->Q2, _tf->Q1, _tf->Q0, D,A, B, C);

    VecD2F(_tf->P2f,_tf->P2,3);
    VecD2F(_tf->P1f,_tf->P1,3);
    VecD2F(_tf->P0f,_tf->P0,3);

    VecD2F(_tf->Q2f,_tf->Q2,3);
    VecD2F(_tf->Q1f,_tf->Q1,3);
    VecD2F(_tf->Q0f,_tf->Q0,3);


    // solve_local((double [3]){0.85979, 0.358, 0.319},_tf->va2,_tf->va1,_tf->va0);
    for(int s=0;s<4;++s)
    {
        _ctd = _cd->temp + p[s];
        for(int k=1; k <= _tf->noz; ++ k)
        {
            for(int i=1; i<= _tf->nox; ++i)
            {
                for(int j=1; j<= _tf->noy; ++j)
                {
                    int nodet = citcoms_offset(i,j,k,_tf->nox,_tf->noy,_tf->noz);
                    double xt[3] = {_ctd->x[nodet], _ctd->y[nodet], _ctd->z[nodet]};

                    // double cAngle = VecDot(xt,xs,3)/ VecLen(xt,3)/ VecLen(xs,3);
                    // double Angle = acos(cAngle)/M_PI * 180.0;
                    // double R = VecLen(xt, 3);
                    double lyt = solve_local(xt, _tf->P2, _tf->P1, _tf->P0);
                    double lxt = solve_local(xt, _tf->Q2, _tf->Q1, _tf->Q0);

                    float xtf[3];
                    VecD2F(xtf,xt,3);
                    float lytf = solve_local_f(xtf, _tf->P2f, _tf->P1f, _tf->P0f);
                    float lxtf = solve_local_f(xtf, _tf->Q2f, _tf->Q1f, _tf->Q0f);
                    if(k==1)
                    {
                        fprintf(stdout,"%d,%d,%d: => (D):%.4f,%.4f, (F): %.4f,%.4f\n",i,j,k,lxt*64+1,lyt*64+1, lxtf*64+1,lytf*64+1 );
                        fflush(stdout);
                    }
                }
            }
        }

    }

}


int citcoms_check_tracer_element(citcoms_dump * _cd)
{
    tracer_finder * ltf = malloc(_cd->nproc* sizeof(tracer_finder));
    assert(NULL != ltf);
    tracer_finder_init(ltf,_cd,(int[4]){1,3,7,5});

    free(ltf);
    return 0;
}


double solve_local(double * x,double * v2, double * v1, double * v0)
{
    double a = VecDot(v2,x,3);
    double b = VecDot(v1,x,3);
    double c = VecDot(v0,x,3);

    b = b/a * 0.5;
    c = c/a * 0.5;
    a = 0.5;

    double delta = b*b - 4.0*a*c;
    assert(delta >= 0.0);
    delta = sqrt(delta);

    // fprintf(stdout,"a=%f, b=%f, c=%f",a,b,c);
    double x1 = -b + delta;
    double x2 = -b - delta;

    if(x1 < 0.0 || x1 > 1.0) return x2;
    return x1;
}

float solve_local_f(float * x,float * v2, float * v1, float * v0)
{
    float tol = 1e-8;
    float a = VecDotF(v2,x,3);
    float b = VecDotF(v1,x,3);
    float c = VecDotF(v0,x,3);

    b = b/a * 0.5f;
    c = c/a * 0.5f;
    a = 0.5f;

    float delta = b*b - 4.0f*a*c;
    assert(delta >= 0.0);
    delta = sqrtf(delta);

    float x1 = -b + delta;
    float x2 = -b - delta;

    // assert(fabs(x1) < 1.0 + tol|| fabs(x2) < 1.0 + tol);
    if(fabs(x1) < 1.0 + tol && fabs(x2) < 1.0 + tol)
    {
        fprintf(stdout,"a=%f, b=%f, c=%f",a,b,c);
        exit(0);
    }

    if(fabsf(x1) > 1.0 + tol)
        return x2;
    else
        return x1;

    if(x1 < 0.0 || x1 > 1.0) return x2;
    return x1;
}


void set_projection_axis(double * n2, double * n1, double * n0, double * A, double * B, double *C, double *D)
{
    VecNormalize(A,3);
    VecNormalize(B,3);
    VecNormalize(C,3);
    VecNormalize(D,3);

    double DA[3], CB[3];
    VecLinear(DA,D, 1.0, A,-1.0,3);
    VecLinear(CB,C, 1.0, B,-1.0,3);

    double AxCB[3], DAxB[3];
    VecCross(AxCB,A,CB,3);
    VecCross(DAxB,DA,B,3);

    VecCross(n2, DA, CB, 3);
    VecLinear(n1, AxCB, 1.0, DAxB, 1.0, 3);
    VecCross(n0, A, B, 3);
}

CitcomsData * init_citcoms_data(const char * input)
{
    InputFile *cfp = OpenInputFile(input);
    char citcoms_input_path[4096];
    char citcoms_data_path[4096];
    char citcoms_output_path[4096];
    GetValueS(cfp,"citcoms.input",citcoms_input_path,"example");
    GetValueS(cfp,"citcoms.data",citcoms_data_path,"a");
    GetValueS(cfp,"citcoms.output",citcoms_output_path,"a");

    char StepOpt[4096];
    GetValueSk(cfp,"citcoms.step",StepOpt,0,"unknown");
    if(0!= strcasecmp("range", StepOpt))
    {
        fprintf(stdout,"%s: unknown step option:%s (shold be range)\n");
        exit(1);
    }

    InputFile * ifp = OpenInputFile(citcoms_input_path);
    CitcomsData * _cdata = malloc(sizeof(CitcomsData));

    _cdata->step0 = GetValueIk(cfp,"citcoms.step",1,"0");
    _cdata->step1 = GetValueIk(cfp,"citcoms.step",2,"1");
    _cdata->step_inc = GetValueIk(cfp,"citcoms.step",3,"1");
    _cdata->nprocx = GetValueI(ifp,"mesh.nprocx","-1");
    _cdata->nprocy = GetValueI(ifp,"mesh.nprocy","-1");
    _cdata->nprocz = GetValueI(ifp,"mesh.nprocz","-1");
    _cdata->nproc_surf = GetValueI(ifp,"mesh.nproc_surf","-1");
    _cdata->nproc = _cdata->nprocx*_cdata->nprocy*_cdata->nprocz*_cdata->nproc_surf;
    _cdata->nproc_surf *= _cdata->nprocx*_cdata->nprocy;

    int mgx = GetValueI(ifp,"mesh.mgunitx","-1");
    int mgy = GetValueI(ifp,"mesh.mgunity","-1");
    int mgz = GetValueI(ifp,"mesh.mgunitz","-1");

    int lev = GetValueI(ifp,"mesh.levels","-1");

    _cdata->nox = (int)(mgx * pow(2.0, lev-1)) + 1;
    _cdata->noy = (int)(mgy * pow(2.0, lev-1)) + 1;
    _cdata->noz = (int)(mgz * pow(2.0, lev-1)) + 1;

    _cdata->nno = _cdata->nox * _cdata->noy * _cdata->noz;
    _cdata->nel = (_cdata->nox-1)*(_cdata->noy-1)*(_cdata->noz-1);

    strcpy(_cdata->OutPrefix,citcoms_output_path);

    char datafile[4096];
    GetValueS(ifp,"mesh.datafile",datafile,"a");
    snprintf(_cdata->VtsPrefix,4096,"%s/%s",citcoms_data_path,datafile);
    CloseInputFile(cfp);
    CloseInputFile(ifp);

    _cdata->VSF = NULL;
    return _cdata;
}

int load_citcoms_step(CitcomsData * _cdata, int step)
{
    _cdata->VSF = malloc(sizeof(VtsInfo)*_cdata->nproc);
    #pragma omp parallel for num_threads(LOADTHREADS) default(shared)
    for(int k=0;k<_cdata->nproc;++k)
    {
        char vts_name[4096];
        snprintf(vts_name,4096,"%s.proc%d.%d.vts",_cdata->VtsPrefix,k,step);
        FILE * fp = fopen(vts_name,"r");
        VtsLoad(_cdata->VSF+k, fp);
        fclose(fp);
    }

    char vtm_name[4096];
    snprintf(vtm_name,4096,"%s.%d.vtm",_cdata->VtsPrefix,step);
    FILE * fp = fopen(vtm_name,"r");
    assert(fp!=NULL);
    char _buffer[4096];
    int k_attach = -1;
    _cdata->len_attach = 0;
    while(fgets(_buffer,4096,fp))
    {
        char * t = strstr(_buffer,"FieldData");
        if(k_attach < 0)
        {
            if(t == NULL)
                continue;
            else
                k_attach = 0;
        }

        strncpy(_cdata->attach[k_attach],_buffer, 4096);

        if(k_attach!=0 && t!=NULL)
        {
            _cdata->len_attach = k_attach+1;
            break;
        }
        else
        {
            k_attach++;
        }
    }

    fclose(fp);
    return _cdata->nproc;
}

int clean_citcoms_data(CitcomsData * _cdata)
{
    for(int k=0;k<_cdata->nproc;++k)
    {
        VtsInfoClean(_cdata->VSF + k);
    }
    free(_cdata->VSF);
    _cdata->VSF = NULL;
    return _cdata->nproc;
}

int close_citcoms_data(CitcomsData * _cdata)
{
    free(_cdata);
    return 0;
}
