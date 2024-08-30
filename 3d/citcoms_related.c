//
// Created by huachengli on 8/26/24.
//

#include "citcoms_related.h"
#include <omp.h>

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

    for(int j=0;j<_ctd->ncaps;++j)
    {
        int dump_offset = j*_ctd->nno;
        fread(_ctd->x + dump_offset + 1, sizeof(double), _ctd->nno, fp);
        fread(_ctd->y + dump_offset + 1, sizeof(double), _ctd->nno, fp);
        fread(_ctd->y + dump_offset + 1, sizeof(double), _ctd->nno, fp);
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

    for(int j=0;j<_ctd->ncaps;++j)
    {
        fread(_ctd->ntracers + j, sizeof(int), 1, fp);
        _ctd->basicq[j] = malloc(sizeof(double)*_ctd->num_basic_q*(_ctd->ntracers[j]+1));
        _ctd->extraq[j] = malloc(sizeof(double)*_ctd->num_extra_q*(_ctd->ntracers[j]+1));
        _ctd->ielement[j] = malloc(sizeof(int)*(_ctd->ntracers[j]+1));
    }

    for(int j=0;j<_ctd->ncaps;++j)
    {
        fread(_ctd->basicq[j], sizeof(double), _ctd->num_basic_q*(_ctd->ntracers[j]+1), fp);
        fread(_ctd->extraq[j], sizeof(double), _ctd->num_extra_q*(_ctd->ntracers[j]+1), fp);
        fread(_ctd->ielement[j], sizeof(int), _ctd->ntracers[j]+1, fp);
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
    _cd->TransformR = GetValueD(ifp,"citcoms.TransformR","1.0");
    GetValueS(ifp,"citcoms.temperature_dump",_cd->temp_prefix,"NONE");
    GetValueS(ifp,"citcoms.tracer_dump",_cd->tracer_prefix,"NONE");
    assert(_cd->nproc > 1);
    _cd->temp = malloc(sizeof(citcoms_temp_dump)* _cd->nproc);
    _cd->tracer = malloc(sizeof(citcoms_tracer_dump)* _cd->nproc);
    assert(_cd->temp != NULL && _cd->tracer != NULL);
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
    LoadVtsData(_sdata,_salec_data_path);
    return _sdata;
}

void CrCloseSALEcData(SALEcData * _sdata)
{
    if(NULL != _sdata)
    {
        CleanSALEcData(_sdata);
    }
    free(_sdata);
}

int UpdateCitcomsTempDump(citcoms_dump * _cd, SALEcData * _sdata)
{
    // locate tem field id
    VtsInfo * _vsf = _sdata->VSF;
    unsigned long fId= find_cellfield("Temperature",_vsf);
    unsigned long fId_vaccum = find_cellfield("VOF-0",_vsf);
    unsigned long NoC = _vsf->CellField[fId].NoC;
    assert(NoC == 1);

    int interpolation_pts = 0;
    int invalid_pts = 0;
    #pragma omp parallel for num_threads(32) default(shared)
    for(int p=0;p<_cd->nproc;++p)
    {
        citcoms_temp_dump * _ctd = _cd->temp + p;
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

                double d0 = sqrt(_pos[0]*_pos[0] + _pos[1]*_pos[1] + (_pos[2]+_cd->TransformR)*(_pos[2]+_cd->TransformR))/_cd->TransformR;
                if(d0 > 1.0 + 1.0e-5)
                {
                    #pragma omp critical
                    {
                        invalid_pts++;
                        fprintf(stdout,"%s:%d:invalid pts for temp field,%f,<%f,%f,%f>\n",__func__,
                                invalid_pts,d0-1.0,_pos[0],_pos[1],_pos[2]);
                    }
                }

                VTSDATAFLOAT * _data = malloc(sizeof(VTSDATAFLOAT)*(NoC+1));
                SALEcGetCDataN(_sdata,(int[]){fId_vaccum,fId},2,_pos,_data);
                double * _data_d = _ctd->data + dump_offset + i*NoC;
                for(int kc=0;kc<NoC;++kc) _data[kc] = (float) _data_d[kc];
            }
        }
    }
    return invalid_pts;
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

    int invalid_tracer = 0;

    #pragma omp parallel for num_threads(32) shared(_cd,IdList,_sdata,invalid_tracer,stdout) default(none)
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
                SALEcGetCDataN(_sdata, IdList, 4, _pos, vof2);

                VTSDATAFLOAT d0 = sqrt(_pos[0]*_pos[0] + _pos[1]*_pos[1] + (_pos[2]+_cd->TransformR)*(_pos[2]+_cd->TransformR));
                if(vof2[0] > 0.05)
                {
                    // padding data from internal
                    VTSDATAFLOAT pad_vec[3] = {_pos[0]/d0, _pos[1]/d0, (_pos[2]+_cd->TransformR)/d0};
                    VTSDATAFLOAT pad_step = _cd->TransformR / 2000.0;
                    for(int k=0;k<1500;++k)
                    {
                        for(int q=0;q<3;++q) _pos[q] = _pos[q] - pad_step*pad_vec[q];
                        SALEcGetCDataN(_sdata, IdList, 4, _pos, vof2);
                        if(vof2[0] <= 0.05){
                            break;
                        }
                    }
                }

                _ctd->extraq[j][i] = vof2[1] + vof2[2] + vof2[3];

                /*_ctd->basicq[j][i + 7*(_ctd->ntracers[j]+1)] = vof2[1];
                _ctd->basicq[j][i + 8*(_ctd->ntracers[j]+1)] = vof2[2];
                _ctd->basicq[j][i + 9*(_ctd->ntracers[j]+1)] = vof2[3];*/

                /*if(_ctd->basicq[j][i + 6*(_ctd->ntracers[j]+1)] > 15)
                {
                    #pragma omp critical
                    {
                        invalid_tracer++;
                        fprintf(stdout,"Error position[%d]:%f <%f,%f,%f>:<%f,%f,%f,%f>\n",
                                invalid_tracer,d0/_cd->TransformR - 1.0,
                                _pos[0],_pos[1],_pos[2] + _cd->TransformR,vof2[0],vof2[1],vof2[2],vof2[3]);
                    }
                }*/
            }
        }
    }

    return invalid_tracer;
}

int UpdateCitcomsDump(citcoms_dump * _cdp, SALEcData * _sdata)
{
    UpdateCitcomsTempDump(_cdp,_sdata);
    UpdateCitcomsTracerDump(_cdp,_sdata);
    citcoms_tracer_dump_pvtp(_cdp,"tracer_dump");
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

int SALEcGetCDataN(SALEcData * _sdata, int *fId, int length, VTSDATAFLOAT * _pos, VTSDATAFLOAT * _data)
{
    /*
     * get cell data at specified position
     * for fId
     */
    assert(length >= 1);
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
    int invalid_tracer = 0;
    for(int j=0;j<_ctd->ncaps;++j)
    {
        for(int i=1;i<=_ctd->ntracers[j];++i)
        {
            tr_pos[tr_index*3 + 0] = (float) _ctd->basicq[j][i + 3*(_ctd->ntracers[j]+1)];
            tr_pos[tr_index*3 + 1] = (float) _ctd->basicq[j][i + 4*(_ctd->ntracers[j]+1)];
            tr_pos[tr_index*3 + 2] = (float) _ctd->basicq[j][i + 5*(_ctd->ntracers[j]+1)];
            tr_mat[tr_index] = (float) (_ctd->extraq[j][i]);
            tr_exq[tr_index*q_len + 0] = (float) _ctd->basicq[j][i + 6*(_ctd->ntracers[j]+1)];
            tr_exq[tr_index*q_len + 1] = (float) _ctd->basicq[j][i + 7*(_ctd->ntracers[j]+1)];
            tr_exq[tr_index*q_len + 2] = (float) _ctd->basicq[j][i + 8*(_ctd->ntracers[j]+1)];
            tr_exq[tr_index*q_len + 3] = (float) _ctd->basicq[j][i + 9*(_ctd->ntracers[j]+1)];
            tr_exq[tr_index*q_len + 4] = (float) _ctd->basicq[j][i + 10*(_ctd->ntracers[j]+1)];
            tr_index++;

            double dr0 = 0.0;
            for(int p=0;p<3;++p)
            {
                dr0 = tr_pos[tr_index*3 + p]*tr_pos[tr_index*3 + p];
            }
            if(dr0 > 1.0 + 1.0e-5)
            {
                fprintf(stdout,"Error[%d]:%f\n",invalid_tracer++,dr0 - 1.0);
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
    for(int j=0;j<_cdp->nproc;++j)
    {
        char temp_dump_name[4097];
        char tracer_dump_name[4097];
        snprintf(temp_dump_name,4096,"%s.%d",_cdp->temp_prefix,j);
        snprintf(tracer_dump_name,4096,"%s.%d",_cdp->tracer_prefix,j);
        write_citcoms_temp_dump(_cdp->temp + j, temp_dump_name);
        write_citcoms_tracer_dump(_cdp->tracer + j,tracer_dump_name);
    }
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
        fwrite(_ctd->y + dump_offset + 1, sizeof(double), _ctd->nno, fp);
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

    for(int j=0;j<_ctd->ncaps;++j)
    {
        fwrite(_ctd->basicq[j], sizeof(double), _ctd->num_basic_q*(_ctd->ntracers[j]+1), fp);
        fwrite(_ctd->extraq[j], sizeof(double), _ctd->num_extra_q*(_ctd->ntracers[j]+1), fp);
        fwrite(_ctd->ielement[j], sizeof(int), _ctd->ntracers[j]+1, fp);
    }
    fclose(fp);
    return _ctd->ncaps;
}