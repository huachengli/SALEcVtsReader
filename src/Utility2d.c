//
// Created by huach on 1/16/2024.
//

#include "Utility2d.h"


void Load2dInpInfo(SALEcData * _sdata,SALEcPlanetInfo * _pdata,const char* _inputName)
{
    InputFile * ifp = OpenInputFile(_inputName);
    _sdata->Npgx[0] = GetValueI(ifp,"processor.npgx","2");
    _sdata->Npgx[1] = GetValueI(ifp,"processor.npgy","2");
    _sdata->Npgx[2] = 1;
    _sdata->Npx[0] = GetValueI(ifp,"mesh.npx","32");
    _sdata->Npx[1] = GetValueI(ifp,"mesh.npy","32");
    _sdata->Npx[2] = 1;
    _sdata->Noffset = GetValueI(ifp,"processor.noffset","2");

    // load planet info
    char TargetType[256];
    GetValueS(ifp,"target.type",TargetType,"none");
    if(0!= strcasecmp(TargetType,"sphere"))
    {
        fprintf(stdout,"%s:target type = %s\n",__func__ ,TargetType);
        exit(0);
    }
    _pdata->r_planet = GetValueD(ifp,"taregt.r_planet","1750.0e3");
    _pdata->center[0] = 0.;
    _pdata->center[1] = - _pdata->r_planet;
    _pdata->center[2] = 0.;
    GetValueS(ifp,"output.prefix",_pdata->prefix,"ParaTest");
    CloseInputFile(ifp);
    _sdata->VtsBlockNum = _sdata->Npgx[0]*_sdata->Npgx[1]*_sdata->Npgx[2];
}

void Reset2dZcoord(SALEcData * _sdata)
{
    _sdata->GCLC[2][0] = 0.f;
    _sdata->GCLV[2][1] = 1.f;
    _sdata->GCLV[2][1] =-1.f;
}

void LoadCitcom2dXX(Citcom2dXX * _cdata, const char * _cname)
{
    FILE * fp = fopen(_cname,"rb");
    if(NULL == fp)
    {
        fprintf(stdout,"%s:cannot open %s\n",__func__,_cname);
        exit(0);
    }
    fread(&(_cdata->nno), sizeof(int),1,fp);
    _cdata->XXs = (float *) malloc(sizeof(float)*_cdata->nno*3);
    _cdata->XX  = (float *) malloc(sizeof(float)*_cdata->nno*3);
    _cdata->T   = (float *) malloc(sizeof(float)*_cdata->nno);

    fread(_cdata->XXs, sizeof(float),3*_cdata->nno,fp);
    fclose(fp);

    char * _bname = strtok(_cname,".");
    strcpy(_cdata->cname,_bname);
}

void ConsistentCitcom2dXX(Citcom2dXX * _cdata, SALEcPlanetInfo * _sinfo)
{
    // convert coordinate in XXs (cartesian 2d nondimensional)
    // to XX (cartesian 2d aligned with SALEc)
    for(int k=0;k<_cdata->nno;++k)
    {
        _cdata->XX[k*3 + 0] = _cdata->XXs[k*3 + 0] * _sinfo->r_planet;
        _cdata->XX[k*3 + 1] = _cdata->XXs[k*3 + 1] * _sinfo->r_planet - _sinfo->r_planet;
        _cdata->XX[k*3 + 2] = _cdata->XXs[k*3 + 2] * _sinfo->r_planet;
        assert(_cdata->XX[k*3 + 2] == 0.);
    }
}


void ExportCitcom2dT(Citcom2dXX * _cdata)
{
    char tname[256];
    strcpy(tname,_cdata->cname);
    strcat(tname,".T");
    FILE * fp = fopen(tname,"wb");
    if(NULL == fp)
    {
        fprintf(stdout,"%s:cannot open %s\n",__func__,tname);
    }
    fwrite(&(_cdata->nno), sizeof(int),1,fp);
    fwrite(_cdata->T, sizeof(float),_cdata->nno,fp);
    fclose(fp);
}

int OffsetSerchC2d(VtsInfo * _vsf, VTSDATAFLOAT * _x, int * _indices, VTSDATAFLOAT * _weight)
{
    /*
     * get the weight and offset in the global mesh, with specified _x (coordinate)
     * weight is the local coordinate in the searched cells,
     */
    for(int k=0;k<2;++k)
    {
        _indices[k] = cbsearch(_vsf->CLC[k],_x[k],_vsf->Nxp[k]-1);
        assert(_indices[k] >= 0); // out of index in offset search; conflict with the BlockId searching
        _weight[k] = (_x[k]-_vsf->CLC[k][_indices[k]])/(_vsf->CLC[k][_indices[k]+1]-_vsf->CLC[k][_indices[k]]);
    }
    _indices[2] = 0; _weight[2] = 0.0; // the z direction
    // the global offset is stored in _indices[VTSDIM]
    _indices[VTSDIM] = _indices[0] + (_vsf->Nxp[0]-1)*(_indices[1] + (_vsf->Nxp[1]-1)*_indices[2]);
    return _indices[0] + (_vsf->Nxp[0]-1)*(_indices[1] + (_vsf->Nxp[1]-1)*_indices[2]);
}

int BlockSearch2d(SALEcData * _sdata,VTSDATAFLOAT * _x)
{
    /*
     * return the block index contain _x
     */
    int _indices[2] = {-1};
    for(int k=0;k<2;++k)
    {
        _indices[k] = cbsearch(_sdata->BCLV[k],_x[k],_sdata->nBCLV[k]);
        if(_indices[k] <0) return -1;
    }

    return _indices[1] + _sdata->Npgx[1]*_indices[0];
}



void InterpolateCitcom2dT(Citcom2dXX * _cdata, SALEcData * _sdata)
{
    // find the vsf id of temperature fiedl
    VtsInfo * _vsf = _sdata->VSF;
    int TemId = -1;
    const char * TemVsfName = "e_tem";
    for(int k=0;k<_vsf->CellNoF;++k)
    {
        if(0== strcasecmp(TemVsfName,_vsf->CellField[k].Name))
        {
            TemId = k;
            break;
        }
    }

    if(TemId==100)
    {
        fprintf(stdout,"%s:cannot find %s in cellfield\n",__func__ ,TemVsfName);
        exit(0);
    }

    // check the Noc of e_tem
    int TemNoc = _vsf[0].CellField[TemId].NoC;
    assert(TemNoc == 1);

    // loop on all the citcom2d coordinates
    for(int k=0;k<_cdata->nno;++k)
    {
        VTSDATAFLOAT * xk = _cdata->XX + 3*k;
        VTSDATAFLOAT local_xk[4] = {0.};
        int indices[4] = {0};

        int BlockId = BlockSearch2d(_sdata,xk);
        assert(BlockId >= 0); // invalid block id
        OffsetSerchC2d(_sdata->VSF + BlockId,xk,indices,local_xk);

        // some grids on the boundary of blocks using nearest interpolation
        for(int j=0;j<2;++j)
        {
            if(indices[j] < _sdata->Noffset)
            {
                local_xk[j] = 1.0;
            } else if(indices[j] >= _vsf[BlockId].Nxp[j] - _sdata->Noffset-2)
            {
                local_xk[j] = 0.0;
            }
        }

        // data id around this point
        int taId[2][2];
        for(int ia=0;ia<2;ia++)
        {
            for(int ja = 0; ja < 2; ++ja)
            {
                taId[ia][ja] = _lId3(indices[0]+ia,indices[1]+ja,0,_vsf[BlockId].Nxp[0]-1,_vsf[BlockId].Nxp[1]-1,1);
            }
        }

        VTSDATAFLOAT t0 = lerp(_vsf[BlockId].CellField[TemId].Data[taId[0][0]],_vsf[BlockId].CellField[TemId].Data[taId[1][0]],local_xk[0]);
        VTSDATAFLOAT t1 = lerp(_vsf[BlockId].CellField[TemId].Data[taId[0][1]],_vsf[BlockId].CellField[TemId].Data[taId[1][1]],local_xk[0]);
        _cdata->T[k] = lerp(t0,t1,local_xk[1]); // the x direction has been skipped.
    }
}

void LoadVts2dData(SALEcData * _sdata,const char * _vtsPrefix)
{
    _sdata->VSF = malloc(sizeof(VtsInfo)*_sdata->VtsBlockNum);
    for(int k=0;k<VTSDIM;++k)
    {
        _sdata->GCLV[k] = (VTSDATAFLOAT*) malloc(sizeof(VTSDATAFLOAT)*(_sdata->Npgx[k]*_sdata->Npx[k]+1));
        _sdata->GCLC[k] = (VTSDATAFLOAT*) malloc(sizeof(VTSDATAFLOAT)*(_sdata->Npgx[k]*_sdata->Npx[k]));
        _sdata->nGCLV[k] = _sdata->Npgx[k]*_sdata->Npx[k]+1;
        _sdata->nGCLC[k] = _sdata->Npgx[k]*_sdata->Npx[k];

        _sdata->BCLV[k] = (VTSDATAFLOAT*) malloc(sizeof(VTSDATAFLOAT)*(_sdata->Npgx[k]+1));
        _sdata->BCLC[k] = NULL;
        _sdata->nBCLV[k] = _sdata->Npgx[k]+1;
    }

    strcpy(_sdata->VtsPrefix,_vtsPrefix);
    int TaskFinished = 0;
#pragma omp parallel for num_threads(12) shared(_sdata,stdout,TaskFinished) default(none)
    for(int k=0;k<_sdata->VtsBlockNum;++k)
    {
        char VtsName[200];
        sprintf(VtsName,_sdata->VtsPrefix,k);
        FILE * fp = fopen(VtsName,"r");
        if(NULL==fp)
        {
            fprintf(stdout,"cannot open %s\n",VtsName);exit(0);
        }
        VtsLoad(_sdata->VSF+k,fp);

#pragma omp critical
        {
            TaskFinished ++;
            if(TaskFinished%10==9)
            {
                fprintf(stdout,"#");
                fflush(stdout);
            }
        };
        fclose(fp);
    }

#define SETGCLV(dim,y) for(int k=0;k<_sdata->Npgx[dim];++k) \
    {\
    memcpy(_sdata->GCLV[dim]+k*_sdata->Npx[dim],_sdata->VSF[k*(y)].CLV[dim]+_sdata->Noffset, sizeof(VTSDATAFLOAT)*(_sdata->Npx[dim]+1));\
    }
#define SETGCLC(dim) for(int k=0;k<_sdata->Npgx[dim]*_sdata->Npx[dim];++k) \
    { \
    _sdata->GCLC[dim][k] = 0.5*(_sdata->GCLV[dim][k]+_sdata->GCLV[dim][k+1]); \
    }
#define SETBCLV(dim) for(int k=0;k<_sdata->Npgx[dim]+1;++k) \
    {\
    _sdata->BCLV[dim][k] = _sdata->GCLV[dim][k*_sdata->Npx[dim]]; \
    }

#define SETGCL(dim,y) SETGCLV(dim,y) SETGCLC(dim) SETBCLV(dim)

    SETGCL(0,_sdata->Npgx[1]); // the sequence of 2d is different !! colum first rank.
    SETGCL(1,1);

    // set the z direction coordinate line value
    _sdata->BCLV[2][0] = -1.;
    _sdata->BCLV[2][1] =  1.;

    // notice the BCLC is not used, Block Id is determined according to the BCLV
}