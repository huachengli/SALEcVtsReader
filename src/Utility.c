//
// Created by huacheng on 10/19/21.
//

#include "Utility.h"
#include "InputParser.h"

void LoadInpInfo(SALEcData * _sdata,const char* _inputName)
{
    InputFile * ifp = OpenInputFile(_inputName);
#define SETNPGX(x,y) _sdata->Npgx[x]=GetValueI(ifp,"processor.npg"#y,"8")
#define SETNX(x,y) _sdata->Npx[x] = GetValueI(ifp,"mesh.np"#y,"32")
#define SETNX_NPGX(x,y) SETNX(x,y); SETNPGX(x,y);
    SETNX_NPGX(0,x)
    SETNX_NPGX(1,y)
    SETNX_NPGX(2,z)

    _sdata->Noffset = GetValueI(ifp,"processor.noffset","8");
    CloseInputFile(ifp);
    _sdata->VtsBlockNum = _sdata->Npgx[0]*_sdata->Npgx[1]*_sdata->Npgx[2];

}

void LoadVtsData(SALEcData * _sdata,const char * _vtsPrefix)
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

    SETGCL(0,1);
    SETGCL(1,_sdata->Npgx[0]);
    SETGCL(2,_sdata->Npgx[0]*_sdata->Npgx[1]);

}

void WriteGCL(SALEcData * _sdata,unsigned int _c,FILE *fp)
{
    for(int k=0;k<_sdata->Npgx[_c]*_sdata->Npx[_c]+1;++k)
    {
        fprintf(fp,"%.9e,",_sdata->GCLV[_c][k]);
    }
}

int cbsearch(VTSDATAFLOAT * _cl, VTSDATAFLOAT _x,unsigned long _len)
{
    if((_x < _cl[0]) || (_x > _cl[_len-1]))
        return -1;

    unsigned _ICL,_JCL,_MCL;
    _ICL = 0;
    _JCL = _len-1;

    do {
        _MCL = (_ICL + _JCL)/2;
        if(_x < _cl[_MCL])
            _JCL = _MCL;
        else
            _ICL = _MCL;
    } while (_ICL + 1< _JCL);

    return (int) _ICL;
}

#define GETNGCL_V(k) (_sdata->Npgx[k]*_sdata->Npx[k]+1)
#define GETNGCL_C(k) (_sdata->Npgx[k]*_sdata->Npx[k])
#define GETNPX_V(k) (_sdata->Npx[k]+1)
#define GETNPX_C(k) (_sdata->Npx[k])

#define SETIndexSearch(t) void IndexSearch##t(SALEcData * _sdata, VTSDATAFLOAT const * _x, int * _indices) \
{\
    for(int k=0;k<VTSDIM;++k) \
    { \
        _indices[k] = cbsearch(_sdata->GCL##t[k],_x[k],GETNGCL_##t(k)); \
    } \
}
#define SETGetBlockIndex(t) int GetBlockIndex##t(SALEcData * _sdata, const int * _indices) \
{ \
    if((-1==_indices[0])||(-1==_indices[1])||(-1==_indices[2])) return -1; \
    int _bindices[VTSDIM]; \
    for(int k=0;k<VTSDIM;k++) \
    { \
        _bindices[k] = _indices[k]/GETNPX_##t(k); \
    } \
    return _bindices[0] + _sdata->Npgx[0]*(_bindices[1] + _sdata->Npgx[1]*_bindices[2]); \
}

SETIndexSearch(C)
SETIndexSearch(V)
SETGetBlockIndex(V)
SETGetBlockIndex(C)



void SetMeshCLV(SALEcData * _sdata,Plane * _out,unsigned int px,unsigned py)
{
    _out->shape[0] = _sdata->nGCLV[px];
    _out->shape[1] = _sdata->nGCLV[py];

    _out->CL[0] = _sdata->GCLV[px];
    _out->CL[1] = _sdata->GCLV[py];

    _out->nCL[0] = _sdata->nGCLV[px];
    _out->nCL[1] = _sdata->nGCLV[py];
}

void SetMeshCLC(SALEcData * _sdata,Plane * _out,unsigned int px,unsigned py)
{
    _out->shape[0] = _sdata->nGCLC[px];
    _out->shape[1] = _sdata->nGCLC[py];

    _out->CL[0] = _sdata->GCLC[px];
    _out->CL[1] = _sdata->GCLC[py];

    _out->nCL[0] = _sdata->nGCLC[px];
    _out->nCL[1] = _sdata->nGCLC[py];
}

void SetPlaneMesh(SALEcData * _sdata, Plane * _out, void (*_set)(SALEcData * ,Plane * ,unsigned int ,unsigned int))
{
    // before this function Plane.n and Plane.d is known
    v_normalize(_out->n);
    VTSDATAFLOAT nx,ny,nz;
    unsigned int px,py,pz;
    if(fabs(_out->n[2]) > fabs(_out->n[1])) {
        nx = _out->n[1];
        px = 1;
        if (fabs(_out->n[2]) > fabs(_out->n[0])) {
            ny = _out->n[0];
            nz = _out->n[2];
            py = 0;
            pz = 2;
        } else {
            ny = _out->n[2];
            nz = _out->n[0];
            py = 2;
            pz = 0;
        }
    } else
    {
        nx = _out->n[2];
        px = 2;
        if(fabs(_out->n[1])> fabs(_out->n[0]))
        {
            ny = _out->n[0];
            nz = _out->n[1];
            py = 0;
            pz = 1;
        } else
        {
            ny = _out->n[1];
            nz = _out->n[0];
            py = 1;
            pz = 0;
        }
    }

    if(px==1 && py==0)
    {
        px = 0;
        py = 1;
        VTSDATAFLOAT tmp = nx;
        nx = ny;
        ny = tmp;
    }

    for(int k=0;k<3;++k)
    {
        _out->X0[k] = _sdata->GCLV[k][0];
    }

    _out->X0[pz] += -(_out->d + vdot3(_out->X0,_out->n))/nz;
    vzero(_out->nX[0]);
    vzero(_out->nX[1]);

    _out->nX[0][px] = 1.0;
    _out->nX[1][py] = 1.0;

    _out->nX[0][pz] += - nx/nz;
    _out->nX[1][pz] += - ny/nz;

    /*
    _out->shape[0] = _sdata->nGCLV[px];
    _out->shape[1] = _sdata->nGCLV[py];

    _out->CL[0] = _sdata->GCLV[px];
    _out->CL[1] = _sdata->GCLV[py];

    _out->nCL[0] = _sdata->nGCLV[px];
    _out->nCL[1] = _sdata->nGCLV[py];
     */
    _set(_sdata,_out,px,py);

}

int BlockSearch(SALEcData * _sdata,VTSDATAFLOAT * _x)
{
    /*
     * return the block index contain _x
     */
    int _indices[VTSDIM];
    for(int k=0;k<VTSDIM;++k)
    {
        _indices[k] = cbsearch(_sdata->BCLV[k],_x[k],_sdata->nBCLV[k]);
        if(_indices[k] <0) return -1;
    }

    return _indices[0] + _sdata->Npgx[0]*(_indices[1] + _sdata->Npgx[1]*_indices[2]);
}

int OffsetSerchV(VtsInfo * _vsf, VTSDATAFLOAT * _x, int * _indices, VTSDATAFLOAT * _weight)
{
    for(int k=0;k<VTSDIM;++k)
    {
        _indices[k] = cbsearch(_vsf->CLV[k],_x[k],_vsf->Nxp[k]);
        if(_indices[k] < 0)
        {
            fprintf(stdout,"out of index in offset search.\n");
            exit(0);
            return -1;
        }
        _weight[k] = (_x[k]-_vsf->CLV[k][_indices[k]])/(_vsf->CLV[k][_indices[k]+1]-_vsf->CLV[k][_indices[k]]);
    }
    _indices[VTSDIM] = _indices[0] + _vsf->Nxp[0]*(_indices[1] + _vsf->Nxp[1]*_indices[2]);

    return _indices[0] + _vsf->Nxp[0]*(_indices[1] + _vsf->Nxp[1]*_indices[2]);
}

int OffsetSerchC(VtsInfo * _vsf, VTSDATAFLOAT * _x, int * _indices, VTSDATAFLOAT * _weight)
{
    for(int k=0;k<VTSDIM;++k)
    {
        _indices[k] = cbsearch(_vsf->CLC[k],_x[k],_vsf->Nxp[k]-1);
        if(_indices[k] < 0)
        {
            fprintf(stdout,"out of index in offset search.\n");
            exit(0);
            return -1;
        }
        _weight[k] = (_x[k]-_vsf->CLC[k][_indices[k]])/(_vsf->CLC[k][_indices[k]+1]-_vsf->CLC[k][_indices[k]]);
    }
    _indices[VTSDIM] = _indices[0] + (_vsf->Nxp[0]-1)*(_indices[1] + (_vsf->Nxp[1]-1)*_indices[2]);

    return _indices[0] + (_vsf->Nxp[0]-1)*(_indices[1] + (_vsf->Nxp[1]-1)*_indices[2]);
}


void SetPlaneMask(SALEcData * _sdata, Plane * _out,int (*_search)(VtsInfo *, VTSDATAFLOAT *, int *,VTSDATAFLOAT*))
{
    _out->mask = (int **) malloc(sizeof(int*)*_out->nCL[0]);
    _out->offset = (int ***) malloc(sizeof(int**)*_out->nCL[0]);
    _out->weight = (VTSDATAFLOAT ***) malloc(sizeof(VTSDATAFLOAT**)*_out->nCL[0]);
    _out->coord = (VTSDATAFLOAT ***) malloc(sizeof(VTSDATAFLOAT**)*_out->nCL[0]);
    for(int k=0;k<_out->nCL[0];++k)
    {
        _out->mask[k] = (int*) malloc(sizeof(int)*_out->nCL[1]);
        _out->offset[k] = (int **) malloc(sizeof(int*)*_out->nCL[1]);
        _out->weight[k] = (VTSDATAFLOAT **) malloc(sizeof(VTSDATAFLOAT*)*_out->nCL[1]);
        _out->coord[k] = (VTSDATAFLOAT **) malloc(sizeof(VTSDATAFLOAT*)*_out->nCL[1]);
        for(int j=0;j<_out->nCL[1];++j)
        {
            _out->offset[k][j] = (int *) malloc(sizeof(int)*(VTSDIM+1));
            _out->weight[k][j] = (VTSDATAFLOAT *) malloc(sizeof(VTSDATAFLOAT)*(VTSDIM));
            _out->coord[k][j] = (VTSDATAFLOAT *) malloc(sizeof(VTSDATAFLOAT)*(VTSDIM));
        }
    }

    _out->scores = (unsigned long *) malloc(sizeof(unsigned long)*_sdata->VtsBlockNum);
    // record number of points need interpolated in every block
    for(int k=0;k<_sdata->VtsBlockNum;++k) _out->scores[k] = 0;

    for(int ix=0;ix<_out->nCL[0];ix++)
    {
        for(int jy=0;jy<_out->nCL[1];jy++)
        {
//            VTSDATAFLOAT ptmp[3];
            VTSDATAFLOAT * ptmp = _out->coord[ix][jy];
            v_copy(ptmp,_out->X0);
            v_add_liner(ptmp,_out->CL[0][ix]-_out->CL[0][0],_out->nX[0],_out->CL[1][jy]-_out->CL[1][0],_out->nX[1]);
            // set the coordinate of points on this plane

            _out->mask[ix][jy] = BlockSearch(_sdata,ptmp);
            if(_out->mask[ix][jy]>=0)
            {
                _out->scores[_out->mask[ix][jy]]++;
                _search(_sdata->VSF+_out->mask[ix][jy],ptmp,_out->offset[ix][jy],_out->weight[ix][jy]);
                // get the position of points plane in specific block (offset + weight)
            }
        }
    }

#if MASKLONG
    fprintf(stdout,"\n");
    for(int k=0;k<_sdata->VtsBlockNum;++k)
    {
        fprintf(stdout,"%5ld,",_out->scores[k]);
        if(k%12==11) fprintf(stdout,"\n");
    }
#endif
}

void GetPlaneDataC(SALEcData * _sdata, Plane * _out)
{
    VtsInfo * _vsf = _sdata->VSF;
    _out->Id = 100;
    for(int k=0;k<_vsf->CellNoF;++k)
    {
        if(0== strcasecmp(_out->Name,_vsf->CellField[k].Name))
        {
            _out->Id = k;
            break;
        }
    }
    if(_out->Id==100)
    {
        fprintf(stdout,"cannot find %s in cellfield\n",_out->Name);
        exit(0);
    }

    _out->NoC = _vsf->CellField[_out->Id].NoC;

    _out->data = (VTSDATAFLOAT ***) malloc(sizeof(VTSDATAFLOAT**)*_out->nCL[0]);
    for(unsigned long ix=0;ix<_out->nCL[0];ix++)
    {
        _out->data[ix] = (VTSDATAFLOAT **) malloc(sizeof(VTSDATAFLOAT*)*_out->nCL[1]);
        for(unsigned long jy=0;jy<_out->nCL[1];jy++)
        {
            _out->data[ix][jy] = (VTSDATAFLOAT*) malloc(sizeof(VTSDATAFLOAT)*(_out->NoC));
        }
    }

    for(unsigned long ix=0;ix<_out->nCL[0];ix++)
    {
        for(unsigned long jy=0;jy<_out->nCL[1];jy++)
        {
            if(_out->mask[ix][jy]<0)
            {
                continue;
            }
            VTSDATAFLOAT * tfield = _vsf[_out->mask[ix][jy]].CellField[_out->Id].Data;
            unsigned long tnx = _vsf[_out->mask[ix][jy]].Nxp[0]-1;
            unsigned long tny = _vsf[_out->mask[ix][jy]].Nxp[1]-1;
            unsigned long tnz = _vsf[_out->mask[ix][jy]].Nxp[2]-1;
            unsigned long tx  = _out->offset[ix][jy][0];
            unsigned long ty  = _out->offset[ix][jy][1];
            unsigned long tz  = _out->offset[ix][jy][2];

            // Apply nearest interpolation in some condition
            for(int kd=0;kd<VTSDIM;kd++)
            {
                if(_out->offset[ix][jy][kd] < _sdata->Noffset)
                {
                    _out->weight[ix][jy][kd] = 1.0;
                } else if(_out->offset[ix][jy][kd] >= _vsf[_out->mask[ix][jy]].Nxp[kd] - _sdata->Noffset-2)
                {
                    _out->weight[ix][jy][kd]=0.0;
                }
            }


            unsigned long taId[2][2][2];
            for(int ia=0;ia<2;ia++)
            {
                for(int ja = 0; ja < 2; ++ja)
                {
                    for(int ka=0;ka<2;ka++)
                        taId[ia][ja][ka] = _lId3(tx+ia,ty+ja,tz+ka,tnx,tny,tnz);
                }
            }

            for(int kc=0;kc<_out->NoC;kc++)
            {
                VTSDATAFLOAT ta[2][2];
                VTSDATAFLOAT tb[2];
                for(int ia=0;ia<2;ia++)
                {
                    for(int ja = 0; ja < 2; ++ja)
                    {
                        ta[ia][ja] = lerp(tfield[_lId2(kc,taId[ia][ja][0],_out->NoC,1)],
                                          tfield[_lId2(kc,taId[ia][ja][1],_out->NoC,1)],
                                          _out->weight[ix][jy][2]);
                    }
                    tb[ia] = lerp(ta[ia][0],ta[ia][1],_out->weight[ix][jy][1]);
                }
                _out->data[ix][jy][kc] = lerp(tb[0],tb[1],_out->weight[ix][jy][0]);
            }
        }
    }
}

void WritePlaneData(Plane * _out,int compoent,const char * fname)
{
    FILE * fp = fopen(fname,"w");
    for(int i=0;i<_out->nCL[0];i++)
    {
        for(int j=0;j<_out->nCL[1];j++)
        {
            fprintf(fp,"%10.6f,",_out->data[i][j][compoent]);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
}

void WritePlaneDataAll(Plane * _out,const char * fname)
{
    FILE * fp = fopen(fname,"w");
    for(int k=0;k<_out->NoC;k++)
    {
        for(int i=0;i<_out->nCL[0];i++)
        {
            for(int j=0;j<_out->nCL[1];j++)
            {
                fprintf(fp,"%10.6f,",_out->data[i][j][k]);
            }
            fprintf(fp,"\n");
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
}

void WritePlaneCoord(Plane * _out, const char * fname)
{
    FILE * fp = fopen(fname,"w");
    if(NULL==fp)
    {
        fprintf(stdout,"can not open %s\n",fname);
        exit(0);
    }
    fprintf(fp,"%ld,%ld,%d\n",_out->nCL[0],_out->nCL[1],VTSDIM);
    for(int i=0;i<_out->nCL[0];i++)
    {
        for(int j=0;j<_out->nCL[1];j++)
        {
            fprintf(fp,"%10.6f,%10.6f,%10.6f,",_out->coord[i][j][0],_out->coord[i][j][1],_out->coord[i][j][2]);
            fprintf(fp,"%d",_out->mask[i][j]);
            fprintf(fp,"\n");
        }
    }
    fclose(fp);
}

void MergeBlockC(SALEcData * _sdata,int * indices)
{
    // ... alternative method for nearest interpolation in some condition
    // in GetPlaneDataC
}

#define SAFEFREE(x) if(NULL!=(x)) free(x);
void CleanSALEcData(SALEcData * _sdata)
{

    for(int k=0;k<VTSDIM;++k)
    {
        SAFEFREE(_sdata->GCLV[k])
        SAFEFREE(_sdata->GCLC[k])
        SAFEFREE(_sdata->BCLV[k])
        SAFEFREE(_sdata->BCLC[k])
    }

    for(int k=0;k<_sdata->VtsBlockNum;++k)
    {
        VtsInfoClean(_sdata->VSF+k);
    }
    free(_sdata->VSF);
}

void CleanPlane(Plane * _out)
{
    for(int k=0;k<_out->nCL[0];++k)
    {

        for(int j=0;j<_out->nCL[1];++j)
        {
            SAFEFREE(_out->offset[k][j])
            SAFEFREE(_out->weight[k][j])
            SAFEFREE(_out->coord[k][j])
        }
        SAFEFREE(_out->mask[k])
        SAFEFREE(_out->offset[k])
        SAFEFREE(_out->weight[k])
        SAFEFREE(_out->coord[k])
    }
    SAFEFREE(_out->scores)
    SAFEFREE(_out->mask)
    SAFEFREE(_out->offset)
    SAFEFREE(_out->weight)
    SAFEFREE(_out->coord)
}


void CleanCache(ProfileCache * _cache,int _n)
{
    for(size_t j=0;j<_n;++j)
    {
        if(_cache[j].CacheState)
        {
            for(size_t ix=0;ix<_cache[j].nCL[0];++ix)
            {
                SAFEFREE(_cache[j].mask[ix])
                SAFEFREE(_cache[j].kz[ix])
                SAFEFREE(_cache[j].Lambda[ix])
            }
            SAFEFREE(_cache[j].mask)
            SAFEFREE(_cache[j].kz)
            SAFEFREE(_cache[j].Lambda)
        }
    }
}

#undef SAFEFREE

VTSDATAFLOAT* VtmGetCellData(SALEcData * _sdata, unsigned long k, unsigned long _i,unsigned long _j, unsigned long _k)
{
    unsigned long BlockId[VTSDIM];
    unsigned long BlockOffset[VTSDIM];

    BlockId[0] = (_i)/_sdata->Npx[0];
    BlockOffset[0] = 2+(_i)%_sdata->Npx[0];

    BlockId[1] = (_j)/_sdata->Npx[1];
    BlockOffset[1] = 2+(_j)%_sdata->Npx[1];

    BlockId[2] = (_k)/_sdata->Npx[2];
    BlockOffset[2] = 2+(_k)%_sdata->Npx[2];

//    unsigned long LId = BlockId[0] + _sdata->Npgx[0]*(BlockId[1] + _sdata->Npgx[1]*BlockId[2]);
    unsigned long LId = BlockId[0] + _sdata->Npgx[0]*(BlockId[1] + _sdata->Npgx[1]*BlockId[2]);
//    fprintf(stdout,"%d,%d,%d,%d\n",LId,BlockOffset[0],BlockOffset[1],BlockOffset[2]);
    return VtsGetCellData(_sdata->VSF+LId,k,BlockOffset[0],BlockOffset[1],BlockOffset[2]);
}

void GetProfileLim(SALEcData * _sdata, Plane * _out, VTSDATAFLOAT _tol)
{
    /*_out->CL[0] = _sdata->GCLC[0]; _out->nCL[0] = _sdata->nGCLC[0];
    _out->CL[1] = _sdata->GCLC[1]; _out->nCL[1] = _sdata->nGCLC[1];
    vzero(_out->nX[0]); _out->nX[0][0] = 1.;
    vzero(_out->nX[1]); _out->nX[1][1] = 1.;
    vzero(_out->n); _out->n[2] = 1.;*/

    VtsInfo * _vsf = _sdata->VSF;
    _out->Id = 100;
    for(int k=0;k<_vsf->CellNoF;++k)
    {
        if(0== strcasecmp(_out->Name,_vsf->CellField[k].Name))
        {
            _out->Id = k;
            break;
        }
    }
    if(_out->Id==100)
    {
        fprintf(stdout,"cannot find %s in cellfield\n",_out->Name);
        exit(0);
    }

    _out->NoC = _sdata->VSF->CellField[_out->Id].NoC;

    _out->data = (VTSDATAFLOAT ***) malloc(sizeof(VTSDATAFLOAT**)*_out->nCL[0]);
    for(unsigned long ix=0;ix<_out->nCL[0];ix++)
    {
        _out->data[ix] = (VTSDATAFLOAT **) malloc(sizeof(VTSDATAFLOAT*)*_out->nCL[1]);
        for(unsigned long jy=0;jy<_out->nCL[1];jy++)
        {
            _out->data[ix][jy] = (VTSDATAFLOAT*) malloc(sizeof(VTSDATAFLOAT)*(_out->NoC));
        }
    }

    unsigned long zColumn = _sdata->nGCLC[2];
    VTSDATAFLOAT * zColData = (VTSDATAFLOAT *) malloc(sizeof(VTSDATAFLOAT)*zColumn);
    for(unsigned long ix=0;ix<_out->nCL[0];ix++)
    {
        for(unsigned long jy=0;jy<_out->nCL[1];jy++)
        {
//            fprintf(stdout,"%d,%d,%d|",ix,jy,0);
            zColData[0] = VtmGetCellData(_sdata,_out->Id,ix,jy,0)[0];
            for(unsigned long kz=1;kz<zColumn;++kz)
            {
//                fprintf(stdout,"%d,%d,%d|",ix,jy,kz);
                zColData[kz] = VtmGetCellData(_sdata,_out->Id,ix,jy,kz)[0];
                if(zColData[kz]>_tol)
                {
                    VTSDATAFLOAT Lambda = (_tol-zColData[kz-1])/(zColData[kz] - zColData[kz-1]);
                    if(Lambda > 1.0) Lambda = 1.0;
                    if(Lambda < 0.0) Lambda = 0.0;
                    _out->data[ix][jy][0] = _sdata->GCLC[2][kz]*Lambda + _sdata->GCLC[2][kz-1]*(1.-Lambda);
                    break;
                }
            }
        }
    }

    free(zColData);
}

void GetProfileVOFWithCache__(SALEcData * _sdata, Plane * _out, ProfileCache * _cache, VTSDATAFLOAT _tol)
{

    // check the id of variable
    VtsInfo * _vsf = _sdata->VSF;
    // the max id used in vsf, 100 is much less than CellNoF
    _out->Id = 100;
    _out->VacuumId = 100;
    // search id of plane name
    _out->VacuumId = find_cellfield(_out->Vacuum,_vsf);
    _out->Id = find_cellfield(_out->Name,_vsf);

    if(_out->Id==100 || _out->VacuumId==100)
    {
        fprintf(stdout,"cannot find %s or %s(set as vacuum name) in cellfield\n",_out->Name,_out->Vacuum);
        exit(0);
    }


    // check the Profile cache, if a cache matched, do not loop over the entire mesh
    ProfileCache * _cur = _cache;
    for(int icahe=0;_cache[icahe].CacheState!=0;icahe++)
    {
        _cur = _cache + icahe;
#define dcmp(a,b) (fabsf((a)-(b))<1.0e-5)
        if(_cur->VacuumId == _out->VacuumId && 0==strcasecmp(_cur->Vacuum,_out->Vacuum) &&
                dcmp(_out->d,_cur->d) && dcmp(_out->n[0],_cur->n[0]) &&
                dcmp(_out->n[1],_cur->n[1]) && dcmp(_out->n[2],_cur->n[2]) &&
                dcmp(_cur->tol,_tol))
        {
            break;
        }
#undef dcmp
    }

    if(_cur->CacheState)
    {
        // get a cache
        GetProfileVOFReadCache(_sdata,_out,_cur);
    } else
    {
        // miss
        GetProfileVOFWriteCache(_sdata,_out,_cur,_tol);
    }
}

void GetProfileVOFReadCache(SALEcData * _sdata, Plane * _out,ProfileCache * _cache)
{
    _out->NoC = _sdata->VSF->CellField[_out->Id].NoC;

    // allocate data array for output
    _out->data = (VTSDATAFLOAT ***) malloc(sizeof(VTSDATAFLOAT**)*_out->nCL[0]);
    for(unsigned long ix=0;ix<_out->nCL[0];ix++)
    {
        _out->data[ix] = (VTSDATAFLOAT **) malloc(sizeof(VTSDATAFLOAT*)*_out->nCL[1]);
        for(unsigned long jy=0;jy<_out->nCL[1];jy++)
        {
            _out->data[ix][jy] = (VTSDATAFLOAT*) malloc(sizeof(VTSDATAFLOAT)*(_out->NoC));
        }
    }

    // using lambda,kz in cache
    for(unsigned long ix=0;ix<_out->nCL[0];ix++)
    {
        for(unsigned long jy=0;jy<_out->nCL[1];jy++)
        {

            unsigned long kz = _cache->kz[ix][jy];
            VTSDATAFLOAT Lambda = _cache->Lambda[ix][jy];
            VTSDATAFLOAT * vright = VtmGetCellData(_sdata,_out->Id,ix,jy,kz);
            VTSDATAFLOAT * vleft = VtmGetCellData(_sdata,_out->Id,ix,jy,kz-1);
            for(unsigned long lc=0;lc<_out->NoC;++lc)
                _out->data[ix][jy][lc] = (VTSDATAFLOAT) (1.-Lambda)*vleft[lc] + Lambda*vright[lc];
        }
    }
}

void GetProfileVOFWriteCache(SALEcData * _sdata, Plane * _out, ProfileCache * _cache, VTSDATAFLOAT _tol)
{
    /*
     * get the variable on upper surface
     */

    VtsInfo * _vsf = _sdata->VSF;
    _out->NoC = _sdata->VSF->CellField[_out->Id].NoC;

    // allocate data array for output
    _out->data = (VTSDATAFLOAT ***) malloc(sizeof(VTSDATAFLOAT**)*_out->nCL[0]);
    for(unsigned long ix=0;ix<_out->nCL[0];ix++)
    {
        _out->data[ix] = (VTSDATAFLOAT **) malloc(sizeof(VTSDATAFLOAT*)*_out->nCL[1]);
        for(unsigned long jy=0;jy<_out->nCL[1];jy++)
        {
            _out->data[ix][jy] = (VTSDATAFLOAT*) malloc(sizeof(VTSDATAFLOAT)*(_out->NoC));
        }
    }
    // allocate data array for cache
    _cache->nCL[0] = _out->nCL[0];
    _cache->nCL[1] = _out->nCL[1];
    _cache->kz = (int **) malloc(sizeof(int*)*_out->nCL[0]);
    _cache->Lambda = (VTSDATAFLOAT **) malloc(sizeof(VTSDATAFLOAT*)*_out->nCL[0]);
    _cache->mask = (int **) malloc(sizeof(int*)*_out->nCL[0]);
    for(unsigned long ix=0;ix<_out->nCL[0];++ix)
    {
        _cache->kz[ix] = (int *) malloc(sizeof(int)*_out->nCL[1]);
        _cache->Lambda[ix] = (VTSDATAFLOAT *) malloc(sizeof(VTSDATAFLOAT)*_out->nCL[1]);
        _cache->mask[ix] = (int *) malloc(sizeof(int)*_out->nCL[1]);
    }
    // copy some information into cache
    _cache->n[0] = _out->n[0];
    _cache->n[1] = _out->n[1];
    _cache->n[2] = _out->n[2];
    _cache->d    = _out->d;
    _cache->VacuumId = _out->VacuumId;
    strcpy(_cache->Vacuum,_out->Vacuum);
    _cache->tol = _tol;
    _cache->CacheState = 1;

    // search the result and write the result into cache
    unsigned long zColumn = _sdata->nGCLC[2];
    VTSDATAFLOAT * zColVacuum = (VTSDATAFLOAT *) malloc(sizeof(VTSDATAFLOAT)*zColumn);
    for(unsigned long ix=0;ix<_out->nCL[0];ix++)
    {
        for(unsigned long jy=0;jy<_out->nCL[1];jy++)
        {
            zColVacuum[0] = VtmGetCellData(_sdata,_out->VacuumId,ix,jy,0)[0];
            for(unsigned long kz=1;kz<zColumn;++kz)
            {

                zColVacuum[kz] = VtmGetCellData(_sdata,_out->VacuumId,ix,jy,kz)[0];

                // write information in cache
                _cache->mask[ix][jy] = _out->mask[ix][jy];
                _cache->kz[ix][jy] = (int) kz;
                if(zColVacuum[kz]>_tol)
                {
                    VTSDATAFLOAT Lambda = (_tol-zColVacuum[kz-1])/(zColVacuum[kz] - zColVacuum[kz-1]);
                    if(Lambda > 1.0) Lambda = 1.0;
                    if(Lambda < 0.0) Lambda = 0.0;
                    _cache->Lambda[ix][jy] = Lambda;

                    VTSDATAFLOAT * vright = VtmGetCellData(_sdata,_out->Id,ix,jy,kz);
                    VTSDATAFLOAT * vleft = VtmGetCellData(_sdata,_out->Id,ix,jy,kz-1);
                    for(unsigned long lc=0;lc<_out->NoC;++lc)
                        _out->data[ix][jy][lc] = (VTSDATAFLOAT) (1.-Lambda)*vleft[lc] + Lambda*vright[lc];

                    break;
                }
            }
        }
    }

    free(zColVacuum);
}

void GetProfileVOFLim(SALEcData * _sdata, Plane * _out, VTSDATAFLOAT _tol)
{
    /*
     * get the variable on upper surface
     */

    VtsInfo * _vsf = _sdata->VSF;

    // the max id used in vsf, 100 is much less than CellNoF
    _out->Id = 100;
    _out->VacuumId = 100;
    // search id of plane name
    for(int k=0;k<_vsf->CellNoF;++k)
    {
        if(0== strcasecmp(_out->Name,_vsf->CellField[k].Name))
        {
            _out->Id = k;
            break;
        }
    }
    if(_out->Id==100)
    {
        fprintf(stdout,"cannot find %s in cellfield\n",_out->Name);
        exit(0);
    }


    for(int k=0;k<_vsf->CellNoF;++k)
    {
        if(0== strcasecmp(_out->Vacuum,_vsf->CellField[k].Name))
        {
            _out->VacuumId = k;
            break;
        }
    }
    if(_out->VacuumId == 100)
    {
        fprintf(stdout,"cannot find %s(set as vacuum name) in cellfield",_out->Vacuum);
    }


    _out->NoC = _sdata->VSF->CellField[_out->Id].NoC;

    // allocate data array for output
    _out->data = (VTSDATAFLOAT ***) malloc(sizeof(VTSDATAFLOAT**)*_out->nCL[0]);
    for(unsigned long ix=0;ix<_out->nCL[0];ix++)
    {
        _out->data[ix] = (VTSDATAFLOAT **) malloc(sizeof(VTSDATAFLOAT*)*_out->nCL[1]);
        for(unsigned long jy=0;jy<_out->nCL[1];jy++)
        {
            _out->data[ix][jy] = (VTSDATAFLOAT*) malloc(sizeof(VTSDATAFLOAT)*(_out->NoC));
        }
    }

    unsigned long zColumn = _sdata->nGCLC[2];
    VTSDATAFLOAT * zColData = (VTSDATAFLOAT *) malloc(sizeof(VTSDATAFLOAT)*zColumn);
    VTSDATAFLOAT * zColVacuum = (VTSDATAFLOAT *) malloc(sizeof(VTSDATAFLOAT)*zColumn);
    for(unsigned long ix=0;ix<_out->nCL[0];ix++)
    {
        for(unsigned long jy=0;jy<_out->nCL[1];jy++)
        {
//            fprintf(stdout,"%d,%d,%d|",ix,jy,0);
            zColData[0]   = VtmGetCellData(_sdata,_out->Id,ix,jy,0)[0];
            zColVacuum[0] = VtmGetCellData(_sdata,_out->VacuumId,ix,jy,0)[0];
            for(unsigned long kz=1;kz<zColumn;++kz)
            {
//                fprintf(stdout,"%d,%d,%d|",ix,jy,kz);
                zColData[kz] = VtmGetCellData(_sdata,_out->Id,ix,jy,kz)[0];
                zColVacuum[kz] = VtmGetCellData(_sdata,_out->VacuumId,ix,jy,kz)[0];

                if(zColVacuum[kz]>_tol)
                {
                    VTSDATAFLOAT Lambda = (_tol-zColVacuum[kz-1])/(zColVacuum[kz] - zColVacuum[kz-1]);
                    if(Lambda > 1.0) Lambda = 1.0;
                    if(Lambda < 0.0) Lambda = 0.0;
                    // the position of surface can also be calculated here
                    // just a little performance loss compared to the file io
                    _out->data[ix][jy][0] = zColData[kz]*Lambda + zColData[kz-1]*(1.-Lambda);
                    break;
                }
            }
        }
    }

    free(zColData);
    free(zColVacuum);
}


/*
 * compare tgt and src , ignoring case
 * if tgt is a prefix of src, return the occurrence position
 * else return NULL
 */
char * head__strcasestr(const char * src, const char *tgt)
{
    if(NULL==src || NULL==tgt) return NULL;
    char * s1 =  (char *) src;
    char * s2 =  (char *) tgt;

    while(*s1 && *s2 && !(tolower(*s1) - tolower(*s2)))
    {
        s1++;s2++;
    }

    if(!(*s2))
        return s1;

    return NULL;
}

unsigned long find_cellfield(const char _src[], VtsInfo * _vsf)
{
    unsigned long result= 100;
    for(int k=0;k<_vsf->CellNoF;++k)
    {
        if(0== strcasecmp(_src,_vsf->CellField[k].Name))
        {
            result = k;
            break;
        }
    }
    return result;
}


void GetProfileWriteCache(SALEcData * _sdata, Plane * _out, ProfileCache * _cache, VTSDATAFLOAT _tol)
{
    /*
     * get the height on upper surface
     * most of this function is same as GetProfileVOFWriteCache
     */

    VtsInfo * _vsf = _sdata->VSF;
    _out->NoC = _sdata->VSF->CellField[_out->Id].NoC;

    // allocate data array for output
    _out->data = (VTSDATAFLOAT ***) malloc(sizeof(VTSDATAFLOAT**)*_out->nCL[0]);
    for(unsigned long ix=0;ix<_out->nCL[0];ix++)
    {
        _out->data[ix] = (VTSDATAFLOAT **) malloc(sizeof(VTSDATAFLOAT*)*_out->nCL[1]);
        for(unsigned long jy=0;jy<_out->nCL[1];jy++)
        {
            _out->data[ix][jy] = (VTSDATAFLOAT*) malloc(sizeof(VTSDATAFLOAT)*(_out->NoC));
        }
    }
    // allocate data array for cache
    _cache->nCL[0] = _out->nCL[0];
    _cache->nCL[1] = _out->nCL[1];
    _cache->kz = (int **) malloc(sizeof(int*)*_out->nCL[0]);
    _cache->Lambda = (VTSDATAFLOAT **) malloc(sizeof(VTSDATAFLOAT*)*_out->nCL[0]);
    _cache->mask = (int **) malloc(sizeof(int*)*_out->nCL[0]);
    for(unsigned long ix=0;ix<_out->nCL[0];++ix)
    {
        _cache->kz[ix] = (int *) malloc(sizeof(int)*_out->nCL[1]);
        _cache->Lambda[ix] = (VTSDATAFLOAT *) malloc(sizeof(VTSDATAFLOAT)*_out->nCL[1]);
        _cache->mask[ix] = (int *) malloc(sizeof(int)*_out->nCL[1]);
    }
    // copy some information into cache
    _cache->n[0] = _out->n[0];
    _cache->n[1] = _out->n[1];
    _cache->n[2] = _out->n[2];
    _cache->d    = _out->d;
    _cache->VacuumId = _out->VacuumId;
    strcpy(_cache->Vacuum,_out->Vacuum);
    _cache->tol = _tol;
    _cache->CacheState = 1;

    // search the result and write the result into cache
    unsigned long zColumn = _sdata->nGCLC[2];
    VTSDATAFLOAT * zColVacuum = (VTSDATAFLOAT *) malloc(sizeof(VTSDATAFLOAT)*zColumn);
    for(unsigned long ix=0;ix<_out->nCL[0];ix++)
    {
        for(unsigned long jy=0;jy<_out->nCL[1];jy++)
        {
            zColVacuum[0] = VtmGetCellData(_sdata,_out->VacuumId,ix,jy,0)[0];
            for(unsigned long kz=1;kz<zColumn;++kz)
            {

                zColVacuum[kz] = VtmGetCellData(_sdata,_out->VacuumId,ix,jy,kz)[0];

                // write information in cache
                _cache->mask[ix][jy] = _out->mask[ix][jy];
                _cache->kz[ix][jy] = (int) kz;
                if(zColVacuum[kz]>_tol)
                {
                    VTSDATAFLOAT Lambda = (_tol-zColVacuum[kz-1])/(zColVacuum[kz] - zColVacuum[kz-1]);
                    if(Lambda > 1.0) Lambda = 1.0;
                    if(Lambda < 0.0) Lambda = 0.0;
                    _cache->Lambda[ix][jy] = Lambda;

                    _out->data[ix][jy][0] = _sdata->GCLC[2][kz]*Lambda + _sdata->GCLC[2][kz-1]*(1.-Lambda);

                    break;
                }
            }
        }
    }

    free(zColVacuum);
}

void GetProfileReadCache(SALEcData * _sdata, Plane * _out,ProfileCache * _cache)
{
    /*
     * Get the height of upper surface
     * use cache
     */
    _out->NoC = _sdata->VSF->CellField[_out->Id].NoC;

    // allocate data array for output
    _out->data = (VTSDATAFLOAT ***) malloc(sizeof(VTSDATAFLOAT**)*_out->nCL[0]);
    for(unsigned long ix=0;ix<_out->nCL[0];ix++)
    {
        _out->data[ix] = (VTSDATAFLOAT **) malloc(sizeof(VTSDATAFLOAT*)*_out->nCL[1]);
        for(unsigned long jy=0;jy<_out->nCL[1];jy++)
        {
            _out->data[ix][jy] = (VTSDATAFLOAT*) malloc(sizeof(VTSDATAFLOAT)*(_out->NoC));
        }
    }

    // using lambda,kz in cache
    for(unsigned long ix=0;ix<_out->nCL[0];ix++)
    {
        for(unsigned long jy=0;jy<_out->nCL[1];jy++)
        {

            unsigned long kz = _cache->kz[ix][jy];
            VTSDATAFLOAT Lambda = _cache->Lambda[ix][jy];
            _out->data[ix][jy][0] = _sdata->GCLC[2][kz]*Lambda + _sdata->GCLC[2][kz-1]*(1.-Lambda);
        }
    }
}

void GetProfileWithCache__(SALEcData * _sdata, Plane * _out, ProfileCache * _cache, VTSDATAFLOAT _tol)
{

    // check the id of variable
    VtsInfo * _vsf = _sdata->VSF;
    // the max id used in vsf, 100 is much less than CellNoF
    _out->Id = 100;
    _out->VacuumId = 100;
    // search id of plane name
    _out->VacuumId = find_cellfield(_out->Vacuum,_vsf);
    _out->Id = find_cellfield(_out->Name,_vsf);
    if(_out->Id==100 || _out->VacuumId==100)
    {
        fprintf(stdout,"cannot find %s in cellfield\n",_out->Name);
        exit(0);
    }


    // check the Profile cache, if a cache matched, do not loop over the entire mesh
    ProfileCache * _cur = _cache;
    for(int icahe=0;_cache[icahe].CacheState!=0;icahe++)
    {
        _cur = _cache + icahe;
#define dcmp(a,b) (fabsf((a)-(b))<1.0e-5)
        if(_cur->VacuumId == _out->VacuumId && 0==strcasecmp(_cur->Vacuum,_out->Vacuum) &&
           dcmp(_out->d,_cur->d) && dcmp(_out->n[0],_cur->n[0]) &&
           dcmp(_out->n[1],_cur->n[1]) && dcmp(_out->n[2],_cur->n[2]) &&
           dcmp(_cur->tol,_tol))
        {
            break;
        }
#undef dcmp
    }

    if(_cur->CacheState)
    {
        // get a cache
        GetProfileReadCache(_sdata,_out,_cur);
    } else
    {
        // miss
        GetProfileWriteCache(_sdata,_out,_cur,_tol);
    }
}


