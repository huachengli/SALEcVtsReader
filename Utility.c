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
    _sdata->VSF = malloc(sizeof(VtsInfo)*_sdata->VtsBlockNum);

    /*_sdata->GCLV = (VTSDATAFLOAT **) malloc(sizeof(VTSDATAFLOAT*)*VTSDIM);
    _sdata->GCLC = (VTSDATAFLOAT **) malloc(sizeof(VTSDATAFLOAT*)*VTSDIM);
    _sdata->BCLV = (VTSDATAFLOAT **) malloc(sizeof(VTSDATAFLOAT*)*VTSDIM);
    _sdata->BCLC = (VTSDATAFLOAT **) malloc(sizeof(VTSDATAFLOAT*)*VTSDIM);*/
    for(int k=0;k<VTSDIM;++k)
    {
        _sdata->GCLV[k] = (VTSDATAFLOAT*) malloc(sizeof(VTSDATAFLOAT)*(_sdata->Npgx[k]*_sdata->Npx[k]+1));
        _sdata->GCLC[k] = (VTSDATAFLOAT*) malloc(sizeof(VTSDATAFLOAT)*(_sdata->Npgx[k]*_sdata->Npx[k]));
        _sdata->nGCLV[k] = _sdata->Npgx[k]*_sdata->Npx[k]+1;
        _sdata->nGCLC[k] = _sdata->Npgx[k]*_sdata->Npx[k];

        _sdata->BCLV[k] = (VTSDATAFLOAT*) malloc(sizeof(VTSDATAFLOAT)*(_sdata->Npgx[k]+1));
        _sdata->nBCLV[k] = _sdata->Npgx[k]+1;
    }
}

void LoadVtsData(SALEcData * _sdata,const char * _vtsPrefix)
{
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
    _sdata->GCLC[dim][k] = 0.5*(_sdata->GCLV[dim][k]+_sdata->GCLV[dim][k]); \
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

void CleanSALEcData(SALEcData * _sdata)
{
    for(int k=0;k<VTSDIM;++k)
    {
        free(_sdata->GCLV[k]);
        free(_sdata->GCLC[k]);
    }

    free(_sdata->GCLC);
    free(_sdata->GCLV);

    for(int k=0;k<_sdata->VtsBlockNum;++k)
    {
        VtsInfoClean(_sdata->VSF+k);
    }
    free(_sdata->VSF);
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
    for(int k=0;k<_out->nCL[0];++k)
    {
        _out->mask[k] = (int*) malloc(sizeof(int)*_out->nCL[1]);
        _out->offset[k] = (int **) malloc(sizeof(int*)*_out->nCL[1]);
        _out->weight[k] = (VTSDATAFLOAT **) malloc(sizeof(VTSDATAFLOAT*)*_out->nCL[1]);
        for(int j=0;j<_out->nCL[1];++j)
        {
            _out->offset[k][j] = (int *) malloc(sizeof(int)*(VTSDIM+1));
            _out->weight[k][j] = (VTSDATAFLOAT *) malloc(sizeof(VTSDATAFLOAT)*(VTSDIM));
        }
    }

    _out->scores = (unsigned long *) malloc(sizeof(unsigned long)*_sdata->VtsBlockNum);
    for(int k=0;k<_sdata->VtsBlockNum;++k) _out->scores[k] = 0;

    for(int ix=0;ix<_out->nCL[0];ix++)
    {
        for(int jy=0;jy<_out->nCL[1];jy++)
        {
            VTSDATAFLOAT ptmp[3];
            v_copy(ptmp,_out->X0);
            v_add_liner(ptmp,_out->CL[0][ix]-_out->CL[0][0],_out->nX[0],_out->CL[1][jy]-_out->CL[1][0],_out->nX[1]);

            _out->mask[ix][jy] = BlockSearch(_sdata,ptmp);
            if(_out->mask[ix][jy]>=0)
            {
                _out->scores[_out->mask[ix][jy]]++;
                _search(_sdata->VSF+_out->mask[ix][jy],ptmp,_out->offset[ix][jy],_out->weight[ix][jy]);
            }
        }
    }

    fprintf(stdout,"\n");
    for(int k=0;k<_sdata->VtsBlockNum;++k)
    {
        fprintf(stdout,"%5ld,",_out->scores[k]);
        if(k%12==11) fprintf(stdout,"\n");
    }
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

void MergeBlockC(SALEcData * _sdata,int * indices)
{
    // ... alternative method for nearest interpolation in some condition
    // in GetPlaneDataC
}



