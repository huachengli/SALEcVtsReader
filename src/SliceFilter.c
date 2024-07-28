//
// Created by lhc on 27/7/2024.
//
/*
 * < slice filter > for 3d SALEc data
 */

#include <stdio.h>
#include "Utility.h"
#include "SliceFilter.h"
#include "VtkWriter.h"

int main(int argc,char * argv[])
{
    char get_plane[200];
    switch (argc) {
        case 1:
            strcpy(get_plane,"planes.plot");
            break;
        case 2:
            strcpy(get_plane,argv[1]);
            break;
        default:
            fprintf(stdout,"error in command line\n");
    }

    InputFile * ifp = OpenInputFile(get_plane);
    char SALEcInpName[200],DataPrefix[200],DataPath[400],OutPrexfix[200],VacuumName[200];

    GetValueS(ifp,"Plane.data",DataPrefix,"Al1100Test");
    GetValueS(ifp,"SALEc.input",SALEcInpName,"SALEc.inp");
    GetValueS(ifp,"Plane.output",OutPrexfix,"SLICE");
    GetValueS(ifp,"Plane.Vacuum",VacuumName,"VOF-0");

    int MinStep = GetValueIk(ifp,"SALEc.step",0,"1");
    int MaxStep = GetValueIk(ifp,"SALEc.step",1,"1");
    int SliceNum = GetValueI(ifp,"Slice.number","1");

    // load slice operation info
    double ** vn,*vd;
    char ** SliceName;
    vn = (double **) malloc(sizeof(double *)*SliceNum);
    vd = (double *) malloc(sizeof(double)*SliceNum);
    SliceName = (char **) malloc(sizeof(char*)*SliceNum);
    for(int k=0;k<SliceNum;++k)
    {
        vn[k] = (double *) malloc(sizeof(double)*VTSDIM);
        vn[k][0] = GetValueDk(ifp, "Slice.nx", k, "0.0");
        vn[k][1] = GetValueDk(ifp, "Slice.ny", k, "0.0");
        vn[k][2] = GetValueDk(ifp, "Slice.nz", k, "0.0");
        vd[k] = GetValueDk(ifp,"Slice.d",k,"0.0");

        SliceNum[k] = (char*) malloc(sizeof(char)*MaxStrLen);
        GetValueSk(ifp,"Slice.name",SliceName[k],k,"XOY");
    }

    int * StepList,StepNum;
    char StepOpt[200];
    GetValueSk(ifp,"Slice.step",StepOpt,0,"range");
    if(strcasecmp(StepOpt,"range")==0)
    {
        int StepList0 = GetValueIk(ifp,"Plane.step",1,"0");
        int StepList1 = GetValueIk(ifp,"Plane.step",2,"0");
        int StepList2 = GetValueIk(ifp,"Plane.step",3,"1");
        StepNum = (StepList1 - StepList0-1)/StepList2 + 1;
        StepList = (int *) malloc(sizeof(int)*StepNum);
        for(int istep=0;istep<StepNum;istep++)
        {
            StepList[istep] = StepList0 + StepList2 * istep;
        }
    }
    else
    {
        fprintf(stdout,"%s:unimplemented step type: (%s)\n",__FILE__,StepOpt);
        exit(1);
    }
    CloseInputFile(ifp);

    // allocate memory for slice
    SliceFilter * OutSlice = malloc(SliceNum * sizeof(SliceFilter));

    SALEcData * SaleData = (SALEcData *) malloc(sizeof(SALEcData));
    LoadInpInfo(SaleData,SALEcInpName);

    for(int istep=0;istep<StepNum;istep++)
    {
        int CurrentStep = StepList[istep];
        if((CurrentStep> MaxStep) || (CurrentStep<MinStep)) continue;
        sprintf(DataPath,"%s.proc%%d.%d.vts",DataPrefix,CurrentStep);
        LoadVtsData(SaleData,DataPath);
        fprintf(stdout,"\n");
        for(int k=0;k<SliceNum;k++)
        {
            OutSlice[k].d = vd[k];
            v_copy(OutSlice[k].n,vn[k]);
            SetPlaneMeshV(SaleData,(Plane *)(OutSlice + k));
            SetSliceMask(SaleData,OutSlice+k,OffsetSerchV);
            GetSliceData(SaleData,OutSlice+k);

            char SliceOutName[MaxStrLen*4];
            sprintf(SliceOutName, "%s.%02d%s.%d.vts", OutPrexfix, k,SliceName[k],CurrentStep);



            WritePlaneDataAll(Out+k, SliceOutName);
            sprintf(SliceOutName, "%s.coord.step.%d.%d", OutPrexfix, CurrentStep, k);
            WritePlaneCoord(Out+k, SliceOutName);
            fprintf(stdout, "==>Write %s\n", SliceOutName);
            CleanSlice(OutSlice+k);
        }
        CleanSALEcData(SaleData);
    }

    free(SaleData);
    free(OutSlice);
    free(vd);
    for(int k=0;k<SliceNum;k++)
    {
        free(vn[k]);
    }
    free(vn);
    return 0;
}

int GetSliceDataC(SALEcData * _sdata, SliceFilter * _out, unsigned long Id)
{
    VtsInfo * _vsf = _sdata->VSF;
    _out->Id = Id;
    if(Id > _vsf->CellNoF) return 0;
    fprintf(stdout,"(%s)",_vsf->CellField[Id].Name);
    _out->NoC = _vsf->CellField[_out->Id].NoC;
    if(NULL != _out->data_vars)
    {
        free(_out->data_vars);
    }
    _out->data_vars = (VTSDATAFLOAT*) malloc(sizeof(VTSDATAFLOAT)*(_out->nCL[0]*_out->nCL[1])*(_out->NoC));
    for(unsigned long ix=0;ix<_out->nCL[0];ix++)
    {
        for(unsigned long jy=0;jy<_out->nCL[1];jy++)
        {
            _out->data[ix][jy] =  _out->data_vars + (ix*_out->nCL[1] + jy)*(_out->NoC);
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

#define SAFEFREE(x) if(NULL!=(x)) free(x);
void CleanSlice(SliceFilter * _out)
{
    for(int k=0;k<_out->nCL[0];++k)
    {

        for(int j=0;j<_out->nCL[1];++j)
        {
            SAFEFREE(_out->offset[k][j])
            SAFEFREE(_out->weight[k][j])
        }
        SAFEFREE(_out->mask[k])
        SAFEFREE(_out->offset[k])
        SAFEFREE(_out->weight[k])
    }
    SAFEFREE(_out->scores)
    SAFEFREE(_out->mask)
    SAFEFREE(_out->offset)
    SAFEFREE(_out->weight)
    SAFEFREE(_out->data_coord)
    SAFEFREE(_out->data_vars)
    SAFEFREE(_out->scores)
}
#undef SAFEFREE

void SetSliceMask(SALEcData * _sdata, SliceFilter * _out,int (*_search)(VtsInfo *, VTSDATAFLOAT *, int *,VTSDATAFLOAT*))
{
    // similar to the old SetPlaneMask,
    // just change the data allocate layout.
    _out->mask = (int **) malloc(sizeof(int*)*_out->nCL[0]);
    _out->offset = (int ***) malloc(sizeof(int**)*_out->nCL[0]);
    _out->weight = (VTSDATAFLOAT ***) malloc(sizeof(VTSDATAFLOAT**)*_out->nCL[0]);
    _out->coord = (VTSDATAFLOAT ***) malloc(sizeof(VTSDATAFLOAT**)*_out->nCL[0]);
    _out->data_coord = malloc(sizeof(VTSDATAFLOAT)*_out->nCL[0]*_out->nCL[1]*(VTSDIM));
    _out->data_vars = NULL;
    _out->data = (VTSDATAFLOAT ***) malloc(sizeof(VTSDATAFLOAT**)*_out->nCL[0]);
    for(int k=0;k<_out->nCL[0];++k)
    {
        _out->mask[k] = (int*) malloc(sizeof(int)*_out->nCL[1]);
        _out->offset[k] = (int **) malloc(sizeof(int*)*_out->nCL[1]);
        _out->weight[k] = (VTSDATAFLOAT **) malloc(sizeof(VTSDATAFLOAT*)*_out->nCL[1]);
        _out->coord[k] = (VTSDATAFLOAT **) malloc(sizeof(VTSDATAFLOAT*)*_out->nCL[1]);
        _out->data[k] = (VTSDATAFLOAT **) malloc(sizeof(VTSDATAFLOAT*)*_out->nCL[1]);
        for(int j=0;j<_out->nCL[1];++j)
        {
            _out->offset[k][j] = (int *) malloc(sizeof(int)*(VTSDIM+1));
            _out->weight[k][j] = (VTSDATAFLOAT *) malloc(sizeof(VTSDATAFLOAT)*(VTSDIM));
            _out->coord[k][j] = _out->data_coord + (k*_out->nCL[1] + j)*VTSDIM;
        }
    }

    _out->scores = (unsigned long *) malloc(sizeof(unsigned long)*_sdata->VtsBlockNum);
    // record number of points need interpolated in every block
    for(int k=0;k<_sdata->VtsBlockNum;++k) _out->scores[k] = 0;

    for(int ix=0;ix<_out->nCL[0];ix++)
    {
        for(int jy=0;jy<_out->nCL[1];jy++)
        {
            // set the coordinate of points on this plane
            VTSDATAFLOAT * ptmp = _out->coord[ix][jy];
            v_copy(ptmp,_out->X0);
            v_add_liner(ptmp,_out->CL[0][ix]-_out->CL[0][0],_out->nX[0],_out->CL[1][jy]-_out->CL[1][0],_out->nX[1]);

            _out->mask[ix][jy] = BlockSearch(_sdata,ptmp); // return the block id of current coordinate
            if(_out->mask[ix][jy]>=0)
            {
                _out->scores[_out->mask[ix][jy]]++;
                _search(_sdata->VSF+_out->mask[ix][jy],ptmp,_out->offset[ix][jy],_out->weight[ix][jy]);
                // get the position of points plane in specific block (offset + weight)
            }
        }
    }
}

int WriteSliceDataAll(SALEcData * _sdata, SliceFilter * _out, const char * _out_name)
{
    FILE * fp = fopen(_out_name,"w");
    if(fp == NULL)
    {
        fprintf(stdout,%s:cannot open %s\n,__func__,_out_name);
        return 0;
    }
    fprintf(stdout, "==>Write %s\n", _out_name);

    char whole_extent[4096], piece_extent[4096];
    snprintf(whole_extent,4096,"%d %lu %d %lu 0 0",1,_out->shape[0]+1,1,_out->shape[1]+1);
    snprintf(piece_extent,4096,"%d %lu %d %lu 0 0",1,_out->shape[0]+1,1,_out->shape[1]+1);
    int len_pointdata = _out->shape[0]*_out->shape[1];
    int len_celldata = (_out->shape[0] - 1)*(_out->shape[1] - 1);
    vts_file_header(fp,whole_extent,piece_extent);
    vtk_point_data_header_with_attr(fp," ");
    // export point data
    for(int k=0;k<_sdata->VSF->CellNoF;++k)
    {
        GetSliceDataC(_sdata,_out,k);
        vtk_dataarray_vecf(fp,_sdata->VSF->CellField[k].Name,"binary",_out->data_vars,len_pointdata,_out->NoC);
    }
    vtk_point_data_trailer(fp);
    vtk_cell_data_header(fp);
    // export cell data
    vtk_cell_data_trailer(fp);
    vtk_output_coordf(fp,"binary",_out->data_coord,len_pointdata);
    vts_file_trailer(fp);

    fclose(fp);
    return 1;
}



