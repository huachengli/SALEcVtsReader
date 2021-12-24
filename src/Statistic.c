//
// Created by huacheng on 12/24/21.
//

#include "Utility.h"

void GetRemnantLim(SALEcData * _sdata, Plane * _out, VTSDATAFLOAT _tol)
{
    // get the id of vacuum and density
    VtsInfo * _vsf = _sdata->VSF;
    _out->Id = 100;
    _out->VacuumId = 100;
    _out->Id = find_cellfield(_out->Name,_vsf);
    _out->VacuumId = find_cellfield(_out->Vacuum,_vsf);

    if(_out->Id==100 || _out->VacuumId==100)
    {
        fprintf(stdout,"cannot find %s or Density in cellfield\n",_out->Name);
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


    // integrate from bottom to top, calculate the target volume
    for(unsigned long ix=0;ix<_out->nCL[0];ix++)
    {
        double ColumnDx = fabs(_sdata->GCLV[0][ix+1] - _sdata->GCLV[0][ix]);
        for(unsigned long jy=0;jy<_out->nCL[1];jy++)
        {
            double ColumnDy = fabs(_sdata->GCLV[1][jy+1] - _sdata->GCLV[1][jy]);
            _out->data[ix][jy][0] = 0.;
            for(unsigned long kz=0;kz<zColumn;++kz)
            {
                double zVacuumVof = VtmGetCellData(_sdata,_out->VacuumId,ix,jy,kz)[0];
                double ColumnDz   = fabs(_sdata->GCLV[2][kz+1] - _sdata->GCLV[2][kz]);
                _out->data[ix][jy][0] += ColumnDx*ColumnDy*ColumnDz*VtmGetCellData(_sdata,_out->Id,ix,jy,kz)[0];
                if(zVacuumVof>_tol)
                {
                    break;
                }
            }
        }
    }

}