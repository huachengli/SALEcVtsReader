//
// Created by huach on 1/16/2024.
//

#ifndef SALECVTSREADER_UTILITY2D_H
#define SALECVTSREADER_UTILITY2D_H

#include "Utility.h"
#include "InputParser.h"


typedef struct
{
    double Asimon;
    double Csimon;
    double Tmelt0;
    double Tfrac;
    double Tdelta;
} MatRef;

typedef struct
{
    char prefix[256];
    double r_planet;
    double center[3];
    int nmat;
    double surf_tem;
    MatRef * ref;
} SALEcPlanetInfo;

typedef struct
{
    int nno;
    char cname[256];
    VTSDATAFLOAT * XXs;
    VTSDATAFLOAT * XXc; // cart
    VTSDATAFLOAT * T; // the temperature need attached

} Citcom2dXX;


void Load2dInpInfo(SALEcData * _sdata,SALEcPlanetInfo * _pdata,const char* _inputName);
void Reset2dZcoord(SALEcData * _sdata);
void LoadCitcom2dXX(Citcom2dXX * _cdata, const char * _cname);
void CleanCitcom2dXX(Citcom2dXX * _cdata);
void ConsistentCitcom2dXX(Citcom2dXX * _cdata, SALEcPlanetInfo * _sinfo);
void ExportCitcom2dT(Citcom2dXX * _cdata);
void ExportCitcom2dTdiff(Citcom2dXX * _cdata, Citcom2dXX * _x);
void ExportCitcom2dTNondim(Citcom2dXX * _cdata, Citcom2dXX * _x);
void InterpolateCitcom2dT(Citcom2dXX * _cdata, SALEcData * _sdata);
void LoadVts2dData(SALEcData * _sdata,const char * _vtsPrefix);
double SimonLiquid(double pre,const MatRef * _s);
int GetVsfCellId(VtsInfo * _vsf, const char * TemVsfName);
void CapLiquidTem(SALEcData * _sdata, SALEcData * _sdata_fill,SALEcPlanetInfo * _pdata);
VTSDATAFLOAT* Vtm2dGetCellData(SALEcData * _sdata, unsigned long k, unsigned long _i,unsigned long _j);
#endif //SALECVTSREADER_UTILITY2D_H
