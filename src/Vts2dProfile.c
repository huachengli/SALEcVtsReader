//
// Created by lhc on 25/4/2024.
//
// read vts data and generate profile automatically
#include "Utility.h"
#include <stdio.h>
#include "Utility2d.h"
#include <unistd.h>

void Vts2dProfile(SALEcData * _sdata, SALEc2dProfile * _spdata);
void Write2dProfile(SALEc2dProfile * _spdata, const char * fname);

int main(int argc,char * argv[])
{

    char SALEcDataDir[256] = ".";
    char SALEcInpName[256] = "sale2d.inp";
    int Expstep = 1;

    int c;
    int f_set_flag = 0;
    opterr = 0; // defined in getopt*.h
    while ((c = getopt (argc, argv, "d:f:r:e:c:")) != -1)
    {
        switch (c)
        {
            case 'd':
                strcpy(SALEcDataDir,optarg);
                break;
            case 'f':
                strcpy(SALEcInpName,optarg);
                f_set_flag = 1;
                break;
            case 'e':
                Expstep = atoi(optarg);
                break;
            case '?':
                fprintf (stderr,"Unknown option character `\\x%x'.\n",optopt);
                return 1;
            default:
                abort();
        }
    }

    if(f_set_flag == 0)
    {
        sprintf(SALEcInpName,"%s/sale2d.inp",SALEcDataDir);
    }

    // load vts data
    // load the reference data
    SALEcData * SaleData = (SALEcData *) malloc(sizeof(SALEcData));
    SALEcPlanetInfo * SaleInfo = malloc(sizeof(SALEcPlanetInfo));
    Load2dInpInfo(SaleData, SaleInfo, SALEcInpName);
    char DataPath[256];
    snprintf(DataPath,256,"%s/vts/%s.proc%%04d.%04d.vts",SALEcDataDir,SaleInfo->prefix,Expstep);
    LoadVts2dData(SaleData, DataPath);

    SALEc2dProfile profile = {.n=0,.pos=NULL,.pos2=NULL};
    Vts2dProfile(SaleData,&profile);
    Write2dProfile(&profile,"profile_test.csv");

    CleanSALEcData(SaleData);
    free(SaleData);
    return 0;
}

void Vts2dProfile(SALEcData * _sdata, SALEc2dProfile * _spdata)
{
    // lets allocate mem
    _spdata->n = _sdata->nGCLC[0];
    _spdata->pos = (VTSDATAFLOAT *) malloc(sizeof(VTSDATAFLOAT)*_spdata->n);
    _spdata->pos2= (VTSDATAFLOAT *) malloc(sizeof(VTSDATAFLOAT)*_spdata->n);

    // get density field id
    unsigned long VofId = find_cellfield("e_vof",_sdata->VSF);
    unsigned long DenId = find_cellfield("e_den",_sdata->VSF);


    for(int i=4;i<_sdata->nGCLC[0];++i)
    {
        _spdata->pos[i] = _sdata->GCLC[0][i];
        _spdata->pos2[i] = _sdata->GCLC[1][0];
        for(int j = 0; j < _sdata->nGCLC[1]; ++j) {
            VTSDATAFLOAT *vof_f = Vtm2dGetCellData(_sdata, VofId, i, j);
            VTSDATAFLOAT *den_f = Vtm2dGetCellData(_sdata, DenId, i, j);

            fprintf(stdout,"(%d,%d):%e,%f\n",i,j,_sdata->GCLC[1][j],*vof_f);
            if(vof_f[0] > 0.5)
            {
                break;
            }
//            if(den_f[0] < 50.0) break;

            _spdata->pos2[i] = _sdata->GCLC[1][j];

        }
        break;
    }
}

void Write2dProfile(SALEc2dProfile * _spdata, const char * fname)
{
    FILE * fp = fopen(fname,"w");

    for(int k=0;k<_spdata->n;++k)
    {
        fprintf(fp,"%10.5e, %10.5e\n",_spdata->pos[k],_spdata->pos2[k]);
    }
    fclose(fp);

    // clean _spdata
    free(_spdata->pos2);
    free(_spdata->pos);
}