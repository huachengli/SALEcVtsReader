//
// Created by lhc on 25/4/2024.
//
// read vts data and generate profile automatically
#include "Utility.h"
#include <stdio.h>
#include "Utility2d.h"
#include <unistd.h>

void Vts2dProfile(SALEcData * _sdata, SALEc2dProfile * _spdata);
void Write2dProfile(SALEc2dProfile * _spdata, const char * fname, const char * format);
void Clean2dProfile(SALEc2dProfile * _spdata);

int main(int argc,char * argv[])
{

    char SALEcDataDir[256] = ".";
    char SALEcInpName[256] = "sale2d.inp";
    int Expstep = 1;

    int c;
    int f_set_flag = 0;
    char f_set_name[256];
    int f_bin_out = 0;
    char o_set_name[256] = "profile_2d";
    opterr = 0; // defined in getopt*.h
    while ((c = getopt (argc, argv, "d:f:e:o:b")) != -1)
    {
        switch (c)
        {
            case 'd':
                strcpy(SALEcDataDir,optarg);
                break;
            case 'f':
                strcpy(f_set_name,optarg);
                f_set_flag = 1;
                break;
            case 'e':
                Expstep = atoi(optarg);
                break;
            case 'b':
                f_bin_out = 1;
                break;
            case 'o':
                strcpy(o_set_name,optarg);
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
    } else
    {
        sprintf(SALEcInpName,"%s/%s",SALEcDataDir,f_set_name);
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

    // write profile to csv files
    strcat(o_set_name,".csv");
    Write2dProfile(&profile,o_set_name,"csv");

    if(f_bin_out)
    {
        strcat(o_set_name,".bin");
        Write2dProfile(&profile,o_set_name,"bin");
    }

    CleanSALEcData(SaleData);
    Clean2dProfile(&profile);
    free(SaleData);

    // logs
    fprintf(stdout,"write profile[%d] to %s\n",Expstep,o_set_name);
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

    /*
     * just used for debug write the vts data to binary files without compression
     */

//    FILE * bfp = fopen("profile_test.bin","wb");
//    for(int i=0;i<_sdata->nGCLC[0];++i) {
//        for (int j = 0; j < _sdata->nGCLC[1]; ++j) {
//            VTSDATAFLOAT *vof_f = Vtm2dGetCellData(_sdata, VofId, i, j);
//            VTSDATAFLOAT *den_f = Vtm2dGetCellData(_sdata, DenId, i, j);
//            fwrite(vof_f, sizeof(VTSDATAFLOAT),1,bfp);
//        }
//    }
//    fclose(bfp);

    for(int i=0;i<_sdata->nGCLC[0];++i)
    {
        _spdata->pos[i] = _sdata->GCLC[0][i];
        _spdata->pos2[i] = _sdata->GCLC[1][0];
        for(int j = 0; j < _sdata->nGCLC[1]; ++j) {
            VTSDATAFLOAT *vof_f = Vtm2dGetCellData(_sdata, VofId, i, j);
            VTSDATAFLOAT *den_f = Vtm2dGetCellData(_sdata, DenId, i, j);

            if(vof_f[0] > 0.5 || den_f[0] < 10.0)
            {
                break;
            }
            _spdata->pos2[i] = _sdata->GCLC[1][j];
        }
    }
}

void Write2dProfile(SALEc2dProfile * _spdata, const char * fname, const char * format)
{
    if(strcasecmp(format,"csv") == 0)
    {
        FILE * fp = fopen(fname,"w");
        for(int k=0;k<_spdata->n;++k)
        {
            fprintf(fp,"%10.5e, %10.5e\n",_spdata->pos[k],_spdata->pos2[k]);
        }
        fclose(fp);
    }

    if(strcasecmp(format,"bin") == 0)
    {
        FILE * fp = fopen(fname,"wb");
        fwrite(_spdata->pos, sizeof(VTSDATAFLOAT),_spdata->n, fp);
        fwrite(_spdata->pos2, sizeof(VTSDATAFLOAT),_spdata->n, fp);
        fclose(fp);
    }
}


void Clean2dProfile(SALEc2dProfile * _spdata)
{
    // clean _spdata
    free(_spdata->pos2);
    free(_spdata->pos);
}