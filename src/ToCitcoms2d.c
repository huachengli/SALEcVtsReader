//
// Created by huach on 1/16/2024.
//
//  read vts data from SALEc2d and export binary data for citcoms2d/saxi
#include "Utility.h"
#include <stdio.h>
#include "Utility2d.h"
#include <unistd.h>

int main(int argc,char * argv[])
{
    char SALEcDataDir[256] = ".";
    char SALEcInpName[256] = "sale2d.inp";
    char Citcom2dCname[256] = "chicxulub.coord";
    int Refstep = 0;
    int Expstep = 1;

    int c;
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
                break;
            case 'r':
                Refstep = atoi(optarg);
                break;
            case 'e':
                Expstep = atoi(optarg);
                break;
            case 'c':
                strcpy(Citcom2dCname,optarg);
                break;
            case '?':
                if(optopt == 'c')
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                else if (isprint(optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf (stderr,"Unknown option character `\\x%x'.\n",optopt);
                return 1;
            default:
                abort();
        }
    }

    // load the reference data
    SALEcData * RefSaleData = (SALEcData *) malloc(sizeof(SALEcData));
    SALEcPlanetInfo * SaleInfo = malloc(sizeof(SALEcPlanetInfo));
    Load2dInpInfo(RefSaleData, SaleInfo, SALEcInpName);
    char DataPath[256];
    snprintf(DataPath,256,"vts/%s.proc%%04d.%04d.vts",SaleInfo->prefix,Refstep);
    LoadVts2dData(RefSaleData, DataPath);

    // load the step need exported, in most time, the metadata of exp and ref step is
    // exactly same. In practice, the ExpSaleData should be duplication of Ref*
    // in following, just load those data again from input file for simplicity
    SALEcData * ExpSaleData = (SALEcData *) malloc(sizeof(SALEcData));
    snprintf(DataPath,256,"vts/%s.proc%%04d.%04d.vts",SaleInfo->prefix,Expstep);
    Load2dInpInfo(ExpSaleData, NULL, SALEcInpName);
    LoadVts2dData(ExpSaleData, DataPath);
    CapLiquidTem(ExpSaleData,RefSaleData,SaleInfo);

    // load coordinates generated by Citcoms2d for ref step
    Citcom2dXX * RefC2dXX = malloc(sizeof(RefC2dXX));
    LoadCitcom2dXX(RefC2dXX, Citcom2dCname);
    ConsistentCitcom2dXX(RefC2dXX, SaleInfo);
    InterpolateCitcom2dT(RefC2dXX, RefSaleData);

    // for exp step, load in interpolate again
    Citcom2dXX * ExpC2dXX = malloc(sizeof(RefC2dXX));
    LoadCitcom2dXX(ExpC2dXX, Citcom2dCname);
    ConsistentCitcom2dXX(ExpC2dXX, SaleInfo);
    InterpolateCitcom2dT(ExpC2dXX, ExpSaleData);

    ExportCitcom2dTNondim(ExpC2dXX,RefC2dXX);

    CleanSALEcData(ExpSaleData);
    CleanSALEcData(RefSaleData);
    free(ExpSaleData);
    free(RefSaleData);

    CleanCitcom2dXX(RefC2dXX);
    CleanCitcom2dXX(ExpC2dXX);
    free(RefC2dXX);
    free(ExpC2dXX);
    free(SaleInfo);
    return 1;
}
