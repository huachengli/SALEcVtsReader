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

    SALEcData * SaleData = (SALEcData *) malloc(sizeof(SALEcData));
    SALEcPlanetInfo * SaleInfo = malloc(sizeof(SALEcPlanetInfo));
    Load2dInpInfo(SaleData,SaleInfo,SALEcInpName);
    char DataPath[256];
    sprintf(DataPath,"vts/%s.proc%%04d.%04d.vts",SaleInfo->prefix,Refstep);
    LoadVts2dData(SaleData,DataPath);
    Reset2dZcoord(SaleData);
    Citcom2dXX * C2dXX = malloc(sizeof(C2dXX));
    LoadCitcom2dXX(C2dXX,Citcom2dCname);
    ConsistentCitcom2dXX(C2dXX,SaleInfo);
    InterpolateCitcom2dT(C2dXX,SaleData);
    ExportCitcom2dT(C2dXX);
    free(C2dXX);
    free(SaleData);
    free(SaleInfo);
    return 1;
}
