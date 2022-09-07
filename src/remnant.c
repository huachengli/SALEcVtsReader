//
// Created by huacheng on 12/24/21.
//

#include <stdio.h>
#include "Utility.h"

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
    GetValueS(ifp,"Plane.output",OutPrexfix,"PLANE");
    GetValueS(ifp,"Plane.Vacuum",VacuumName,"VOF-0");

    int MinStep = GetValueIk(ifp,"SALEc.step",0,"1");
    int MaxStep = GetValueIk(ifp,"SALEc.step",1,"1");
    int PlaneNum = GetValueI(ifp,"Plane.number","1");

    double ** vn,*vd;
    char ** PlaneName;
    vn = (double **) malloc(sizeof(double *)*PlaneNum);
    vd = (double *) malloc(sizeof(double)*PlaneNum);
    PlaneName = (char **) malloc(sizeof(char*)*PlaneNum);
    for(int k=0;k<PlaneNum;++k)
    {
        vn[k] = (double *) malloc(sizeof(double)*VTSDIM);
        PlaneName[k] = (char*) malloc(sizeof(char)*200);
#define GETVNK(dim,y) vn[k][dim] = GetValueDk(ifp,"Plane.n"#y,k,"0.0");
        GETVNK(0,x)
        GETVNK(1,y)
        GETVNK(2,z)
#undef GETVNK

        vd[k] = GetValueDk(ifp,"Plane.d",k,"0.0");
        GetValueSk(ifp,"Plane.name",PlaneName[k],k,"Density");
    }

    int * StepList,StepNum;
    char StepOpt[200];
    GetValueSk(ifp,"Plane.step",StepOpt,0,"range");
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

    CloseInputFile(ifp);

    Plane * Out = (Plane *) malloc(PlaneNum* sizeof(Plane));
    ProfileCache * PCache = (ProfileCache *) malloc((PlaneNum+1)* sizeof(ProfileCache));


    SALEcData * SaleData = (SALEcData *) malloc(sizeof(SALEcData));
    LoadInpInfo(SaleData,SALEcInpName);

    for(int istep=0;istep<StepNum;istep++)
    {

        int CurrentStep = StepList[istep];
        if((CurrentStep> MaxStep) || (CurrentStep<MinStep)) continue;
        sprintf(DataPath,"%s.proc%%d.%d.vts",DataPrefix,CurrentStep);

        LoadVtsData(SaleData,DataPath);
        for(int icache=0;icache<PlaneNum+1;++icache) PCache[icache].CacheState = 0; // set the initial cache state

        fprintf(stdout,"\n");
        for(int k=0;k<PlaneNum;k++)
        {
            Out[k].d = vd[k];
            v_copy(Out[k].n,vn[k]);
            if(head__strcasestr(PlaneName[k],"Remnant_") != NULL)
            {
                strcpy(Out[k].Name,PlaneName[k]+8);
                strcpy(Out[k].Vacuum,VacuumName);

                SetPlaneMeshV(SaleData,Out+k);
                SetPlaneMaskV(SaleData,Out+k);
                // Get the remnant_ distribution
                SumRemant(SaleData,Out+k);
            } else
            {
                // for other parameters in planes.plot do nothing
                // without 'continue' statement, Segmentation fault occurs.
                continue;
            }

            char PlaneOutName[400];
            sprintf(PlaneOutName,"%s.data.step.%d.%d",OutPrexfix,CurrentStep,k);
            fprintf(stdout,"==>Write %s\n",PlaneOutName);
            WritePlaneDataAll(Out+k,PlaneOutName);
            sprintf(PlaneOutName,"%s.coord.step.%d.%d",OutPrexfix,CurrentStep,k);
            WritePlaneCoord(Out+k,PlaneOutName);
            fprintf(stdout,"==>Write %s\n",PlaneOutName);
            CleanPlane(Out+k);
        }
        CleanSALEcData(SaleData);
        CleanCache(PCache,PlaneNum+1);
    }

    free(SaleData);
    free(PCache);
    free(Out);
    free(vd);
    for(int k=0;k<PlaneNum;k++)
    {
        free(vn[k]);
        free(PlaneName[k]);
    }
    free(vn);
    free(PlaneName);
    return 0;
}


