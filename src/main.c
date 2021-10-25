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
    char SALEcInpName[200],DataPrefix[200],DataPath[400],OutPrexfix[200];

    GetValueS(ifp,"Plane.data",DataPrefix,"Al1100Test");
    GetValueS(ifp,"SALEc.input",SALEcInpName,"SALEc.inp");
    GetValueS(ifp,"Plane.output",OutPrexfix,"PLANE");

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
        StepNum = (StepList1 - StepList0)/StepList2;
        StepList = (int *) malloc(sizeof(int)*StepNum);
        for(int istep=0;istep<StepNum;istep++)
        {
            StepList[istep] = StepList0 + (StepList1 - StepList0)/StepList2 * istep;
        }
    }

    CloseInputFile(ifp);


    Plane * Out = (Plane *) malloc(PlaneNum* sizeof(Plane));
    SALEcData * SaleData = (SALEcData *) malloc(sizeof(SALEcData));
    LoadInpInfo(SaleData,SALEcInpName);

    for(int istep=0;istep<StepNum;istep++)
    {
        int CurrentStep = StepList[istep];
        if((CurrentStep> MaxStep) || (CurrentStep<MinStep)) continue;
        sprintf(DataPath,"%s.proc%%d.%d.vts",DataPrefix,CurrentStep);

        LoadVtsData(SaleData,DataPath);
        fprintf(stdout,"\n");
        for(int k=0;k<PlaneNum;k++)
        {
            Out[k].d = vd[k];
            v_copy(Out[k].n,vn[k]);
            if(strcasecmp(PlaneName[k],"Profile")==0)
            {
                strcpy(Out[k].Name,"VOF-0");
                SetPlaneMeshV(SaleData,Out+k);
                SetPlaneMaskV(SaleData,Out+k);
                GetProfile(SaleData,Out+k);
            } else
            {
                strcpy(Out[k].Name,PlaneName[k]);
                SetPlaneMeshV(SaleData,Out+k);
                SetPlaneMaskV(SaleData,Out+k);
                GetPlaneDataC(SaleData,Out+k);
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
    }

    free(SaleData);
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



int TestPlane()
{
    const char TestInpFilePath[] = "/home/huacheng/Documents/Github/data/pdata/SALEc_20.inp";
    const char TestDataPrefix[] = "/home/huacheng/Documents/Github/data/pdata/Al1100Test.proc%d.1.vts";
    Plane Out;
    SALEcData SaleData;
    LoadInpInfo(&SaleData,TestInpFilePath);
    LoadVtsData(&SaleData,TestDataPrefix);

    for(int k=0;k<3;++k)
        print_vec(SaleData.BCLV[k],SaleData.nBCLV[k],6);

    vzero(Out.X0);

    Out.n[0] = 0.;
    Out.n[1] = 1.;
    Out.n[2] = 0.;
    Out.d = 0.;

    SetPlaneMeshV(&SaleData,&Out);
    SetPlaneMaskV(&SaleData,&Out);

    strcpy(Out.Name,"Density");
    GetPlaneDataC(&SaleData,&Out);
    WritePlaneData(&Out,0,"Planetmp.txt");
    WritePlaneCoord(&Out,"Planemask.txt");
    CleanSALEcData(&SaleData);
    CleanPlane(&Out);
    return 0;
}


int TestSALEcData()
{
    const char TestInpFilePath[] = "/home/huacheng/Documents/Github/data/pdata/SALEc_20.inp";
    const char TestDataPrefix[] = "/home/huacheng/Documents/Github/data/pdata/Al1100Test.proc%d.1.vts";
    SALEcData SaleData;
    LoadInpInfo(&SaleData,TestInpFilePath);
    LoadVtsData(&SaleData,TestDataPrefix);
    WriteGCL(&SaleData,2,stdout);
    CleanSALEcData(&SaleData);
    return 1;
}



int TestUtility()
{
    const char TestInpFilePath[] = "/home/huacheng/Documents/Github/data/pdata/SALEc_20.inp";
    InputFile * ifp = OpenInputFile(TestInpFilePath);

    int Npgx[VTSDIM],Noffset,Npx[VTSDIM];
    Npgx[0] = GetValueI(ifp,"processor.npgx","8");
    Npgx[1] = GetValueI(ifp,"processor.npgy","8");
    Npgx[2] = GetValueI(ifp,"processor.npgz","8");
    Noffset = GetValueI(ifp,"processor.noffset","8");
    Npx[0] = GetValueI(ifp,"mesh.npx","32");
    Npx[1] = GetValueI(ifp,"mesh.npy","32");
    Npx[2] = GetValueI(ifp,"mesh.npz","32");
    CloseInputFile(ifp);
    unsigned VtsBlockNum = Npgx[0]*Npgx[1]*Npgx[2];
    VtsInfo * VSF = malloc(sizeof(VtsInfo)*VtsBlockNum);

    const char TestDataPrefix[] = "/home/huacheng/Documents/Github/data/pdata/Al1100Test.proc%d.1.vts";
    unsigned int TaskFinished = 0;
#pragma omp parallel for num_threads(12) default(shared)
    for(int k=0;k<VtsBlockNum;k++)
    {
        char VtsFileName[200];
        sprintf(VtsFileName,TestDataPrefix,k);
        FILE *fp = fopen(VtsFileName,"r");
        if(NULL==fp)
        {
            fprintf(stdout,"cannot open %s\n",VtsFileName);
            exit(0);
        }
        VtsLoad(VSF+k,fp);
        fclose(fp);

#pragma omp critical
        {
            TaskFinished ++;
            if(TaskFinished%10==9)
            {
                fprintf(stdout,"#");
                fflush(stdout);
            }
        };
    }

    VTSDATAFLOAT ** GCL; // global coordinate line
    GCL = (VTSDATAFLOAT **) malloc(sizeof(VTSDATAFLOAT*)*VTSDIM);
    for(int k=0;k<VTSDIM;++k)
    {
        GCL[k] = (VTSDATAFLOAT*) malloc(sizeof(VTSDATAFLOAT)*(Npgx[k]*Npx[k]+1));
    }

    // Set GCL X-direction
    for(int k=0;k<Npgx[0];++k)
    {
        memcpy(GCL[0]+k*Npx[0], VSF[k].CLV[0] + Noffset, sizeof(VTSDATAFLOAT) * (Npx[0] + 1));
    }
    // Set GCL Y-direction
    for(int k=0;k<Npgx[1];++k)
    {
        memcpy(GCL[1]+k*Npx[1], VSF[k*Npgx[0]].CLV[1] + Noffset, sizeof(VTSDATAFLOAT) * (Npx[1] + 1));
    }
    // Set GCL Z-direction
    for(int k=0;k<Npgx[2];++k)
    {
        memcpy(GCL[2]+k*Npx[2], VSF[k*Npgx[0]*Npgx[1]].CLV[2] + Noffset, sizeof(VTSDATAFLOAT) * (Npx[2] + 1));
    }

    FILE *fp = fopen("GCLtmp.csv","w");
    if(NULL==fp)
    {
        fprintf(stdout,"cannot open GCLtmp file!\n");
        exit(0);
    }

    for(int k=0;k<(Npx[2]*Npgx[2]+1);++k)
    {
        fprintf(fp,"%10.6f,",GCL[2][k]);
        fprintf(stdout,"%10.6f,",GCL[2][k]);
        if(k%15==14)
        {
            fprintf(stdout,"\n");
//            fprintf(fp,"\n");
        }

    }
    fclose(fp);

    for(int k=0;k<VTSDIM;++k)
        free(GCL[k]);
    free(GCL);


    for(int k=0;k<VtsBlockNum;++k)
    {
        VtsInfoClean(VSF+k);
    }
    free(VSF);
    return 1;
}


int TestOmpVtsLoad()
{
    const char TestDataPrefix[] = "/home/huacheng/Documents/Github/data/pdata/Al1100Test.proc%d.1.vts";
    VtsInfo VSF[288];
    int TaskFinished = 0;
#pragma omp parallel for num_threads(6) shared(TaskFinished,TestDataPrefix,stdout,VSF) default(none)
    for(int k=0;k<288;k++)
    {
        char VtsFileName[200];
        sprintf(VtsFileName,TestDataPrefix,k);
        FILE *fp = fopen(VtsFileName,"r");
        VtsLoad(VSF+k,fp);
        fclose(fp);

#pragma omp critical
        {
            TaskFinished ++;
            if(TaskFinished%10==9)
            {
                fprintf(stdout,"#");
                fflush(stdout);
            }
        };

    }

    for(int k=0;k<288;k++)
    {
        VtsInfoClean(VSF+k);
    }
    return 1;
}


int TestLoadVts()
{
    const char TestDataFile[] = "/home/huacheng/Documents/Github/data/pdata/Al1100Test.proc189.1.vts";

    FILE * fp = fopen(TestDataFile,"r");
    if(NULL==fp)
    {
        fprintf(stdout,"can not open %s\n",TestDataFile);
        exit(0);
    }

    VtsInfo VSF;

    VtsLoad(&VSF,fp);

    fprintf(stdout,"\n");
    for(unsigned long k=0;k<69;k+=2)
    {
        VTSDATAFLOAT * _tpoint = VtsGetPoint(&VSF,k,0,0);
        fprintf(stdout,"%04ld=[%20.10f,%20.10f,%20.10f]\n",k,_tpoint[0],_tpoint[1],_tpoint[2]);
    }

    fclose(fp);
    VtsInfoClean(&VSF);
    return 1;
}


int TestBinaryRead()
{
    const char TestDataFile[] = "/home/huacheng/Documents/Github/data/pdata/Al1100Test.proc188.1.vts";
    int DataLine = 7;

    FILE * fp = fopen(TestDataFile,"r");
    if(NULL==fp)
    {
        fprintf(stdout,"can not open %s\n",TestDataFile);
        exit(0);
    }

    char LineBuffer[1024];
    for(int k=0;k<DataLine-1;++k)
    {
        fscanf(fp,"%[^\n]",LineBuffer);
        fgetc(fp);
        fprintf(stdout,"%s\n",LineBuffer);
    }

    float ** FdataArray;
    unsigned long * FdataLen;
    ReadVtsBinaryF32(FdataArray,FdataLen,fp);
    fclose(fp);
    return 0;
}


int TestTypeCast()
{
    int TestIntLen = 50;
    int TestCharLen = sizeof(int) / sizeof(unsigned char) * TestIntLen + 1;

    fprintf(stdout,"%ld\n", sizeof(unsigned char));

    int *IntArray = (int *) malloc(sizeof(int) * TestIntLen);
    unsigned char *ChrArray1 = (unsigned char *) malloc(sizeof(unsigned char) * TestCharLen);
    unsigned char *ChrArray2 = (unsigned char *) malloc(sizeof(unsigned char) * TestCharLen);

    int *IntArray1 = (int *) malloc(sizeof(int) * TestIntLen);
    int *IntArray2 = (int *) malloc(sizeof(int) * TestIntLen);

    ChrArray1[TestCharLen] = 0;
    ChrArray2[TestCharLen] = 0;
    for (int k = 0; k < TestIntLen; k++)
    {
        IntArray[k] = rand();
    }
    IntToUnsignedChar(IntArray,TestIntLen,ChrArray1);
    memcpy(ChrArray2,IntArray, sizeof(int)*TestIntLen);

    fprintf(stdout,"Cmp result:%d\n", strcmp(ChrArray2,ChrArray1));

    UnsignedCharToInt(ChrArray1,TestCharLen-1,IntArray1);
    memcpy(IntArray2,ChrArray2,sizeof(int)*TestIntLen);

    for(int k=0;k<TestIntLen;k++)
    {
        fprintf(stdout,"%d|%d|%d, ",IntArray1[k],IntArray2[k],IntArray[k]);
        if(k%5==4) fprintf(stdout,"\n");
    }
    return 0;
}

int TestCompress()
{
    const unsigned char * TestStr = "My guess is that there is some issue with memory. But how can I know more about? Did I not ";
    const unsigned int TestStrLen = strlen(TestStr) + 1;

    unsigned char CompStr[200],UncompStr[200];
    unsigned long CompLen = sizeof(CompStr)/ sizeof(CompStr[0]);
    unsigned long UncompLen = sizeof(UncompStr)/ sizeof(UncompStr[0]);
    int CompErr = Z_OK;

    unsigned char CompStr_I[200];
    unsigned int CompLen_I = sizeof(CompStr_I)/ sizeof(CompStr_I[0]);

//    CompErr = compress(CompStr,&CompLen,(const Bytef*)TestStr,TestStrLen);

    unsigned char * tCompStr;
    zlibcompress((unsigned char*)TestStr,TestStrLen,&tCompStr,&CompLen_I);

//    fprintf(stdout,"%s\n%s\n%d\n",tCompStr,CompStr,strcmp(tCompStr,CompStr));

    strcpy(CompStr,tCompStr);
    CompLen = CompLen_I;

    if(Z_OK!=CompErr)
    {
        fprintf(stdout,"compress error:%s\n",CompStr);
        exit(0);
    } else
    {
        fprintf(stdout,"orignal size: %d,compressd size:%ld\n",TestStrLen,CompLen);
    }

    fprintf(stdout,"rLen:%d\n",TestStrLen);
    CompErr = uncompress(UncompStr,&UncompLen,CompStr,CompLen);
    fprintf(stdout,"rLen:%d\n",TestStrLen);

    if(Z_OK!=CompErr)
    {
        fprintf(stdout,"Uncompress error!\n");
    } else
    {
        fprintf(stdout,"orignal size: %d,uncompressd size:%ld\n",TestStrLen,UncompLen);
    }

    if(strcmp(UncompStr,TestStr))
    {
        fprintf(stdout,">%s\n>%s\n",TestStr,UncompStr);
    } else
    {
        fprintf(stdout,"uncompress succeed: >%s\n",UncompStr);
    }
    return 0;
}

int TestBase64()
{
    unsigned char bStr[300];
    unsigned char bCode[300];
    unsigned char bStr2[300];

    unsigned int ErrorNum = 0;
    for(int i=0;i<2000000;i++)
    {
        unsigned int TestLen = rand()%200;
//        unsigned int TestLen = 5;
        RandomStr(bStr,TestLen);
        bStr[TestLen] = '\0';
        unsigned int sLen = strlen(bStr);

        Base64Encode(bStr,sLen,bCode);
        unsigned int cLen = (sLen/3 + ((sLen%3)>0))*4;
        bCode[cLen] = '\0';

        Base64Decode(bCode,cLen,bStr2);

        if(strcmp(bStr,bStr2))
        {
            fprintf(stdout,"  %d|%d|%s\n=>%s\n=>%s\n", 0!=strcmp(bStr,bStr2),TestLen,bStr,bCode,bStr2);
            ErrorNum++;
        }

        if(i%10000 == 1)
        {
            fprintf(stdout,"#");
        }
    }
    fprintf(stdout,"\n ErrorNum:%d",ErrorNum);
    return 0;
}
