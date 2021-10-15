#include <stdio.h>
#include "VtsReader.h"
#include <math.h>

int TestTypeCast()
{
    int TestIntLen = 50;
    int TestCharLen = sizeof(int) / sizeof(unsigned char) * TestIntLen + 1;

    fprintf(stdout,"%d\n", sizeof(unsigned char));

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
        fprintf(stdout,"orignal size: %d,compressd size:%d\n",TestStrLen,CompLen);
    }

    fprintf(stdout,"rLen:%ld\n",TestStrLen);
    CompErr = uncompress(UncompStr,&UncompLen,CompStr,CompLen);
    fprintf(stdout,"rLen:%ld\n",TestStrLen);

    if(Z_OK!=CompErr)
    {
        fprintf(stdout,"Uncompress error!\n");
    } else
    {
        fprintf(stdout,"orignal size: %d,uncompressd size:%d\n",TestStrLen,UncompLen);
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
