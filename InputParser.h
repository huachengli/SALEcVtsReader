//
// Created by huacheng on 10/19/21.
//
/*
 * A minium version of InputParser From SALEc
 */

#ifndef SALECVTSREADER_INPUTPARSER_H
#define SALECVTSREADER_INPUTPARSER_H
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>

#define MaxStrLen   100
#define MaxKeyId    200

typedef struct {
    char Key[MaxKeyId][MaxStrLen];
    char Value[MaxKeyId][MaxStrLen];
    int  Len;
} InputFile;

InputFile * OpenInputFile(const char fname[]);
InputFile * ParseInputFile(FILE * fp);
void * CloseInputFile(InputFile * ifp);
int SearchInput(InputFile * ifp,const char key[]);

int GetValueI(InputFile * ifp,const char key[], char dvalue[]);
double GetValueD(InputFile * ifp,const char key[], char dvalue[]);
int GetValueS(InputFile * ifp,const char key[], char value[], char dvalue[]);
int GetStrk(char strlist[],int pos, char value[]);
double GetValueDk(InputFile * ifp,const char key[], int k,char dvalue[]);
int GetValueSk(InputFile * ifp,const char key[], char value[],int k,char dvalue[]);
int GetValueIk(InputFile * ifp,const char key[], int k,char dvalue[]);

/*
 * basic string processing
 */
char *trim(char *str);
int Strok(const char _str[],const char _delim[], char value[]);
int InSubset(char _c,const char _set[]);
int ReadLineTrim(unsigned char _buffer[],FILE *fp);

#endif //SALECVTSREADER_INPUTPARSER_H
