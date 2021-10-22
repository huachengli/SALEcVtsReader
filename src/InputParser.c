//
// Created by huacheng on 10/19/21.
//

#include "InputParser.h"


char * rtrim(char *str)
{
    if (str == NULL || *str == '\0')
    {
        return str;
    }

    int len = strlen(str);
    char *p = str + len - 1;
    while (p >= str  && isspace(*p))
    {
        *p = '\0';
        --p;
    }
    return str;
}

char * ltrim(char *str)
{
    if (str == NULL || *str == '\0')
    {
        return str;
    }

    int len = 0;
    char *p = str;
    while (*p != '\0' && (isspace(*p) || (*p=='.')))
    {
        ++p;
        ++len;
    }
    memmove(str, p, strlen(str) - len + 1);
    return str;
}

char *trim(char *str)
{
    str = rtrim(str);
    str = ltrim(str);
    return str;
}

int InSubset(char _c,const char _set[])
{
    const char * p = _set;
    for(;*p!='\0';p++)
    {
        if(*p==_c)
            return 1;
    }
    return 0;
}

int Strok(const char _str[],const char _delim[], char value[])
{
    const char * p = _str;
    int k=0,i=0;
    for(;*p!='\0';p++)
    {
        if(InSubset(*p,_delim))
        {
            if(k==0)
                continue;
            else
                break;
        } else
        {
            value[k++] = *p;
        }
        i++;
    }
    for(;*p!='\0';p++)
    {
        if(!InSubset(*p,_delim)) break;
        i++;
    }
    value[k] = '\0';
    return i;
}

int ReadLineTrim(unsigned char _buffer[],FILE *fp)
{
    (void )fscanf(fp,"%[^\n]",_buffer);
    fgetc(fp);
    trim(_buffer);
    return strlen(_buffer);
}

int IsComment(const char c[])
{
    if(strlen(c)==0)
        return 1;
    char CommentHead[5]   = "\"#/-*";
    int IsComm = 0;
    for(int k=0;k<5;k++)
    {
        if(c[0] == CommentHead[k])
        {
            IsComm = 1;
            break;
        }
    }
    return IsComm;
}

int InField(const char c[])
{
    const char * p = c;
    for(;*p!='\0';p++)
    {
        if(*p=='.') return 0;
    }
    return 1;
}

void Merge(InputFile * ifp, int start, int mid, int end)
{
    char (*tKey)[MaxStrLen]   = (char (*)[MaxStrLen]) malloc((end+1-start)*sizeof(char [MaxStrLen]));
    char (*tValue)[MaxStrLen] = (char (*)[MaxStrLen]) malloc((end+1-start)*sizeof(char [MaxStrLen]));

    int i = start;
    int j = mid + 1;
    int k = 0;

    while (i <= mid && j <= end)
    {
        if(strcasecmp(ifp->Key[i],ifp->Key[j]) <= 0)
        {
            strcpy(tKey[k],ifp->Key[i]);
            strcpy(tValue[k],ifp->Value[i]);
            k++;i++;
        } else
        {
            strcpy(tKey[k],ifp->Key[j]);
            strcpy(tValue[k],ifp->Value[j]);
            k++;j++;
        }
    }

    while (i <= mid)
    {
        strcpy(tKey[k],ifp->Key[i]);
        strcpy(tValue[k],ifp->Value[i]);
        k++;i++;
    }
    while (j <= end)
    {
        strcpy(tKey[k],ifp->Key[j]);
        strcpy(tValue[k],ifp->Value[j]);
        k++;j++;
    }

    for (i = 0; i < k; i++)
    {
        strcpy(ifp->Key[i+start],tKey[i]);
        strcpy(ifp->Value[i+start],tValue[i]);
    }

    free(tKey);
    free(tValue);
}

void SortInputFile(InputFile * ifp, int start, int end)
{
    if(start >= end)
        return;
    else
    {
        int mid = start + (end - start) / 2;
        SortInputFile(ifp,start,mid);
        SortInputFile(ifp,mid+1,end);
        Merge(ifp,start,mid,end);
    }
}

void Swap(char a[],char b[])
{
    char t[MaxStrLen];
    strcpy(t,a);
    strcpy(a,b);
    strcpy(b,t);
}

InputFile * OpenInputFile(const char fname[])
{
    FILE * fp = fopen(fname,"r");
    InputFile * ifp = (InputFile *) malloc(sizeof(InputFile));
    if(NULL == ifp)
    {
        fprintf(stdout,"Cannot open input file * %s * \n",fname);
        exit(0);
    }
    ifp->Len = 0;

    char LineBuffer[300];
    char MainDelimiter[] = "= ";
    char SubDelimiter[]  = "#*\"";
    char Field[MaxStrLen] = "mesh";

    while(fscanf(fp,"%[^\n]",LineBuffer)!=EOF)
    {
        fgetc(fp);
        trim(LineBuffer);
        int LineBufferLen = strlen(LineBuffer);
        if((LineBufferLen==0) || IsComment(LineBuffer)) continue;
        if((LineBuffer[0]=='[') && (LineBuffer[LineBufferLen-1]==']'))
        {
            LineBuffer[LineBufferLen-1] = '\0';
            strcpy(Field,LineBuffer+1);
            continue;
        }

        char tKey[MaxStrLen], tValue[MaxStrLen];
        int r = Strok(LineBuffer,MainDelimiter,tKey);
        Strok(LineBuffer+r,SubDelimiter,tValue);
        if((0== strlen(tKey)) || (0== strlen(tValue)))
        {
            LineBuffer[0] = '\0';
            continue;
        }

        if(InField(tKey))
        {
            sprintf(ifp->Key[ifp->Len],"%s.%s",Field,tKey);
            strcpy(ifp->Value[ifp->Len], tValue);
        } else
        {
            strcpy(ifp->Key[ifp->Len], tKey);
            strcpy(ifp->Value[ifp->Len], tValue);
        }

        trim(ifp->Key[ifp->Len]);
        trim(ifp->Value[ifp->Len]);
        ifp->Len++;
        LineBuffer[0] = '\0';
    }
    fclose(fp);
    SortInputFile(ifp,0,ifp->Len-1);
    return ifp;
}

InputFile * ParseInputFile(FILE * fp)
{
    InputFile * ifp = (InputFile *) malloc(sizeof(InputFile));
    if(NULL == ifp)
    {
        fprintf(stdout,"Cannot open input file with fp \n");
        exit(0);
    }
    ifp->Len = 0;

    char LineBuffer[300];
    char MainDelimiter[] = "= ";
    char SubDelimiter[]  = "/#*\"";
    char Field[MaxStrLen] = "mesh";

    while(fscanf(fp,"%[^\n]",LineBuffer)!=EOF)
    {
        fgetc(fp);
        trim(LineBuffer);
        int LineBufferLen = strlen(LineBuffer);
        if((LineBufferLen==0) || IsComment(LineBuffer)) continue;
        if((LineBuffer[0]=='[') && (LineBuffer[LineBufferLen-1]==']'))
        {
            LineBuffer[LineBufferLen-1] = '\0';
            strcpy(Field,LineBuffer+1);
            continue;
        }

        char tKey[MaxStrLen], tValue[MaxStrLen];
        int r = Strok(LineBuffer,MainDelimiter,tKey);
        Strok(LineBuffer+r,SubDelimiter,tValue);
        if((0== strlen(tKey)) || (0== strlen(tValue)))
        {
            LineBuffer[0] = '\0';
            continue;
        }

        if(InField(tKey))
        {
            sprintf(ifp->Key[ifp->Len],"%s.%s",Field,tKey);
            strcpy(ifp->Value[ifp->Len], tValue);
        } else
        {
            strcpy(ifp->Key[ifp->Len], tKey);
            strcpy(ifp->Value[ifp->Len], tValue);
        }

        trim(ifp->Key[ifp->Len]);
        trim(ifp->Value[ifp->Len]);
        ifp->Len++;
        LineBuffer[0] = '\0';
    }
    fclose(fp);
    SortInputFile(ifp,0,ifp->Len-1);
    return ifp;
}

void * CloseInputFile(InputFile * ifp)
{
    free(ifp);
}

int SearchInput(InputFile * ifp,const char key[])
{
    int il=0,ir=ifp->Len-1;
    if(strcasecmp(ifp->Key[il],key)>0 || strcasecmp(ifp->Key[ir],key) <0)
        return -1;
    while(il + 1 < ir)
    {
        int mid = (il+ir)/2;
        int caseresult = strcasecmp(ifp->Key[mid],key);
        if(caseresult > 0)
        {
            ir = mid;
        } else if(caseresult < 0)
        {
            il = mid;
        } else
        {
            return mid;
        }
    }
    if(0==strcasecmp(ifp->Key[il],key))
        return il;
    else if(0==strcasecmp(ifp->Key[ir],key))
        return ir;
    else
        return -1;
}

int GetValueI(InputFile * ifp,const char key[], char dvalue[])
{
    int r = SearchInput(ifp,key);
    if(r < 0)
        return strtol(dvalue,NULL,10);
    else
        return strtol(ifp->Value[r],NULL,10);
}

double GetValueD(InputFile * ifp,const char key[], char dvalue[])
{
    int r = SearchInput(ifp,key);
    if(r < 0)
        return strtod(dvalue,NULL);
    else
        return strtod(ifp->Value[r],NULL);
}

int GetValueS(InputFile * ifp,const char key[], char value[], char dvalue[])
{
    int r = SearchInput(ifp,key);
    if(r < 0)
        return strlen(strcpy(value,dvalue));
    else
        return strlen(strcpy(value,ifp->Value[r]));
}

int GetStrk(char strlist[],int pos, char value[])
{
    int len = strlen(strlist);
    if(('['!=strlist[0]) || (']'!=strlist[len-1]))
    {
        return -1;
    }
    int lpos = 0,rpos=0,k=0;
    for(int i=1;i<len-1;i++)
    {
        if(','==strlist[i])
        {
            lpos++;
            continue;
        }
        if(lpos>pos) break;
        if(lpos==pos)
        {
            value[k++] = strlist[i];
        }
    }
    value[k]='\0';
    return k;
}

double GetValueDk(InputFile * ifp,const char key[], int k,char dvalue[])
{
    int r = SearchInput(ifp,key);
    if(r < 0)
        return strtod(dvalue,NULL);
    else
    {
        char tmp[MaxStrLen];
        if(GetStrk(ifp->Value[r],k,tmp))
        {
            return strtod(tmp,NULL);
        } else
        {
            return strtod(dvalue,NULL);
        }
    }
}

int GetValueSk(InputFile * ifp,const char key[], char value[],int k,char dvalue[])
{
    int r = SearchInput(ifp,key);
    if(r<0)
        return strlen(strcpy(value,dvalue));
    else
    {
        char tmp[MaxStrLen];
        GetStrk(ifp->Value[r],k,tmp);
        trim(tmp);
        if(strlen(tmp))
        {
            return strlen(strcpy(value,tmp));
        } else
        {
            return strlen(strcpy(value,dvalue));
        }
    }
}

int GetValueIk(InputFile * ifp,const char key[], int k,char dvalue[])
{
    int r = SearchInput(ifp,key);
    if(r < 0)
        return strtod(dvalue,NULL);
    else
    {
        char tmp[MaxStrLen];
        if(GetStrk(ifp->Value[r],k,tmp))
        {
            return strtol(tmp,NULL,10);
        } else
        {
            return strtol(dvalue,NULL,10);
        }
    }
}

