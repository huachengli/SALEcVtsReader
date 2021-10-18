//
// Created by huacheng on 10/14/21.
//

#include "VtsReader.h"

void Base64Encode(unsigned char * _str, unsigned int _slen, unsigned char * _code)
{
    unsigned char Base64Table[] ="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    int i = 0;

    for(;i+2<_slen;i+=3)
    {
        _code[i/3*4 + 0] = Base64Table[_str[i] >> 2];
        _code[i/3*4 + 1] = Base64Table[((_str[i]&0x03)<<4)|((_str[i+1]&0xf0)>>4)];
        _code[i/3*4 + 2] = Base64Table[((_str[i+1]&0x0f)<<2)|((_str[i+2]&0xc0)>>6)];
        _code[i/3*4 + 3] = Base64Table[_str[i+2]&0x3f];
    }

    switch (_slen-i)
    {
        case 1:
            _code[i/3*4 + 0] = Base64Table[_str[i] >> 2];
            _code[i/3*4 + 1] = Base64Table[(_str[i]&0x03)<<4];
            _code[i/3*4 + 2] = '=';
            _code[i/3*4 + 3] = '=';
            break;
        case 2:
            _code[i/3*4 + 0] = Base64Table[_str[i] >> 2];
            _code[i/3*4 + 1] = Base64Table[((_str[i]&0x03)<<4)|((_str[i+1]&0xf0)>>4)];
            _code[i/3*4 + 2] = Base64Table[((_str[i+1]&0x0f)<<2)];
            _code[i/3*4 + 3] = '=';
            break;
        case 0:
            break;
        default:
            fprintf(stdout,"error in Base64Encode!\n");
            exit(0);
            break;
    }
}

int Base64Decode(unsigned char * _code, unsigned int _clen, unsigned char * _str)
{
    unsigned int BaseStrTable[] ={0,0,0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,62,0,0,0,
                 63,52,53,54,55,56,57,58,
                 59,60,61,0,0,0,0,0,0,0,0,
                 1,2,3,4,5,6,7,8,9,10,11,12,
                 13,14,15,16,17,18,19,20,21,
                 22,23,24,25,0,0,0,0,0,0,26,
                 27,28,29,30,31,32,33,34,35,
                 36,37,38,39,40,41,42,43,44,
                 45,46,47,48,49,50,51};

    unsigned _slen = 0;
    if(strstr(_code,"=="))
    {
        _slen = _clen/4*3 - 2;
    } else if(strstr(_code,"="))
    {
        _slen = _clen/4*3 - 1;
    } else
    {
        _slen = _clen/4*3;
    }

//    _str = (unsigned char *) malloc((_slen+1)*sizeof(unsigned char));
    _str[_slen] = '\0';

    int i = 0;
    for(;i*3/4+2<_slen;i+=4)
    {
        _str[i*3/4 + 0] = (BaseStrTable[_code[i]]<<2)|(BaseStrTable[_code[i+1]]>>4);
        _str[i*3/4 + 1] = (BaseStrTable[_code[i+1]]<<4)|(BaseStrTable[_code[i+2]]>>2);
        _str[i*3/4 + 2] = (BaseStrTable[_code[i+2]]<<6)|(BaseStrTable[_code[i+3]]);
    }

    switch (_slen%3)
    {
        case 1:
            _str[_slen-1] = (BaseStrTable[_code[i]]<<2)|(BaseStrTable[_code[i+1]]>>4);
            break;
        case 2:
            _str[_slen-2] = (BaseStrTable[_code[i]]<<2)|(BaseStrTable[_code[i+1]]>>4);
            _str[_slen-1] = (BaseStrTable[_code[i+1]]<<4)|(BaseStrTable[_code[i+2]]>>2);
            break;
        case 0:
            break;
        default:
            fprintf(stdout,"error in Base64Decode!\n");
            exit(0);
            break;
    }
    return _slen;
}

void RandomStr(unsigned char * _str,int _slen)
{
    unsigned char SourceStr[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz,./;\"'<>?";
    unsigned int SourceStrLen = strlen(SourceStr);
    for(int k=0;k<SourceStrLen;++k)
    {
        _str[k] = SourceStr[rand()%SourceStrLen];
    }
    return;
}


#define CHUNK 16384
void zlibcompress(unsigned char* in, int nn, unsigned char** out, int *nn2)
/* function to compress "in" to "out" reducing size from nn to nn2 */
{
    int ntemp=0;

    /* in and out of z-stream */
    unsigned char inz[CHUNK];
    unsigned char outz[CHUNK];

    /* compression level */
    int level = Z_DEFAULT_COMPRESSION;
    int ret,flush;
    int i,j,k;

    /* zlib compression stream */
    z_stream strm;

    /* hope compressed data will be <= uncompressed */
    *out = malloc(sizeof(unsigned char)*nn);

    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;

    /* zlib init */
    ret = deflateInit(&strm, level);
    if (ret == Z_OK){
        i=0;     // position in "in" array
        do{
            j=0; // position in "inz"
            do{
                inz[j++]=in[i++];
            } while((j<CHUNK) && (i<nn)); // stopps if "inz"-buffer is full or "in" array empty
            strm.avail_in=j;              // set number of input chars

            flush = (i==nn) ? Z_FINISH : Z_NO_FLUSH; // done?
            strm.next_in = inz;           // set input buffer

            do{
                strm.avail_out = CHUNK;   // set number of max output chars
                strm.next_out = outz;     // set output buffer

                /* zlib compress */
                ret = deflate(&strm, flush);
                assert(ret != Z_STREAM_ERROR);

                /* zlib changed strm.avail_out=CHUNK
                 to the number of chars we can NOT use
                 in outz */

                for (k=0;k<CHUNK-strm.avail_out;k++){
                    (*out)[ntemp+k]=outz[k];
                }

                /* increase position in "out" */
                ntemp+=(CHUNK-strm.avail_out);
            }while(strm.avail_out==0);
            assert(strm.avail_in == 0);

        }while (flush != Z_FINISH);
    }
    else{fprintf(stderr,"Error during compression init\n");}

    // now we know how short "out" should be!
    *nn2=ntemp;
    *out = realloc(*out,sizeof(unsigned char)*ntemp);

    (void)deflateEnd(&strm);
    return;
}


void IntToUnsignedChar(int * _iarray, int _inn, unsigned char * _carray)
{
    union IntToUnsignedChars
    {
        int Input;
        unsigned char Output[4];
    } _inttransform;

    for(int k=0;k<_inn;k++)
    {
        _inttransform.Input = _iarray[k];
        _carray[k*4+0] = _inttransform.Output[0];
        _carray[k*4+1] = _inttransform.Output[1];
        _carray[k*4+2] = _inttransform.Output[2];
        _carray[k*4+3] = _inttransform.Output[3];
    }
}

void UnsignedCharToInt(unsigned char * _carray, int _cnn, int * _iarray)
{
    union IntToUnsignedChars
    {
        int Output;
        unsigned char Input[4];
    } _inttransform;

    for(int k=0;k<_cnn;k+=4)
    {
        _inttransform.Input[0]=_carray[k+0];
        _inttransform.Input[1]=_carray[k+1];
        _inttransform.Input[2]=_carray[k+2];
        _inttransform.Input[3]=_carray[k+3];
        _iarray[k/4] = _inttransform.Output;
    }
}

void FloatToUnsignedChar(float * _farray, int _fnn, unsigned char * _carray)
{

    union FloatToUnsignedChars
    {
        float Input;
        unsigned char Output[4];
    } _floattransform;

    for(int k=0;k<_fnn;k++)
    {
        _floattransform.Input = _farray[k];
        _carray[4 * k + 0] = _floattransform.Output[0];
        _carray[4 * k + 1] = _floattransform.Output[1];
        _carray[4 * k + 2] = _floattransform.Output[2];
        _carray[4 * k + 3] = _floattransform.Output[3];
    }
}

void ReadVtsBinaryF32(float ** _data,unsigned long * _dlen,FILE * fp)
{
    unsigned char HeadB64[24],HeadChar[16];
    fread(HeadB64, sizeof(char),24,fp);
    Base64Decode(HeadB64,24,HeadChar);
    int HeadInt[4];
    /*
     * HeadInt[0]= [const int]
     * HeadInt[1] = HeadInt[2] length of data uncompressed
     * HeadInt[3] = length of data compressed
     */
    memcpy(HeadInt,HeadChar,16);

    unsigned char *BodyB64,*BodyComp,*BodyUncomp;
    unsigned long BodyB64Len = 4*ceil((double) HeadInt[3]/3.0);
    unsigned long BodyCompLen = HeadInt[3];
    unsigned long BodyUncompLen = HeadInt[2];

    BodyB64 = (unsigned char*) malloc(sizeof(unsigned char)*BodyB64Len);
    fread(BodyB64, sizeof(unsigned char),BodyB64Len,fp);

    BodyComp = (unsigned char*) malloc(sizeof(unsigned char)*BodyCompLen);
    Base64Decode(BodyB64,BodyB64Len,BodyComp);
    free(BodyB64);

    BodyUncomp = (unsigned char*) malloc(sizeof(unsigned char)*BodyUncompLen);
    uncompress(BodyUncomp,&BodyUncompLen,BodyComp,BodyCompLen);
    free(BodyComp);

    *_dlen = BodyUncompLen/4;
    *_data = (float *) malloc(sizeof(float)*(*_dlen));
    memcpy(*_data,BodyUncomp,BodyUncompLen);
    free(BodyUncomp);
    fgetc(fp);
}


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
    fscanf(fp,"%[^\n]",_buffer);
    fgetc(fp);
    trim(_buffer);
    return strlen(_buffer);
}


void VtsLoad(VtsInfo * _vfp,FILE * fp)
{
    _vfp->VtsStack = (VtsStackFrame *) malloc(sizeof(VtsStackFrame)*MaxStackDepth);
    _vfp->StackPos = 0;
    _vfp->CellField = (VtsData *) malloc(sizeof(VtsData)*MaxStackDepth);
    _vfp->CellNoF = 0;
    _vfp->PointField = (VtsData *) malloc(sizeof(VtsData)*MaxStackDepth);
    _vfp->PointNoF = 0;

    VtsFrameHeadLoad(_vfp,fp);
    while(VtsFrameLoad(_vfp,fp))
    {
        fprintf(stdout,"$\n");
        VtsStackFrame * _vsf = _vfp->StackPos - 1 + _vfp->VtsStack;
        if(_vsf->Tag == SALEC_VTS_DATAARRAY)
        {
            ReadVtsBinaryF32(&(_vfp->ActiveVtsData->Data),&(_vfp->ActiveVtsData->DataLen),fp);
        }
    }

    VtsCoordinateReshape(_vfp);
}

int VtsFrameHeadLoad(VtsInfo * _vfp,FILE *fp)
{
    /*
     * Read the first line of vts file
     */

    unsigned char LineBuffer[1024];
    ReadLineTrim(LineBuffer,fp);
    VtsStackFrame * _vsf = _vfp->VtsStack + _vfp->StackPos;

    unsigned char SALEcVtsHead[] = "<?xml version=\"1.0\"?>";
    if(0!= strcmp(SALEcVtsHead,LineBuffer))
    {
        fprintf(stdout,"Wranning/the header of vts is not consistent with SALEc!\n");
        exit(0);
    }

    strcpy(_vsf->Name,LineBuffer);
    _vsf->Tag = SALEC_VTS_HEADER;
    _vfp->StackPos ++;
    return 1;
}

int VtsFrameLoad(VtsInfo * _vsf,FILE *fp)
{
    unsigned char LineBuffer[1024];
    ReadLineTrim(LineBuffer,fp);

    char tKey[100];
    int r = Strok(LineBuffer+1," <>",tKey)+1;
//    VtsTagParaser(tKey,LineBuffer+r,_vsf);
    int k = 0;
    for(;k<SALEC_VTS_TAG_TYPES;++k)
    {
        if(0==strcmp(tKey,TagName[k])) break;
    }
    if(k==SALEC_VTS_TAG_TYPES)
    {
        fprintf(stdout,"Undefined TagName:%s\n",tKey);
        exit(0);
    } else
    {
        TagNameP[k](LineBuffer+r,_vsf);
    }
    return _vsf->StackPos - 1;
}

int VtsCoordinateReshape(VtsInfo * _vsf)
{
    int k = 0;
    for(;k<_vsf->PointNoF;++k)
    {
        if(0== strcasecmp("coordinate",_vsf->PointField[k].Name)) break;
    }

    if(k==_vsf->PointNoF)
    {
        fprintf(stdout,"No coordinate information in the vts!\n");
    }
    VTSDATAFLOAT * PointData = _vsf->PointField[k].Data;
    unsigned long PointDataLen = _vsf->PointField[k].DataLen;


    unsigned int nx0 = _vsf->PieceExtent[0][1]-_vsf->PieceExtent[0][0]+1;
    unsigned int nx1 = _vsf->PieceExtent[1][1]-_vsf->PieceExtent[1][0]+1;
    unsigned int nx2 = _vsf->PieceExtent[2][1]-_vsf->PieceExtent[2][0]+1;

    if(nx0*nx1*nx2*VTSDIM!=PointDataLen)
    {
        fprintf(stdout,"number of coordinates is %ld, but %d is wanted\n",PointDataLen,nx1*nx0*nx2*VTSDIM);
        exit(0);
    }

    _vsf->Point = (VTSDATAFLOAT ****) malloc(sizeof(VTSDATAFLOAT ***)*nx0);
    for(int xi0=0;xi0<nx0;++xi0)
    {
        _vsf->Point[xi0] = (VTSDATAFLOAT ***) malloc(sizeof(VTSDATAFLOAT **)*nx1);
        for(int xi1=0;xi1<nx1;++xi1)
        {
            _vsf->Point[xi0][xi1] = (VTSDATAFLOAT **) malloc(sizeof(VTSDATAFLOAT *)*nx2);
            for(unsigned long xi2=0;xi2<nx2;++xi2)
            {
                _vsf->Point[xi0][xi1][xi2] = PointData + (xi2 + xi1*nx0 + xi0*nx0*nx1)*VTSDIM;
            }
        }
    }

    return (int) PointDataLen;
}

void VtsInfoClean(VtsInfo * _vsf)
{
    for(int k=0;k<_vsf->PointNoF;++k)
        free(_vsf->PointField[k].Data);
    free(_vsf->PointField);

    for(int k=0;k<_vsf->CellNoF;++k)
        free(_vsf->CellField[k].Data);
    free(_vsf->CellField);

    free(_vsf->VtsStack);
    unsigned int nx0 = _vsf->PieceExtent[0][1]-_vsf->PieceExtent[0][0]+1;
    unsigned int nx1 = _vsf->PieceExtent[1][1]-_vsf->PieceExtent[1][0]+1;
    unsigned int nx2 = _vsf->PieceExtent[2][1]-_vsf->PieceExtent[2][0]+1;

    for(int xi0=0;xi0<nx0;++xi0)
    {
        for(int xi1=0;xi1<nx1;++xi1)
        {
            free(_vsf->Point[xi0][xi1]);
        }
        free(_vsf->Point[xi0]);
    }
    free(_vsf->Point);
}

