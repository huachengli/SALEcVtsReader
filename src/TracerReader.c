//
// Created by huacheng on 3/23/22.
//

#include "TracerReader.h"

void LoadTracerInfo(TracerInfo * _tinfo,const char * fname)
{
    InputFile * ifp = OpenInputFile(fname);
    int npgx,npgy,npgz;
    npgx = GetValueI(ifp,"processor.npgx","1");
    npgy = GetValueI(ifp,"processor.npgy","1");
    npgz = GetValueI(ifp,"processor.npgz","1");

    _tinfo->nproc = npgz*npgy*npgx;
    _tinfo->maxstep = GetValueI(ifp,"cycle.maxstep","1");
    GetValueS(ifp,"tracer.type",_tinfo->type,"NONE");
    GetValueS(ifp,"tracer.prefix",_tinfo->prefix,"bm");
    CloseInputFile(ifp);
}

void TracerName(TracerInfo * _tinfo,int rank,int step,char * _tname)
{
    // set as a marco ?
    sprintf(_tname,"%s/%s.proc%d.%04d.tracer",_tinfo->path,_tinfo->prefix,rank,step);
}

int LoadTracerFile(const char * fname, Tracer * _tracer)
{
    FILE * fp = fopen(fname,"r");
    if(NULL == fp)
    {
        _tracer->num = 0;
        return 0;
    }

    int TracerNum = 0;
    fscanf(fp,"    %d\n",&TracerNum);
    _tracer->num = TracerNum;
    _tracer->data = (TracerData *) malloc(sizeof(TracerData)*TracerNum);
    _tracer->node = (AVLNode *) malloc(sizeof(AVLNode)*TracerNum);

    for(int k=0;k<TracerNum;++k)
    {
        TracerData * pdata = _tracer->data + k;
        fscanf(fp,"%d,%d,%lf,%lf,%lf,%le,%lf\n",&(pdata->Id),&(pdata->eId),
               pdata->position+0,pdata->position+1,pdata->position+2,
               pdata->data+0,pdata->data+1);
        // attach the data to AVLNode
        _tracer->node[k].attach = (uintptr_t) pdata;
        BalanceInsert(&_tracer->root,_tracer->node+k,TracerIdCmp);
    }

    fclose(fp);
    return TracerNum;
}

Tracer * OpenTracerFile(const char * fname)
{
    Tracer * _tracer = (Tracer *) malloc(sizeof(Tracer));
    FILE * fp = fopen(fname,"r");
    if(NULL == fp)
    {
        _tracer->num = 0;
        free(_tracer);
        return NULL;
    }

    int TracerNum = 0;
    fscanf(fp,"    %d\n",&TracerNum);
    _tracer->num = TracerNum;
    _tracer->data = (TracerData *) malloc(sizeof(TracerData)*TracerNum);
    _tracer->node = (AVLNode *) malloc(sizeof(AVLNode)*TracerNum);
    _tracer->root = NULL;

    for(int k=0;k<TracerNum;++k)
    {
        TracerData * pdata = _tracer->data + k;
        fscanf(fp,"%d,%d,%lf,%lf,%lf,%le,%lf\n",&(pdata->Id),&(pdata->eId),
               pdata->position+0,pdata->position+1,pdata->position+2,
               pdata->data+0,pdata->data+1);
        // attach the data to AVLNode
        _tracer->node[k].attach = (uintptr_t) pdata;
        BalanceInsert(&_tracer->root,_tracer->node+k,TracerIdCmp);
    }

    fclose(fp);
    return _tracer;
}



void CloseTracerFile(Tracer * _tracer)
{
    free(_tracer->node);
    free(_tracer->data);
    free(_tracer);
}

double CompareTracerFile(Tracer * a,Tracer * b)
{
    int node_matched = 0;
    for(int k=0;k<b->num;++k)
    {
        if(NULL != AVLQuery(a->root,b->node+k,TracerIdCmp)) ++node_matched;
    }
    return 2.0*node_matched/(a->num + b->num);
}