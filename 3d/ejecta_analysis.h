//
// Created by li on 8/20/24.
//

#ifndef SALECVTSREADER_EJECTA_ANALYSIS_H
#define SALECVTSREADER_EJECTA_ANALYSIS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include "InputParser.h"
#include "VtkWriter.h"


typedef struct EjectaImpl
{
    double t;
    double pos[3];
    double init_pos[3];
    double land_pos[3];
    double vel[3];
    double maxpre;
    double maxtem;
    double pre;
    double tem;

    int matid;
    int rank;
    int id;
    int NeE;
} ejecta_t;

typedef struct EjectaCollect{
    int cur_step;
    int min_step;
    int max_step;
    int min_step_e;
    int max_step_e;
    int nproc;
    char prefix[4096];
    char output[4096];
    ejecta_t * data;
    int len;
    int len_allocated;
} ejecta_collect;

int ejecta_collect_test_init(ejecta_collect * _ec);
int ejecta_collect_init(ejecta_collect * _ec, InputFile * ifp);
void ejecta_collect_test_clean(ejecta_collect * _ec);
int load_ejecta_collect(ejecta_collect * _ec, int step);
int load_ejecta_collect_single_file(ejecta_collect * _ec, const char * _tmp_name);
int ejecta_collect_push(ejecta_collect * _ec, ejecta_t * _e);
int ejecta_collect_to_vtp(ejecta_collect * _ec, const char * vtp_name);
void numerical_ejecta_orbit_moon(ejecta_collect * _ec, double dt);

#endif //SALECVTSREADER_EJECTA_ANALYSIS_H
