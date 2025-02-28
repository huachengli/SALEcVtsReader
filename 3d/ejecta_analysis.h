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
#include "lmath.h"
#include "citcoms_related.h"

typedef struct EjectaImpl
{
    double t;
    double theta;
    double a;
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
    double v0;
    char prefix[4096];
    char output[4096];
    double R;
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
double approximate_ejecta(double *x, double *v, double R, double g0);
void analytical_ejecta_orbit_moon(ejecta_collect * _ec, double R, double g0);
void calculate_ejecta_thickness(citcoms_sphere * _cs, ejecta_collect * _ec, double R);

#endif //SALECVTSREADER_EJECTA_ANALYSIS_H
