//
// Created by li on 22-12-1.
//
/*
 * first edition is copy from citcoms,
 * https://github.com/geodynamics/citcoms
 * this edition is revised for sale2d_rebuild.
 */


#ifndef SALE_REBUILD_OUTPUT_H
#define SALE_REBUILD_OUTPUT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define USE_GZDIR 1

void write_ascii_array(int nn, int perLine, float *array, FILE *fp);
void zlibcompress(unsigned char* in, int nn, unsigned char** out, int *nn2);
void FloatToUnsignedChar(float * floatarray, int nn, unsigned char * chararray);
void base64plushead(unsigned char * in, int nn, int orinn, unsigned char* out);
void IntToUnsignedChar(int * intarray, int nn, unsigned char * chararray);
void base64(unsigned char * in, int nn, unsigned char* out);
void write_binary_array(int nn, float* array, FILE * f);

void vts_file_header(FILE *fp, const char * Pextent, const char * Wextent);
void vts_file_trailer(FILE *fp);
void vtk_point_data_header(FILE *fp);
void vtk_point_data_trailer(FILE *fp);
void vtk_cell_data_header(FILE *fp);
void vtk_cell_data_trailer(FILE *fp);
void pvtk_point_data_header(FILE *fp);
void pvtk_point_data_trailer(FILE *fp);
void pvtk_cell_data_header(FILE *fp);
void pvtk_cell_data_trailer(FILE *fp);
void vtk_dataarray(FILE *fp, const char * name, const char * vtk_format, const double * dataarray, int len_dataarray);
void vtk_dataarrayf(FILE *fp, const char * name, const char * vtk_format, const float * dataarray, int len_dataarray);
void vtk_dataarray_vec(FILE *fp, const char * name, const char * vtk_format, const double * dataarray, int len_dataarray, int n_component);
void vtk_dataarray_vecf(FILE *fp, const char * name, const char * vtk_format, const float * dataarray, int len_dataarray, int n_component);
void vtk_output_vec(FILE *fp, const char * array_name, const char *vtk_format, const double * dataarray, int len_dataarray);
void vtk_output_coord(FILE *fp, const char * vtk_format,const double * dataarray, int len_dataarray);
void vtk_output_coordf(FILE *fp, const char * vtk_format,const float * dataarray, int len_dataarray);
void vtk_output_coord2d(FILE *fp, const char * vtk_format,const double * dataarray, int len_dataarray);
void vtk_output_coord2df(FILE *fp, const char * vtk_format,const float * dataarray, int len_dataarray);

void vtp_file_header(FILE *fp, int num_pts);
void vtp_file_trailer(FILE *fp);
void vtk_point_header(FILE *fp);
void vtk_point_trailer(FILE *fp);
void vtk_point_data_header_with_attr(FILE *fp, const char * attr);
#endif //SALE_REBUILD_OUTPUT_H
