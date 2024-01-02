//
// Created by li on 22-12-1.
//

#include "VtkWriter.h"

#ifdef USE_GZDIR
#include <zlib.h>
#include <assert.h>
#define CHUNK 16384
#endif
void vts_file_header(FILE *fp, const char * Pextent, const char * Wextent)
{
    const char format[] =
            "<?xml version=\"1.0\"?>\n"
            "<VTKFile type=\"StructuredGrid\" version=\"0.1\" compressor=\"vtkZLibDataCompressor\" byte_order=\"LittleEndian\">\n"
            "  <StructuredGrid WholeExtent=\"%s\">\n"
            "    <Piece Extent=\"%s\">\n";

    char header[1024];

    snprintf(header, 1024, format, Wextent, Pextent);
    fputs(header, fp);
}

void vts_file_trailer(FILE *fp)
{
    const char trailer[] =
            "    </Piece>\n"
            "  </StructuredGrid>\n"
            "</VTKFile>\n";

    fputs(trailer, fp);
}

void vtk_point_data_trailer(FILE *fp)
{
    fputs("      </PointData>\n", fp);
}

void vtk_cell_data_header(FILE *fp)
{
    fputs("      <CellData>\n", fp);
}

void vtk_cell_data_trailer(FILE *fp)
{
    fputs("      </CellData>\n", fp);
}

void pvtk_point_data_header(FILE *fp)
{
    fputs("      <PPointData>\n", fp);
}

void pvtk_point_data_trailer(FILE *fp)
{
    fputs("      </PPointData>\n", fp);
}

void pvtk_cell_data_header(FILE *fp)
{
    fputs("      <PCellData>\n", fp);
}

void pvtk_cell_data_trailer(FILE *fp)
{
    fputs("      </PCellData>\n", fp);
}


void write_ascii_array(int nn, int perLine, float *array, FILE *fp)
{
    int i;
    switch(perLine) {
        case 1:
            for(i=0; i<nn; i++)
                fprintf(fp, "%.4e\n", array[i]);
            break;
        case 2:
            for(i=0; i < nn/2; i++)
                fprintf(fp,"%.4e %.4e\n",array[2*i],array[2*i+1]);
            break;
        case 3:
            for(i=0; i < nn/3; i++)
                fprintf(fp,"%.4e %.4e %.4e\n",array[3*i],array[3*i+1],array[3*i+2]);
            break;
        case 6:
            for(i=0; i < nn/6; i++)
                fprintf(fp,"%.4e %.4e %.4e %.4e %.4e %.4e\n",
                        array[6*i],array[6*i+1],array[6*i+2],
                        array[6*i+3],array[6*i+4],array[6*i+5]);
            break;
        default:
            break;
    }
}

void zlibcompress(unsigned char* in, int nn, unsigned char** out, int *nn2)
/* function to compress "in" to "out" reducing size from nn to nn2 */
{
#ifdef USE_GZDIR
    int ntemp=0;

    /* in and out of z-stream */
    unsigned char inz[CHUNK];
    unsigned char outz[CHUNK];

    /* compression level */
    int level = Z_NO_COMPRESSION; //Z_BEST_SPEED;// Z_DEFAULT_COMPRESSION;
    int ret,flush;
    int i,j,k;

    /* zlib compression stream */
    z_stream strm;

    /* hope compressed data will be <= uncompressed */
    *out = malloc(sizeof(unsigned char)*(nn*1.02+64));

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
#endif
    return;
}

void FloatToUnsignedChar(float * floatarray, int nn, unsigned char * chararray)
{
    /* simple float to unsigned chararray routine via union
    nn=length(intarray) chararray has to be BIG ENOUGH! */
    int i;
    union FloatToUnsignedChars
    {
        float input;
        unsigned char output[4];
    } floattransform;

    for (i=0; i<nn; i++){
        floattransform.input=floatarray[i];
        chararray[4*i]=floattransform.output[0];
        chararray[4*i+1]=floattransform.output[1];
        chararray[4*i+2]=floattransform.output[2];
        chararray[4*i+3]=floattransform.output[3];
    }
}

void base64plushead(unsigned char * in, int nn, int orinn, unsigned char* out)
{
    /* writing vtk compatible zlib compressed base64 encoded data to "out" */
    int i;
    unsigned char * b64head;
    int b64bodylength;
    unsigned char * b64body;
    /* header of data */
    unsigned char * charhead = malloc(sizeof(unsigned char)*16);
    /* - consists of "1" (number of pieces) */
    /* - original datalength in byte */
    /* - original datalength in byte */
    /* - new datalength after z-lib compression */
    int * headInts= malloc(sizeof(int)*4);
    headInts[0]=1;
    headInts[1]=orinn;
    headInts[2]=orinn;
    headInts[3]=nn;
    // transform to unsigned char
    IntToUnsignedChar(headInts,4,charhead);

    // base64: 16byte -> 24byte
    b64head =  malloc(sizeof(unsigned char)*24);
    // fills b64head
    base64(charhead, 16, b64head);

    // base64 data
    b64bodylength = 4*ceil((double) nn/3.0);
    b64body = malloc(sizeof(unsigned char)*b64bodylength);
    // writes base64 data to b64body
    base64(in,nn,b64body);

    // combines header and body
    for (i=0; i<24 ; i++){
        out[i]=b64head[i];
    }

    for (i=0; i<b64bodylength ; i++){
        out[24+i]=b64body[i];
    }

    if(b64body){free(b64body);}
    if(b64head){free(b64head);}
    if(headInts){free(headInts);}
    if(charhead){free(charhead);}
}

void IntToUnsignedChar(int * intarray, int nn, unsigned char * chararray)
{
    /* simple int - to unsigned chararray routine via union
    nn=length(intarray) chararray has to be BIG ENOUGH! */
    int i;
    union IntToUnsignedChars
    {
        int input;
        unsigned char output[4];
    } inttransform;

    for (i=0; i<nn; i++){
        inttransform.input=intarray[i];
        chararray[4*i]=inttransform.output[0];
        chararray[4*i+1]=inttransform.output[1];
        chararray[4*i+2]=inttransform.output[2];
        chararray[4*i+3]=inttransform.output[3];
    }
}

void base64(unsigned char * in, int nn, unsigned char* out)
{
    /*takes *in*-array and "in"-length-"nn" and fills "out"-array
    with base64(in) "out" needs to be big enough!!!
    length(out) >= 4* |^ nn/3.0 ^| */
    char cb64[]="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    int len;
    int i;

    for (i=0; i < nn; i+=3){

        len = (3 < nn-i ? 3 : nn-i);
        if (len >= 3){
            /* normal base64 encoding */
            out[i/3*4+0] = cb64[ in[i] >> 2 ];
            out[i/3*4+1] = cb64[ ((in[i] & 0x03) << 4) | ((in[i+1] & 0xf0) >> 4) ];
            out[i/3*4+2] = cb64[ ((in[i+1] & 0x0f) << 2) | ((in[i+2] & 0xc0) >> 6)];
            out[i/3*4+3] = cb64[ in[i+2] & 0x3f ];
        } else if (len == 2){
            /* at the end of array fill up with '=' */
            out[i/3*4+0] = cb64[ in[i] >> 2 ];
            out[i/3*4+1] = cb64[ ((in[i] & 0x03) << 4) | ((in[i+1] & 0xf0) >> 4) ];
            out[i/3*4+2] = cb64[((in[i+1] & 0x0f) << 2)];
            out[i/3*4+3] = (unsigned char) '=';
        } else if (len == 1){
            /* at the end of array fill up with '=' */
            out[i/3*4+0] = cb64[ in[i] >> 2 ];
            out[i/3*4+1] = cb64[ ((in[i] & 0x03) << 4) ];
            out[i/3*4+2] = (unsigned char) '=';
            out[i/3*4+3] = (unsigned char) '=';
        }
    }
}


void write_binary_array(int nn, float* array, FILE * f)
{
    /* writes vtk-data array of floats and performs zip and base64 encoding */
    int chararraylength=4*nn;	/* nn floats -> 4*nn unsigned chars */
    unsigned char * chararray = malloc (chararraylength * sizeof(unsigned char));
    int compressedarraylength = 0;
    unsigned char * compressedarray;
    unsigned char ** pointertocompressedarray= &compressedarray;
    int base64plusheadlength;
    unsigned char * base64plusheadarray;

    FloatToUnsignedChar(array,nn,chararray);

    /* compression routine */
    zlibcompress(chararray,chararraylength,pointertocompressedarray,&compressedarraylength);

    /* special header for zip compressed and bas64 encoded data
    header needs 4 int32 = 16 byte -> 24 byte due to base64 (4*16/3) */
    base64plusheadlength = 24 + 4*ceil((double) compressedarraylength/3.0);
    base64plusheadarray  = malloc(sizeof(unsigned char)* base64plusheadlength);

    /* fills base64plusheadarray with everything ready for simple writing */
    base64plushead(compressedarray,compressedarraylength, chararraylength, base64plusheadarray);

    fwrite(base64plusheadarray,sizeof(unsigned char),base64plusheadlength,f);
    fprintf(f,"\n");
    free(chararray);
    free(base64plusheadarray);
    free(compressedarray);
}

void vtk_dataarray(FILE *fp, const char * name, const char * vtk_format, const double * dataarray, int len_dataarray)
{
    float* floattemp = (float *)malloc(len_dataarray*sizeof(float));
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"%s\">\n",name,vtk_format);
    for(int i=0;i < len_dataarray;i++)
        floattemp[i] =  (float) dataarray[i];
    if(0==strcmp(vtk_format,"binary"))
        write_binary_array(len_dataarray,floattemp,fp);
    else
        write_ascii_array(len_dataarray,1,floattemp,fp);
    fputs("        </DataArray>\n", fp);
    free(floattemp);
}

void vtk_dataarrayf(FILE *fp, const char * name, const char * vtk_format, const float * dataarray, int len_dataarray)
{
    float* floattemp = (float *)malloc(len_dataarray*sizeof(float));
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"%s\">\n",name,vtk_format);
    for(int i=0;i < len_dataarray;i++)
        floattemp[i] =  (float) dataarray[i];
    if(0==strcmp(vtk_format,"binary"))
        write_binary_array(len_dataarray,floattemp,fp);
    else
        write_ascii_array(len_dataarray,1,floattemp,fp);
    fputs("        </DataArray>\n", fp);
    free(floattemp);
}


void vtk_dataarray_vec(FILE *fp, const char * name, const char * vtk_format, const double * dataarray, int len_dataarray, int n_component)
{
    float* floattemp = (float *)malloc(len_dataarray*sizeof(float)*n_component);
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"%s\">\n", name,n_component,vtk_format);
    for(int i=0;i < len_dataarray*n_component;i++)
        floattemp[i] =  (float) dataarray[i];

    if(0==strcmp(vtk_format,"binary"))
        write_binary_array(len_dataarray*n_component,floattemp,fp);
    else
        write_ascii_array(len_dataarray*n_component,n_component,floattemp,fp);
    fputs("        </DataArray>\n", fp);
    free(floattemp);
}

void vtk_dataarray_vecf(FILE *fp, const char * name, const char * vtk_format, const float * dataarray, int len_dataarray, int n_component)
{
    float* floattemp = (float *)malloc(len_dataarray*sizeof(float)*n_component);
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"%s\">\n", name,n_component,vtk_format);
    for(int i=0;i < len_dataarray*n_component;i++)
        floattemp[i] =  (float) dataarray[i];

    if(0==strcmp(vtk_format,"binary"))
        write_binary_array(len_dataarray*n_component,floattemp,fp);
    else
        write_ascii_array(len_dataarray*n_component,n_component,floattemp,fp);
    fputs("        </DataArray>\n", fp);
    free(floattemp);
}

void vtk_output_vec(FILE *fp, const char * array_name, const char *vtk_format, const double * dataarray, int len_dataarray)
{
    int n_component = 3;
    int nodes = len_dataarray;
    float* floatpos = (float*) malloc(nodes*n_component*sizeof(float));

    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"%s\">\n",
            array_name,n_component,vtk_format);

    for(int i=0;i<nodes;i++)
    {
        floatpos[3*i    ] = (float) dataarray[2*i];
        floatpos[3*i + 1] = (float) dataarray[2*i + 1];
        floatpos[3*i + 2] = 0.0f;
    }

    if(0==strcmp(vtk_format,"binary"))
        write_binary_array(nodes*n_component,floatpos,fp);
    else
        write_ascii_array(nodes*n_component,n_component,floatpos,fp);
    fputs("        </DataArray>\n", fp);
    free(floatpos);
}


void vtk_output_coord2d(FILE *fp, const char * vtk_format,const double * dataarray, int len_dataarray)
{
    /*
     * Output Cartesian coordinates as most VTK visualization softwares
     * assume it.
     */

    int n_component = 3;
    int nodes = len_dataarray;
    float* floatpos = (float*) malloc(nodes*n_component*sizeof(float));

    fputs("      <Points>\n", fp);
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"coordinate\" NumberOfComponents=\"%d\" format=\"%s\">\n", n_component,vtk_format);

    for(int i=0;i<nodes;i++)
    {
        floatpos[3*i    ] = (float) dataarray[2*i];
        floatpos[3*i + 1] = (float) dataarray[2*i + 1];
        floatpos[3*i + 2] = 0.0f;
    }

    if(0==strcmp(vtk_format,"binary"))
        write_binary_array(nodes*n_component,floatpos,fp);
    else
        write_ascii_array(nodes*n_component,n_component,floatpos,fp);
    fputs("        </DataArray>\n", fp);
    fputs("      </Points>\n", fp);
    free(floatpos);
}

void vtk_output_coord2df(FILE *fp, const char * vtk_format,const float * dataarray, int len_dataarray)
{
    /*
     * Output Cartesian coordinates as most VTK visualization softwares
     * assume it.
     */

    int n_component = 3;
    int nodes = len_dataarray;
    float* floatpos = (float*) malloc(nodes*n_component*sizeof(float));

    fputs("      <Points>\n", fp);
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"coordinate\" NumberOfComponents=\"%d\" format=\"%s\">\n", n_component,vtk_format);

    for(int i=0;i<nodes;i++)
    {
        floatpos[3*i    ] = (float) dataarray[2*i];
        floatpos[3*i + 1] = (float) dataarray[2*i + 1];
        floatpos[3*i + 2] = 0.0f;
    }

    if(0==strcmp(vtk_format,"binary"))
        write_binary_array(nodes*n_component,floatpos,fp);
    else
        write_ascii_array(nodes*n_component,n_component,floatpos,fp);
    fputs("        </DataArray>\n", fp);
    fputs("      </Points>\n", fp);
    free(floatpos);
}

void vtk_output_coord(FILE *fp, const char * vtk_format,const double * dataarray, int len_dataarray)
{
    /*
     * Output Cartesian coordinates as most VTK visualization softwares
     * assume it.
     */

    int n_component = 3;
    int nodes = len_dataarray;
    float* floatpos = (float*) malloc(nodes*n_component*sizeof(float));

    fputs("      <Points>\n", fp);
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"coordinate\" NumberOfComponents=\"%d\" format=\"%s\">\n", n_component,vtk_format);

    for(int i=0;i<nodes;i++)
    {
        floatpos[3*i    ] = (float) dataarray[3*i];
        floatpos[3*i + 1] = (float) dataarray[3*i + 1];
        floatpos[3*i + 2] = (float) dataarray[3*i + 2];
    }

    if(0==strcmp(vtk_format,"binary"))
        write_binary_array(nodes*n_component,floatpos,fp);
    else
        write_ascii_array(nodes*n_component,n_component,floatpos,fp);
    fputs("        </DataArray>\n", fp);
    fputs("      </Points>\n", fp);
    free(floatpos);
}

void vtk_output_coordf(FILE *fp, const char * vtk_format,const float * dataarray, int len_dataarray)
{
    /*
     * Output Cartesian coordinates as most VTK visualization softwares
     * assume it.
     */

    int n_component = 3;
    int nodes = len_dataarray;
    float* floatpos = (float*) malloc(nodes*n_component*sizeof(float));

    fputs("      <Points>\n", fp);
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"coordinate\" NumberOfComponents=\"%d\" format=\"%s\">\n", n_component,vtk_format);

    for(int i=0;i<nodes;i++)
    {
        floatpos[3*i    ] = (float) dataarray[3*i];
        floatpos[3*i + 1] = (float) dataarray[3*i + 1];
        floatpos[3*i + 2] = (float) dataarray[3*i + 2];
    }

    if(0==strcmp(vtk_format,"binary"))
        write_binary_array(nodes*n_component,floatpos,fp);
    else
        write_ascii_array(nodes*n_component,n_component,floatpos,fp);
    fputs("        </DataArray>\n", fp);
    fputs("      </Points>\n", fp);
    free(floatpos);
}


void vtp_file_header(FILE *fp, int num_pts)
{
    const char format[] =
            "<?xml version=\"1.0\"?>\n"
            "<VTKFile type=\"PolyData\" version=\"1.0\" compressor=\"vtkZLibDataCompressor\" byte_order=\"LittleEndian\">\n"
            "  <PolyData>\n"
            "    <Piece NumberOfPoints=\"%d\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";

    char header[1024];

    snprintf(header, 1024, format, num_pts);
    fputs(header, fp);
}

void vtp_file_trailer(FILE *fp)
{
    const char trailer[] =
            "    </Piece>\n"
            "  </PolyData>\n"
            "</VTKFile>\n";

    fputs(trailer, fp);
}

void vtk_point_data_header_with_attr(FILE *fp, const char * attr)
{
    fprintf(fp,"      <PointData %s >\n",attr);
}