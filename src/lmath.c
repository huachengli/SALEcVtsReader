//
// Created by huachengli on 9/1/24.
//

#include "lmath.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

void VecAddF(float *x, float *y, float p, int n)
{
    assert(n >= 1);
    for(int k=0;k<n;++k) x[k] += y[k]*p;
}

void VecZeroF(float *x, int n)
{
    assert(n >= 1);
    for(int k=0;k<n;++k) x[k] = 0.f;
}

void VecZero(double *x, int n)
{
    assert(n >= 1);
    for(int k=0;k<n;++k) x[k] = 0.;
}

int VecMaxArgF(const float *x, int n)
{
    assert(n >= 1);
    int rst = 0;
    for(int k=0;k<n;++k)
    {
        if(x[rst] < x[k]) rst = k;
    }
    return rst;
}

int VecMinArgF(const float *x, int n)
{
    assert(n >= 1);
    int rst = 0;
    for(int k=0;k<n;++k)
    {
        if(x[rst] > x[k]) rst = k;
    }
    return rst;
}

float VecMaxF(const float *x, int n)
{
    assert(n >= 1);
    float rst = x[0];
    for(int k=0;k<n;++k)
    {
        if(rst < x[k]) rst = x[k];
    }
    return rst;
}

float VecMinF(const float *x, int n)
{
    assert(n >= 1);
    float rst = x[0];
    for(int k=0;k<n;++k)
    {
        if(rst > x[k]) rst = x[k];
    }
    return rst;
}

float VecDisF(const float *x, const float * y, int n)
{
    assert(n>=1);
    float rst = 0;
    for(int k=0;k<n;++k)
    {
        rst += (x[k] - y[k])*(x[k] - y[k]);
    }
    return sqrtf(rst);
}

double VecDis(const double *x, const double *y, int n)
{
    assert(n>=1);
    double rst = 0;
    for(int k=0;k<n;++k)
    {
        rst += (x[k] - y[k])*(x[k] - y[k]);
    }
    return sqrt(rst);
}

double VecScaler(double *x, double p,int n)
{
    for(int k=0;k<n;++k) x[k] *= p;
}

double VecLen(const double *x, int n)
{
    double rst = 0.;
    for(int k=0; k<n; ++k) rst += x[k]*x[k];
    return sqrt(rst);
}

float VecLenF(const float *x, int n)
{
    float rst = 0.f;
    for(int k=0; k<n; ++k) rst += x[k]*x[k];
    return sqrtf(rst);
}

double VecDot(const double *x, const double *y, int n)
{
    double rst = 0.;
    for(int k=0; k<n; ++k) rst += x[k]*y[k];
    return rst;
}

float VecDotF(const float *x, const float *y, int n)
{
    float rst = 0.f;
    for(int k=0; k<n; ++k) rst += x[k]*y[k];
    return rst;
}

void VecCross(double *z, const double *x, const double *y, int n) {
    assert(3 == n);
    z[0] = x[1] * y[2] - x[2] * y[1];
    z[1] = x[2] * y[0] - x[0] * y[2];
    z[2] = x[0] * y[1] - x[1] * y[0];
}

void VecLinear(double *z, const double *x, double px, const double *y, double py, int n)
{
    for(int k=0;k<n;++k) z[k] = x[k]*px + y[k]*py;
}

void VecAdd(double *x , const double *y, double p,int n)
{
    assert(n >= 1);
    for(int k=0;k<n;++k) x[k] += y[k]*p;
}

void VecScale(double *x, double px, int n)
{
    assert(n>=1 && x!=NULL);
    for(int k=0;k<n;++k) x[k] *= px;
}

void VecCopy(double * y, double * x, int n)
{
    for(int k=0;k<n;++k)
    {
        y[k] = x[k];
    }
}

void VecNormalize(double *x, int n)
{
    double Ln = VecLen(x,n);
    VecScaler(x, 1.0/Ln, n);
}


void VecD2F(float *y, const double * x, int n)
{
    for(int k=0;k<n;++k)
        y[k] = (float) x[k];
}

void VecF2D(double *y, const float * x, int n)
{
    for(int k=0;k<n;++k)
        y[k] = (float) x[k];
}

void VecRotate(double * x, double ro, double fo, int n)
{
    assert(n==3);

    double Rm[4][4], t[3];

    Rm[1][1] = cos(ro) * cos(fo);
    Rm[1][2] = cos(ro) * sin(fo);
    Rm[1][3] = -sin(ro);
    Rm[2][1] = -sin(fo);
    Rm[2][2] = cos(fo);
    Rm[2][3] = 0.0;
    Rm[3][1] = sin(ro) * cos(fo);
    Rm[3][2] = sin(ro) * sin(fo);
    Rm[3][3] = cos(ro);

    for(int k=0;k<3;++k)
    {
        t[k] = 0;
        for(int j=0;j<3;++j)
        {
            t[k] += Rm[k+1][j+1]*x[j];
        }
    }
    for(int k=0;k<3;++k)
        x[k] = t[k];
}


double Clock(int i){
    static struct timeval start = {0, 0};

    if(i == 0)
    {
        gettimeofday(&start, NULL);
        return 0.;
    }
    else
    {
        struct timeval current;
        gettimeofday(&current,NULL);
        double elapsed = ( - start.tv_sec + current.tv_sec) + ( - start.tv_usec + current.tv_usec)/1000000.0;
        return elapsed;
    }
}

/*
   Derivation from the fortran version of CONREC by Paul Bourke
   d               ! matrix of data to contour
   ilb,iub,jlb,jub ! index bounds of data matrix
   x               ! data matrix column coordinates
   y               ! data matrix row coordinates
   nc              ! number of contour levels
   z               ! contour levels in increasing order
   COPY　from paulbourke.net/papers/conrec/conrec.c

    EDITED by huachengli 20/5/2025
    change (nc,z) to single level (z)
    change func ConrecLine to (con_pts, num_pts)
*/
void Contour(double **d, int ilb, int iub, int jlb, int jub, double *x, double *y, double z,
             double ** pcon_pts, int * num_pts)
{
    #define xsect(p1, p2) (h[p2]*xh[p1]-h[p1]*xh[p2])/(h[p2]-h[p1])
    #define ysect(p1, p2) (h[p2]*yh[p1]-h[p1]*yh[p2])/(h[p2]-h[p1])

    int m1, m2, m3, case_value;
    double x1 = 0, x2 = 0, y1 = 0, y2 = 0;
    int i, j, k, m;
    double h[5];
    int sh[5];
    double xh[5], yh[5];
    int im[4] = {0, 1, 1, 0}, jm[4] = {0, 0, 1, 1};
    const int castab[3][3][3] = {
            {{0, 0, 8}, {0, 2, 5}, {7, 6, 9}},
            {{0, 3, 4}, {1, 3, 1}, {4, 3, 0}},
            {{9, 6, 7}, {5, 2, 0}, {8, 0, 0}}
    };

    // check potential contour cells
    int num_con_cells = 0;
    for(j = (jub - 1); j >= jlb; j--)
    {
        for(i = ilb; i <= iub - 1; i++)
        {
            // check z-level in local rect (i,j)->(i+1,j+1)
            double temp1 = Min(d[i][j], d[i][j + 1]);
            double temp2 = Min(d[i + 1][j], d[i + 1][j + 1]);
            double dminl = Min(temp1, temp2);
            temp1 = Max(d[i][j], d[i][j + 1]);
            temp2 = Max(d[i + 1][j], d[i + 1][j + 1]);
            double dmaxl = Max(temp1, temp2);

            if(z < dminl || z > dmaxl)
                continue;
            num_con_cells ++;
        }
    }

    // alloc mem for con_pts
    int num_line_alloc = num_con_cells + 8;
    double * con_pts = malloc(sizeof(double)*4*num_line_alloc);
    int k_line = 0;

    for(j = (jub - 1); j >= jlb; j--)
    {
        for(i = ilb; i <= iub - 1; i++)
        {
            // check z-level in local rect (i,j)->(i+1,j+1)
            double temp1 = Min(d[i][j], d[i][j + 1]);
            double temp2 = Min(d[i + 1][j], d[i + 1][j + 1]);
            double dminl = Min(temp1, temp2);
            temp1 = Max(d[i][j], d[i][j + 1]);
            temp2 = Max(d[i + 1][j], d[i + 1][j + 1]);
            double dmaxl = Max(temp1, temp2);

            if(z < dminl || z > dmaxl)
                continue;

            for(m = 4; m >= 0; m--)
            {
                if(m > 0)
                {
                    h[m] = d[i + im[m - 1]][j + jm[m - 1]] - z;
                    xh[m] = x[i + im[m - 1]];
                    yh[m] = y[j + jm[m - 1]];
                }
                else
                {
                    h[0] = 0.25 * (h[1] + h[2] + h[3] + h[4]);
                    xh[0] = 0.50 * (x[i] + x[i + 1]);
                    yh[0] = 0.50 * (y[j] + y[j + 1]);
                }
                if(h[m] > 0.0)
                    sh[m] = 1;
                else if(h[m] < 0.0)
                    sh[m] = -1;
                else
                    sh[m] = 0;
            }

            /*
               Note: at this stage the relative heights of the corners and the
               centre are in the h array, and the corresponding coordinates are
               in the xh and yh arrays. The centre of the box is indexed by 0
               and the 4 corners by 1 to 4 as shown below.
               Each triangle is then indexed by the parameter m, and the 3
               vertices of each triangle are indexed by parameters m1,m2,and m3.
               It is assumed that the centre of the box is always vertex 2
               though this is important only when all 3 vertices lie exactly on
               the same contour level, in which case only the side of the box
               is drawn.
                  vertex 4 +-------------------+ vertex 3
                           | \               / |
                           |   \    m-3    /   |
                           |     \       /     |
                           |       \   /       |
                           |  m=2    X   m=2   |       the centre is vertex 0
                           |       /   \       |
                           |     /       \     |
                           |   /    m=1    \   |
                           | /               \ |
                  vertex 1 +-------------------+ vertex 2
            */
            /* Scan each triangle in the box */
            for(m = 1; m <= 4; m++)
            {
                m1 = m;
                m2 = 0;
                if(m != 4)
                    m3 = m + 1;
                else
                    m3 = 1;
                if((case_value = castab[sh[m1] + 1][sh[m2] + 1][sh[m3] + 1]) == 0)
                    continue;
                switch(case_value)
                {
                    case 1: /* Line between vertices 1 and 2 */
                        x1 = xh[m1];
                        y1 = yh[m1];
                        x2 = xh[m2];
                        y2 = yh[m2];
                        break;
                    case 2: /* Line between vertices 2 and 3 */
                        x1 = xh[m2];
                        y1 = yh[m2];
                        x2 = xh[m3];
                        y2 = yh[m3];
                        break;
                    case 3: /* Line between vertices 3 and 1 */
                        x1 = xh[m3];
                        y1 = yh[m3];
                        x2 = xh[m1];
                        y2 = yh[m1];
                        break;
                    case 4: /* Line between vertex 1 and side 2-3 */
                        x1 = xh[m1];
                        y1 = yh[m1];
                        x2 = xsect(m2, m3);
                        y2 = ysect(m2, m3);
                        break;
                    case 5: /* Line between vertex 2 and side 3-1 */
                        x1 = xh[m2];
                        y1 = yh[m2];
                        x2 = xsect(m3, m1);
                        y2 = ysect(m3, m1);
                        break;
                    case 6: /* Line between vertex 3 and side 1-2 */
                        x1 = xh[m3];
                        y1 = yh[m3];
                        x2 = xsect(m1, m2);
                        y2 = ysect(m1, m2);
                        break;
                    case 7: /* Line between sides 1-2 and 2-3 */
                        x1 = xsect(m1, m2);
                        y1 = ysect(m1, m2);
                        x2 = xsect(m2, m3);
                        y2 = ysect(m2, m3);
                        break;
                    case 8: /* Line between sides 2-3 and 3-1 */
                        x1 = xsect(m2, m3);
                        y1 = ysect(m2, m3);
                        x2 = xsect(m3, m1);
                        y2 = ysect(m3, m1);
                        break;
                    case 9: /* Line between sides 3-1 and 1-2 */
                        x1 = xsect(m3, m1);
                        y1 = ysect(m3, m1);
                        x2 = xsect(m1, m2);
                        y2 = ysect(m1, m2);
                        break;
                    default:
                        break;
                }

                /* Finally draw the line */
                // ConrecLine(x1, y1, x2, y2, z);
                // store line into con_pts

                if(k_line+1 >= num_line_alloc - 1)
                {
                    con_pts = realloc(con_pts, 2*num_line_alloc*sizeof(double)*4);
                    num_line_alloc = 2*num_line_alloc;
                }

                con_pts[4*k_line + 0] = x1;
                con_pts[4*k_line + 1] = y1;
                con_pts[4*k_line + 2] = x2;
                con_pts[4*k_line + 3] = y2;
                k_line += 1;
            } /* m */
        } /* i */
    } /* j */

    num_pts[0] = 2*k_line;
    num_pts[1] = 2*num_line_alloc;
    *pcon_pts = con_pts;
}

/*
 * sort and connect lines generated in Contour
 */
void Contour_pts_sort(double * con_pts, int * num_pts, int ** pcon_seg, double tol)
{
    // get minimum length and upper limit of tol
    double mtol = fabs(con_pts[2] - con_pts[0]) + fabs(con_pts[3] - con_pts[1]);
    for(int k=0;k<num_pts[0]/2;++k)
    {
        double k_dist = fabs(con_pts[4*k] - con_pts[4*k+2]) + fabs(con_pts[4*k+1] - con_pts[4*k+3]);
        for(int j=k+1;j<num_pts[0];++j)
        {
            if(k_dist < mtol)
            {
               mtol = k_dist;
            }
        }
    }
    if(tol <= 0)
    {
        tol = (fabs(tol) + 1e-6)*mtol;
    }
    else
    {
        tol = Min(tol, 0.5*mtol);
    }

    int * nei_tab = malloc(sizeof(int) * num_pts[0]);
    for(int k=0;k<num_pts[0];++k)
    {
        nei_tab[k] = -1;
    }

    // find pairs of distance < tol
    for(int k=0;k<num_pts[0];++k)
    {
        if(nei_tab[k] >= 0)
            continue;
        for(int j=k+1;j<num_pts[0];++j)
        {
            double kj_dist = fabs(con_pts[2*k] - con_pts[2*j]) + fabs(con_pts[2*k+1] - con_pts[2*j+1]);
            if(kj_dist < tol)
            {
                nei_tab[k] = j;
                nei_tab[j] = k;
                break;
            }
        }
    }

    int * contour_id = malloc(sizeof(int) * num_pts[0]);
    int * poly_record = malloc(sizeof(int) * num_pts[0]);
    int k_poly_record = 0;

    for(int k=0;k<num_pts[0];++k)
    {
        contour_id[k] = -1;
        poly_record[k] = -1;
    }

    int cur_id = 1;
    for(int k=0;k<num_pts[0]/2;++k)
    {
        int head = 2*k;
        int tail = 2*k+1;
        if(contour_id[head] >= 0 || contour_id[tail] >= 0)
            continue;
        contour_id[head] = contour_id[tail] = cur_id;
        int len_poly = 2;
        while(1)
        {
            if(nei_tab[tail] == -1)
                break;
            int p0 = nei_tab[tail];
            int p1 = p0%2==0? p0+1: p0-1;

            if(contour_id[p0] == -1 && contour_id[p1] == -1)
            {
                contour_id[p0] = contour_id[p1] = contour_id[tail];
                tail = p1;
                len_poly++;
            }
            else
            {
                break;
            }
        }

        while(1)
        {
            if(nei_tab[head] == -1)
                break;
            int p0 = nei_tab[head];
            int p1 = p0%2==0? p0+1: p0-1;
            if(contour_id[p0] == -1 && contour_id[p1] == -1)
            {
                contour_id[p0] = contour_id[p1] = contour_id[head];
                head = p1;
                len_poly++;
            }
            else
            {
                break;
            }
        }
        cur_id++;

        poly_record[k_poly_record] = len_poly;
        poly_record[k_poly_record + 1] = head;
        int poly_tail = head;
        for(int j=2;j<=len_poly;++j)
        {
            int p0 = poly_tail;
            if(p0==-1)
            {
                poly_record[k_poly_record + j] = -1;
                continue;
            }

            int p1 = p0%2==0? p0+1: p0-1;
            poly_record[k_poly_record + j] = p1;
            poly_tail = nei_tab[p1];
        }
        k_poly_record += 1 + len_poly;
    }

    int * con_pts_seg = malloc(sizeof(int)*cur_id);
    int k_con_pts_seg = 1;
    double * con_pts_sort = malloc(sizeof(double)*num_pts[0]*2);
    int k_con_pts_sort = 0;

    int k = 0;
    while(k < num_pts[0])
    {
        int poly_counter = poly_record[k];
        if(-1==poly_counter)
            break;
        k++;

        con_pts_seg[k_con_pts_seg++] = poly_counter;
        int ph = poly_record[k]; // head
        int pt = poly_record[k + poly_counter - 1]; //tail
        int reversed = 0;
        if(ph>=0 && pt >= 0 && con_pts[2*ph] > con_pts[2*pt])
            reversed = 1;

        for(int j=0;j<poly_counter;++j)
        {
            int p0 = poly_record[reversed?k+poly_counter-1-j:k+j];
            if(-1 == p0)
            {
                con_pts_sort[k_con_pts_sort*2 + 0] = con_pts_sort[k_con_pts_sort*2 + 1] = 0.0;
            }
            else
            {
                con_pts_sort[k_con_pts_sort*2 + 0] = con_pts[2*p0 + 0];
                con_pts_sort[k_con_pts_sort*2 + 1] = con_pts[2*p0 + 1];
            }
            k_con_pts_sort++;
        }
        k += poly_counter;
    }

    memcpy(con_pts,con_pts_sort,k_con_pts_sort* sizeof(double)*2);

    con_pts_seg[0] = k_con_pts_seg - 1;

    fputs("\n",stdout);
    for(int k=0;k<k_con_pts_seg;++k)
    {
        fprintf(stdout,"seg[%3d]=%d\n",k,con_pts_seg[k]);
    }

    num_pts[0] = k_con_pts_sort;
    *pcon_seg = con_pts_seg;
    free(con_pts_sort);
    free(contour_id);
    free(poly_record);
    free(nei_tab);
}
