/*
 * This file is part of OpenModelica.
 *
 * Copyright (c) 1998-CurrentYear, Linköping University,
 * Department of Computer and Information Science,
 * SE-58183 Linköping, Sweden.
 *
 * All rights reserved.
 *
 * THIS PROGRAM IS PROVIDED UNDER THE TERMS OF GPL VERSION 3
 * AND THIS OSMC PUBLIC LICENSE (OSMC-PL).
 * ANY USE, REPRODUCTION OR DISTRIBUTION OF THIS PROGRAM CONSTITUTES RECIPIENT'S
 * ACCEPTANCE OF THE OSMC PUBLIC LICENSE.
 *
 * The OpenModelica software and the Open Source Modelica
 * Consortium (OSMC) Public License (OSMC-PL) are obtained
 * from Linköping University, either from the above address,
 * from the URLs: http://www.ida.liu.se/projects/OpenModelica or
 * http://www.openmodelica.org, and in the OpenModelica distribution.
 * GNU version 3 is obtained from: http://www.gnu.org/copyleft/gpl.html.
 *
 * This program is distributed WITHOUT ANY WARRANTY; without
 * even the implied warranty of  MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE, EXCEPT AS EXPRESSLY SET FORTH
 * IN THE BY RECIPIENT SELECTED SUBSIDIARY LICENSE CONDITIONS
 * OF OSMC-PL.
 *
 * See the full OSMC Public License conditions for more details.
 *
 */

/*
 * Developed by:
 * FH-Bielefeld
 * Developer: Vitalij Ruge
 * Contact: vitalij.ruge@fh-bielefeld.de
 */

#include "../ipoptODEstruct.h"
#include "../localFunction.h"

#ifdef WITH_IPOPT

/* static  int jac_struc(Index *iRow, Index *iCol, long int nx, long int nv, int nsi); */
static  int radauJac1(double *a, double *J, double dt, double scalRes, double * values, int nv, int *k, int j,IPOPT_DATA_ *iData);
static  int lobattoJac1(double *a, double *J, double *J0, double dt, double scalRes, double * values, int nv, int *k, int j, long double tmp,IPOPT_DATA_ *iData);
static  int radauJac2(double *a, double *J, double dt, double scalRes, double * values, int nv, int *k, int j,IPOPT_DATA_ *iData);
static  int lobattoJac2(double *a, double *J, double *J0, double dt, double scalRes, double * values, int nv, int *k, int j, long double tmp,IPOPT_DATA_ *iData);
static  int radauJac3(double *a, double *J, double dt, double scalRes, double * values, int nv, int *k, int j,IPOPT_DATA_ *iData);
static  int lobattoJac3(double *a, double *J, double *J0, double dt, double scalRes, double * values, int nv, int *k, int j, long double tmp,IPOPT_DATA_ *iData);
static int jac_struc(IPOPT_DATA_ *iData,int *iRow, int *iCol);
static int conJac(double *J, double * values, int nv, int *k, int j,IPOPT_DATA_ *iData);


/*!
 *  eval derivation of s.t.
 *  author: Vitalij Ruge
 **/
Bool evalfDiffG(Index n, double * x, Bool new_x, Index m, Index njac, Index *iRow, Index *iCol, Number *values, void * useData)
{
  int i,j,k,l;
  IPOPT_DATA_ *iData;
  k = 0;

  iData = (IPOPT_DATA_ *) useData;
  if(values == NULL)
  {
    /* jac_struc(iRow, iCol, iData->nx, iData->nv, iData->nsi); */
    jac_struc(iData, iRow, iCol);

   /*
    printf("\n m = %i , %i",m ,iData->NRes);
    printf("\nk = %i , %i" ,k ,njac);
    for(i = 0; i< njac; ++i)
      printf("\nJ(%i,%i) = 1; i= %i",iRow[i]+1, iCol[i]+1,i);

    assert(0);
    */

  }
  else
  {
    ipoptDebuge(iData,x);
    long double tmp[3];
    int id;
    int nng = iData->nx + iData->nc;
    double *a1, *a2, *a3;
    double *d1, *d2, *d3;

    a1 = iData->a1;
    a2 = iData->a2;
    a3 = iData->a3;

    d1 = iData->d1;
    d2 = iData->d2;
    d3 = iData->d3;

    tmp[0] = iData->dt[0]*iData->d1[4];
    tmp[1] = iData->dt[0]*iData->d2[4];
    tmp[2] = iData->dt[0]*iData->d3[4];

    for(i = 0, id = iData->nv; i<1; ++i)
    {
      diff_functionODE0(x, iData->t0 , iData);
      for(l=0; l<iData->deg; ++l, id += iData->nv)
      {
        diff_functionODE(x+id , iData->time[i*iData->deg + l] , iData, iData->J);
        for(j=0; j<iData->nx; ++j)
        {
          switch(l)
          {
          case 0:
            lobattoJac1(d1, iData->J[j], iData->J0[j], iData->dt[i], 1.0, values, iData->nv, &k, j, tmp[0],iData);
            break;
          case 1:
            lobattoJac2(d2, iData->J[j], iData->J0[j], iData->dt[i], 1.0, values, iData->nv, &k, j, tmp[1],iData);
            break;
          case 2:
            lobattoJac3(d3, iData->J[j], iData->J0[j], iData->dt[i], 1.0, values, iData->nv, &k, j, tmp[2],iData);
            break;
          }
        }
        for(;j<nng; ++j)
        {
          conJac(iData->J[j], values, iData->nv, &k, j,iData);
        }
      }


    }

    for(; i<iData->nsi; ++i)
    {
      for(l=0; l<iData->deg; ++l, id += iData->nv)
      {
        diff_functionODE(x+id, iData->time[i*iData->deg + l], iData, iData->J);
        for(j=0; j<iData->nx; ++j)
        {
          switch(l)
          {
          case 0:
            radauJac1(a1, iData->J[j], iData->dt[i], 1.0, values, iData->nv, &k, j,iData);
            break;
          case 1:
            radauJac2(a2, iData->J[j], iData->dt[i], 1.0, values, iData->nv, &k, j,iData);
            break;
          case 2:
            radauJac3(a3, iData->J[j], iData->dt[i], 1.0, values, iData->nv, &k, j,iData);
            break;
          }
        }
        for(;j<nng; ++j)
          conJac(iData->J[j], values, iData->nv, &k, j,iData);
      }
    }
     /*assert(k == njac);*/
  }
  return TRUE;
}

/* static  int jac_struc(Index *iRow, Index *iCol, long int nx, long int nv, int nsi) */

/*!
 *  special jacobian struct
 *  author: Vitalij Ruge
 **/
static  int radauJac1(double *a, double *J, double dt, double scalRes, double * values, int nv, int *k, int j,IPOPT_DATA_ *iData)
{
  int l;
  values[(*k)++] = a[0];
  values[(*k)-1] *= scalRes;
  /*1*/
  for(l=0; l<nv; ++l)
  {
    if(iData->knowedJ[j][l] == 1)
    {
      values[(*k)++] = (j == l) ? dt*J[l]-a[1] : dt*J[l];
      values[(*k)-1] *= scalRes;
    }
  }

  /*2*/
  values[(*k)++] = -a[2];
  values[(*k)-1] *= scalRes;

  /*3*/
  values[(*k)++] = a[3];
  values[(*k)-1] *= scalRes;
  return 0;
}

/*!
 *  special jacobian struct
 *  author: Vitalij Ruge
 **/
static  int lobattoJac1(double *a, double *J, double *J0, double dt, double scalRes, double * values, int nv, int *k, int j, long double tmp,IPOPT_DATA_ *iData)
{
  int l;
  /*0*/
  for(l = 0; l< nv; ++l)
  {
    if(j == l) {
      values[(*k)++] = tmp*J0[l] + a[0];
      values[(*k)-1] *= scalRes;
    } else if(iData->knowedJ[j][l] == 1) {
      values[(*k)++] = tmp*J0[l];
      values[(*k)-1] *= scalRes;
    }
  }
  /*1*/
  for(l = 0; l< nv; ++l)
  {
    if(iData->knowedJ[j][l] == 1)
    {
      values[(*k)++] = ((j == l)? dt*J[l] - a[1] : dt*J[l]);
      values[(*k)-1] *= scalRes;
    }
  }
  /*2*/
  values[(*k)++] = -a[2];
  values[(*k)-1] *= scalRes;

  /*3*/
  values[(*k)++] = a[3];
  values[(*k)-1] *= scalRes;
  return 0;
}


/*!
 *  special jacobian struct
 *  author: Vitalij Ruge
 **/
static  int radauJac2(double *a, double *J, double dt, double scalRes, double * values, int nv, int *k, int j,IPOPT_DATA_ *iData)
{
  int l;
  /*0*/
  values[(*k)++] = -a[0];
  values[(*k)-1] *= scalRes;

  /*1*/
  values[(*k)++] = a[1];
  values[(*k)-1] *= scalRes;

  /*2*/
  for(l = 0; l< nv; ++l)
  {
    if(iData->knowedJ[j][l] == 1)
    {
      values[(*k)++] = ((j == l)? dt*J[l] - a[2] : dt*J[l]);
      values[(*k)-1] *= scalRes;
    }
  }

  /*3*/
  values[(*k)++] = -a[3];
  values[(*k)-1] *= scalRes;
  return 0;
}

/*!
 *  special jacobian struct
 *  author: Vitalij Ruge
 **/
static  int lobattoJac2(double *a, double *J, double *J0, double dt, double scalRes, double * values, int nv, int *k, int j, long double tmp,IPOPT_DATA_ *iData)
{
  int l;
  /*0*/
  for(l = 0; l< nv; ++l)
  {
    if( j==l){
      values[(*k)++] = -(tmp*J0[l] + a[0]);
      values[(*k)-1] *= scalRes;
    } else if(iData->knowedJ[j][l] == 1) {
      values[(*k)++] = -tmp*J0[l];
      values[(*k)-1] *= scalRes;
    }
  }
  /*1*/
  values[(*k)++] = a[1];
  values[(*k)-1] *= scalRes;

  /*2*/
  for(l = 0; l< nv; ++l)
  {
    if(iData->knowedJ[j][l] == 1)
    {
      values[(*k)++] = ((j == l)? dt*J[l]-a[2] : dt*J[l]);
      values[(*k)-1] *= scalRes;
    }
  }
  /*3*/
  values[(*k)++] = -a[3];
  values[(*k)-1] *= scalRes;
  return 0;
}

/*!
 *  special jacobian struct
 *  author: Vitalij Ruge
 **/
static  int radauJac3(double *a, double *J, double dt, double scalRes, double * values, int nv, int *k, int j,IPOPT_DATA_ *iData)
{
  int l;
  /*0*/
  values[(*k)++] = a[0];
  values[(*k)-1] *= scalRes;
  /*1*/
  values[(*k)++] = -a[1];
  values[(*k)-1] *= scalRes;
  /*2*/
  values[(*k)++] = a[2];
  values[(*k)-1] *= scalRes;
  /*3*/
  for(l = 0; l< nv; ++l)
  {
    if(iData->knowedJ[j][l] == 1)
    {
      values[(*k)++] = ((j == l)? dt*J[l] - a[3] : dt*J[l]);
      values[(*k)-1] *= scalRes;
    }
  }
  return 0;
}

/*!
 *  special jacobian struct
 *  author: Vitalij Ruge
 **/
static  int lobattoJac3(double *a, double *J, double *J0, double dt, double scalRes, double * values, int nv, int *k, int j, long double tmp,IPOPT_DATA_ *iData)
{
  int l;
  /*0*/
  for(l=0; l<nv; ++l)
  {
    if(j==l){
      values[(*k)++] = tmp*J0[l] + a[0];
      values[(*k)-1] *= scalRes;
    } else if(iData->knowedJ[j][l] == 1) {
      values[(*k)++] = tmp*J0[l];
      values[(*k)-1] *= scalRes;
    }
  }
  /*1*/
  values[(*k)++] = -a[1];
  values[(*k)-1] *= scalRes;
  /*2*/
  values[(*k)++] = a[2];
  values[(*k)-1] *= scalRes;
  /*3*/
  for(l=0; l<nv; ++l)
  {
    if(iData->knowedJ[j][l] == 1)
    {
      values[(*k)++] = ((j == l)? dt*J[l] - a[3] : dt*J[l]);
      values[(*k)-1] *= scalRes;
    }
  }
  return 0;
}

static int conJac(double *J, double *values, int nv, int *k, int j,IPOPT_DATA_ *iData)
{
  int l;
  for(l=0; l<nv; ++l)
    if(iData->knowedJ[j][l] == 1)
      values[(*k)++] = J[l];
  return 0;
}

/*!
 *  special jacobian struct
 *  author: Vitalij Ruge
 **/
static int jac_struc(IPOPT_DATA_ *iData, int *iRow, int *iCol)
{
  int nr, nc,r,c, nv,nsi,nx, nJ;
  int i,j,k=0,l;

  nv = iData->nv;
  nx = iData->nx;
  nsi = iData->nsi;
  nJ = iData->nc + nx;
  r = 0;
  c = nv;
  for(i=0; i<nsi; ++i)
  {
    /*1*/
    for(j=0; j<nx; ++j)
    {
      /*0*/
      if(i > 0)
      {
         iRow[k] = r + j;
         iCol[k++] = c - nv + j;
      }
      else
      {
        for(l=0; l<nv; ++l)
        {
          if(iData->knowedJ[j][l] == 1)
          {
            iRow[k] = j;
            iCol[k++] = l;
          }
        }
      }

      /*1*/
      for(l=0; l<nv; ++l)
      {
        if(iData->knowedJ[j][l] == 1)
        {
          iRow[k] = r + j;
          iCol[k++] = c + l;
        }
      }

      /*2*/
      iRow[k] = r + j;
      iCol[k++] = c + nv + j;
      /*3*/
      iRow[k] = r + j;
      iCol[k++] = c + 2*nv + j;
    }
    for(;j<nJ; ++j)
    {
        /*1*/
        for(l=0; l<nv; ++l)
        {
          if(iData->knowedJ[j][l] == 1)
          {
            iRow[k] = r + j;
            iCol[k++] = c + l;
          }
        }
    }

    /*2*/
    r += nJ;
    c += nv;

    for(j=0; j<nx; ++j)
    {
      /*0*/
      if(i > 0)
      {
        iRow[k] = r + j;
        iCol[k++] = c - 2*nv + j;
      }
      else
      {
        for(l=0; l<nv; ++l)
        {
          if(iData->knowedJ[j][l] == 1)
          {
            iRow[k] = r + j;
            iCol[k++] = l;
          }
        }
      }

      /*1*/
      iRow[k] = iRow[k-1];
      iCol[k++] =  c - nv + j;

      /*2*/
      for(l=0; l<nv; ++l)
      {
        if(iData->knowedJ[j][l] == 1)
        {
          iRow[k] = iRow[k-1];
          iCol[k++] = c + l;
        }
      }

      /*3*/
       iRow[k] = iRow[k-1];
       iCol[k++] = c + nv + j;
    }
    for(;j<nJ; ++j)
    {
        /*2*/
        for(l=0; l<nv; ++l)
        {
          if(iData->knowedJ[j][l] == 1)
          {
            iRow[k] = r + j;
            iCol[k++] = c + l;
          }
        }
    }

    /*3*/
    r += nJ;
    c += nv;

    for(j=0; j<nx; ++j)
    {
      /*0*/
      if(i > 0)
      {
        iRow[k] = r + j;
        iCol[k++] = c - 3*nv + j;
      }
      else
      {
        for(l=0; l<nv; ++l)
        {
          if(iData->knowedJ[j][l] == 1)
          {
            iRow[k] = r + j;
            iCol[k++] = l;
          }
        }
      }

      /*1*/
      iRow[k] = iRow[k-1];
      iCol[k++] = c - 2*nv + j;

      /*2*/
      iRow[k] = iRow[k-1];
      iCol[k++] = c - nv + j;

      /*3*/
      for(l=0; l<nv; ++l)
      {
        if(iData->knowedJ[j][l] == 1)
        {
          iRow[k] = r + j;
          iCol[k++] = c + l;
        }
      }
    }
    for(;j<nJ; ++j)
    {
        /*3*/
        for(l=0; l<nv; ++l)
        {
          if(iData->knowedJ[j][l] == 1)
          {
            iRow[k] = r + j;
            iCol[k++] = c + l;
          }
        }
    }

    r += nJ;
    c += nv;
  }

    /*
    printf("\n\n%i = %i",iData->njac,k);
    assert(0);
    */


  return 0;
}

#endif
