/*********************************************************************
Functions for parametric fitting.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
     Giacomo Lorenzetti <glorenzetti@cefca.es>
Copyright (C) 2022-2026 Free Software Foundation, Inc.

Gnuastro is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

Gnuastro is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with Gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <config.h>

#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <string.h>
#include <stdlib.h>

#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>

#include <gnuastro/fit.h>
#include <gnuastro/blank.h>
#include <gnuastro/pointer.h>

#include <gnuastro-internal/checkset.h>





/**********************************************************************/
/****************              Identifiers             ****************/
/**********************************************************************/
int
gal_fit_name_robust_to_id(char *name)
{
  /* In case 'name' is NULL, then return the invalid type. */
  if(name==NULL) return GAL_FIT_ROBUST_INVALID;

  /* Match the name. */
  if(      !strcmp(name, "bisquare") ) return GAL_FIT_ROBUST_BISQUARE;
  else if( !strcmp(name, "cauchy")   ) return GAL_FIT_ROBUST_CAUCHY;
  else if( !strcmp(name, "fair")     ) return GAL_FIT_ROBUST_FAIR;
  else if( !strcmp(name, "huber")    ) return GAL_FIT_ROBUST_HUBER;
  else if( !strcmp(name, "ols")      ) return GAL_FIT_ROBUST_OLS;
  else if( !strcmp(name, "welsch")   ) return GAL_FIT_ROBUST_WELSCH;
  else                                 return GAL_FIT_ROBUST_INVALID;

  /* If control reaches here, there was a bug! */
  error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at '%s' to "
        "find a fix it. Control should not have reached this point",
        __func__, PACKAGE_BUGREPORT);
  return GAL_FIT_ROBUST_INVALID;
}





char *
gal_fit_name_robust_from_id(uint8_t robustid)
{
  switch(robustid)
    {
    case GAL_FIT_ROBUST_BISQUARE:   return "bisquare";
    case GAL_FIT_ROBUST_CAUCHY:     return "cauchy";
    case GAL_FIT_ROBUST_FAIR:       return "fair";
    case GAL_FIT_ROBUST_HUBER:      return "huber";
    case GAL_FIT_ROBUST_OLS:        return "ols";
    case GAL_FIT_ROBUST_WELSCH:     return "welsch";
    default:                        return NULL;
    }

  /* If control reaches here, there was a bug! */
  error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at '%s' to "
        "find a fix it. Control should not have reached this point",
        __func__, PACKAGE_BUGREPORT);
  return NULL;
}





static char *
fit_name_matrix_from_id(uint8_t matrixid)
{
  switch(matrixid)
    {
    case GAL_FIT_MATRIX_POLYNOMIAL_1D: return "polynomial-1d";
    case GAL_FIT_MATRIX_POLYNOMIAL_2D: return "polynomial-2d";
    case GAL_FIT_MATRIX_POLYNOMIAL_2D_TPV: return "polynomial-2d-tpv";
    case GAL_FIT_MATRIX_POLYNOMIAL_2D_TPV_NO_RADIAL:
      return "polynomial-2d-tpv-no-radial";
    default: return NULL;
    }

  /* If control reaches here, there was a bug! */
  error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at '%s' to "
        "find a fix it. Control should not have reached this point",
        __func__, PACKAGE_BUGREPORT);
  return NULL;
}




















/**********************************************************************/
/****************            Common to all             ****************/
/**********************************************************************/
static gal_data_t *
fit_sanity_check_col(gal_data_t *in, gal_data_t *ref, const char *func)
{
  gal_data_t *out;

  /* Make sure the input is 1-dimensional. */
  if(in->ndim!=1)
    error(EXIT_FAILURE, 0, "%s: inputs must have one dimension", func);

  /* Make sure the input has the same size as the reference. */
  if(in->size != ref->size)
    error(EXIT_FAILURE, 0, "%s: all inputs must have the same size",
          func);

  /* Make sure output has a double type. */
  out = ( in->type==GAL_TYPE_FLOAT64
          ? in
          : gal_data_copy_to_new_type(in, GAL_TYPE_FLOAT64) );

  /* If there are blank values, print a warning, then return. */
  if(gal_blank_present(out, 1))
    error(EXIT_SUCCESS, 0, "%s: at least one of the input columns "
          "have a blank value; the fit will become NaN. Within the "
          "Gnuastro, you can use 'gal_blank_remove_rows' to remove "
          "all rows that have at least one blank value in any column",
          func);
  return out;
}




















/**********************************************************************/
/****************              Linear fit              ****************/
/**********************************************************************/
static gal_data_t *
fit_linear_1d_base(gal_data_t *xin, gal_data_t *yin, gal_data_t *ywht,
                   double *redchisq, int withconstant)
{
  size_t osize;
  double *o, sqs, nparam=NAN;
  gal_data_t *x=NULL, *y=NULL, *w=NULL, *out;

  /* Basic sanity checks. */
  if(xin==NULL || xin->size==0 || yin==NULL || yin->size==0)
    error(EXIT_FAILURE, 0, "%s: the inputs are either NULL or "
          "do not contain any data (have a length of zero)", __func__);
  x=fit_sanity_check_col(xin, xin, __func__);
  y=fit_sanity_check_col(yin, xin, __func__);
  if(ywht) w=fit_sanity_check_col(ywht, xin, __func__);

  /* Allocate the output dataset. */
  osize = withconstant ? 5 : 2;
  out=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &osize, NULL, 0,
                     -1, 1, NULL, NULL, NULL);

  /* For a check.
  {
    size_t i;
    double *xa=x->array, *ya=y->array;
    for(i=0;i<x->size;++i)
      printf("%-15f %-15f\n", xa[i], ya[i]);
  } //*/

  /* Do the fitting. */
  o=out->array;
  if(withconstant)
    {
      nparam=2;
      if(ywht)
        gsl_fit_wlinear(x->array, 1, w->array, 1, y->array, 1, x->size,
                        o, o+1, o+2, o+3, o+4, &sqs);
      else
        gsl_fit_linear(x->array, 1, y->array, 1, x->size, o, o+1,
                       o+2, o+3, o+4, &sqs);
    }
  else
    {
      nparam=1;
      if(ywht)
        gsl_fit_wmul(x->array, 1, w->array, 1, y->array, 1, x->size,
                     o, o+1, &sqs);
      else
        gsl_fit_mul(x->array, 1, y->array, 1, x->size, o, o+1, &sqs);
    }

  /* For a check.
  {
    printf("c0: %f\nc1: %f\n"
           "cov00: %f\ncov01: %f\ncov11: %f\nsumsq: %f\n",
           o[0], o[1], o[2], o[3], o[4], o[5]);
  } //*/

  /* Calculate the reduced chi^2: As mentioned in [1], in case we have the
     chi^2, then it is simply the chi^2 divided by the degrees of
     freedom. GSL returns the chi^2 for weighted fits and the sum of
     squares for non-weighted fits [2]. This is because without weights,
     the chi^2 is the same as the sum of squares [1].

     The number of degrees of freedom is defined by the number of
     observations subtracted from the number of fitted parameters.

     [1] https://en.wikipedia.org/wiki/Reduced_chi-squared_statistic
     [2] https://www.gnu.org/software/gsl/doc/html/lls.html */
  *redchisq = sqs / (x->size - nparam);

  /* Clean up and return. */
  if(x!=xin) gal_data_free(x);
  if(y!=yin) gal_data_free(y);
  if(ywht && w!=ywht) gal_data_free(w);
  return out;
}





gal_data_t *
gal_fit_linear_1d(gal_data_t *xin, gal_data_t *yin, gal_data_t *ywht,
                  double *redchisq)
{
  return fit_linear_1d_base(xin, yin, ywht, redchisq, 1);
}





gal_data_t *
gal_fit_linear_no_constant_1d(gal_data_t *xin, gal_data_t *yin,
                              gal_data_t *ywht, double *redchisq)
{
  return fit_linear_1d_base(xin, yin, ywht, redchisq, 0);
}





static gal_data_t *
fit_estimate_prepare(gal_data_t *xin, gal_data_t *fit, gal_data_t **xd,
                     const char *func)
{
  gal_data_t *ftmp, *out=NULL;

  /* The Fit arrays should be double precision. 1D fits produce a single
     'gal_data_t' and 2D arrays produce a list of two 'gal_data_t's. */
  for(ftmp=fit; ftmp!=NULL; ftmp=ftmp->next)
    if(ftmp->type!=GAL_TYPE_FLOAT64
       || (ftmp->next && fit->next->type!=GAL_TYPE_FLOAT64) )
      error(EXIT_FAILURE, 0, "%s: the 'fit' argument should only "
            "contain double precision floating point types", func);
  if(fit->ndim!=1 || (fit->next && fit->next->ndim!=2) )
    error(EXIT_FAILURE, 0, "%s: the 'fit' argument should only "
          "contain single-dimensional data", func);
  if(fit->next && (fit->next->dsize[0]!=fit->next->dsize[1]))
    error(EXIT_FAILURE, 0, "%s: the secont dataset of the 'fit' "
          "argument should be square (same size in both "
          "dimensions)", func);
  if(xin->ndim>2)
    error(EXIT_FAILURE, 0, "%s: currenly only 1D and 2D datasets "
          "are supported, but the input has %zu dimensions", __func__,
          xin->ndim);
  if(xin->ndim==2 && xin->next)
    error(EXIT_FAILURE, 0, "%s: the 'xin' input has %zu dimentions, "
          "as well as a 'next' element: this is not an expected "
          "situation here. A 2D input can be a 2D 'xin', but without "
          "a 'next' element, or a list of two 1D datasets connected by "
          "'next'", __func__, xin->ndim);
  if(xin->ndim==2 && xin->array)
    error(EXIT_FAILURE, 0, "%s: the 'xin' argument has a non-NULL "
          "'array' and has %zu dimensions. This argument should "
          "either be a list of 1D arrays/columns, or a single 2D "
          "dataset with 'array==NULL'", __func__, xin->ndim);
  if(xin->next)
    {
      if(xin->ndim!=1)
        error(EXIT_FAILURE, 0, "%s: when the 'xin' argument is a list, "
              "each node should only have a single dimension, but the "
              "first node has %zu dimensions", __func__, xin->ndim);
      if(gal_dimension_is_different(xin, xin->next))
        error(EXIT_FAILURE, 0, "%s: when the 'xin' argument is a list, "
              "all nodes should only have the same size, but that is not "
              "the case", __func__);
      if(xin->next->next)
        error(EXIT_FAILURE, 0, "%s: when 'xin' is a list of 1D "
              "datasets, there should only be two nodes in the list",
              __func__);
    }

  /* Make sure the input X values are in double precision. It can happen
     that 'xin->array==NULL': in which case, we just need the size of the
     array later and its type is irrelevant. */
  *xd = ( xin->array
          ? ( xin->type==GAL_TYPE_FLOAT64
              ? xin
              : gal_data_copy_to_new_type(xin, GAL_TYPE_FLOAT64) )
          : xin );
  if(xin->next)
    (*xd)->next = ( xin->next->array
                    ? ( xin->next->type==GAL_TYPE_FLOAT64
                        ? xin->next
                        : gal_data_copy_to_new_type(xin->next,
                                                    GAL_TYPE_FLOAT64) )
                    : xin->next );

  /* Allocate the output datasets. */
  gal_list_data_add_alloc(&out, NULL, GAL_TYPE_FLOAT64, xin->ndim,
                          xin->dsize, xin->wcs, 1, xin->minmapsize,
                          xin->quietmmap, NULL, NULL, NULL);
  gal_list_data_add_alloc(&out, NULL, GAL_TYPE_FLOAT64, xin->ndim,
                          xin->dsize, xin->wcs, 1, xin->minmapsize,
                          xin->quietmmap, NULL, NULL, NULL);
  gal_list_data_reverse(&out);

  /* Return the output. */
  return out;
}





gal_data_t *
gal_fit_linear_estimate_1d(gal_data_t *fit, gal_data_t *xin)
{
  size_t i;
  gal_data_t *out=NULL, *xd;
  double *x, *y, *yerr, *f=fit->array;

  /* Do the basic preparations. */
  out=fit_estimate_prepare(xin, fit, &xd, __func__);

  /* Set the pointers. */
  x    = xd->array;
  y    = out->array;
  yerr = out->next->array;

  /* Estimate the values. */
  switch(fit->size)
    {
    case 6:                     /* Linear with constant. */
      for(i=0;i<out->size;++i)
        gsl_fit_linear_est(x[i], f[0], f[1], f[2], f[3], f[4],
                           y+i, yerr+i);
      break;

    case 3:                     /* Linear WITHOUT constant. */
      for(i=0;i<out->size;++i)
        gsl_fit_mul_est(x[i], f[0], f[1], y+i, yerr+i);
      break;

    default:                    /* Un-recognized situation! */
      error(EXIT_FAILURE, 0, "%s: the 'fit' argument should "
            "either have 6 or 3 elements (be an output of "
            "'gal_fit_linear_1d' or 'gal_fit_linear_1d_no_constant'"
            "respectively), but it has %zu elements", __func__,
            fit->size);
    }

  /* Clean up. */
  if(xd!=xin) gal_data_free(xd);
  return out;
}




















/**********************************************************************/
/****************           Polynomial fits            ****************/
/**********************************************************************/
static size_t
fit_polynomial_nconst(uint8_t degree, uint8_t matrixid)
{
  size_t nconst=GAL_BLANK_SIZE_T;

  switch(matrixid)
    {
    case GAL_FIT_MATRIX_POLYNOMIAL_1D:
      nconst=degree+1;
      break;
    case GAL_FIT_MATRIX_POLYNOMIAL_2D:
      nconst=(degree+1)*(degree+2)/2;
      break;
    case GAL_FIT_MATRIX_POLYNOMIAL_2D_TPV:
      /* TPV has odd radial terms. The radial terms involve a mixture of
         the two dimensions, so when a robust fit is requested, solving
         using a big 'block diagonal' matrix is better than 2 small
         independent matrices. */
      nconst=(degree+1)*(degree+2)/2 + (degree/2 + 1);
      nconst*=2;
      break;
    case GAL_FIT_MATRIX_POLYNOMIAL_2D_TPV_NO_RADIAL:
      /* TPV without radial terms is equivalent to a 2D polynomial in a
         block diagonal matrix. */
      nconst=(degree+1)*(degree+2)/2;
      nconst*=2;
      break;
    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at "
            "'%s' to fix the problem. The value '%u' is not "
            "recognized for 'matrixid'", __func__, PACKAGE_BUGREPORT,
            matrixid);
    }

  return nconst;
}





static void
fit_polynomial_row_1d(double *xi, size_t i, size_t nconst, double *xo)
{
  size_t j;

  /* The first column (constant) of this row doesn't depend on X. So we'll
     give it a value of 1.0. */
  xo[0] = 1.0f;

  /* Column j is the multiplication of column j-1 with the input horizontal
     value. This will make it a polynomial. */
  for(j=1;j<nconst;++j)
    xo[ j ] = xo[ j-1 ] * xi[ i ];
}





static void
fit_polynomial_matrix_fill_1d(gal_data_t *xin, int nconst, double *xo)
{
  size_t i;
  double *xi;

  /* Fill in the X matrix. */
  xi=xin->array;
  for(i=0;i<xin->size;++i)
    fit_polynomial_row_1d(xi, i, nconst, &xo[i*nconst]);

  /* For a check.
  {
    size_t checki=5;
    printf("Row %zu: ", checki);
    for(j=0;j<degree;++j)
      printf("%.3f ", xo[ checki*degree + j ]);
    printf("\n");
    exit(0);
  } //*/
}





/* Compute the powers of the elements that are then combined to build the
   polynomial. */
static void
fit_polynomial_row_init_2d(double xi1, double xi2, size_t degree,
                           double *xi1_pow, double *xi2_pow,
                           double *r_pow)
{
  size_t deg;
  double rsq = r_pow ? (xi1*xi1 + xi2*xi2) : NAN;

  /* Fill the first radial term. */
  if(r_pow) r_pow[1]=sqrt(rsq);

  /* Go over all the degrees. */
  for(deg=0; deg<=degree; deg++)
    {
      xi1_pow[deg] = deg>0 ? xi1*xi1_pow[deg-1] : 1.0;
      xi2_pow[deg] = deg>0 ? xi2*xi2_pow[deg-1] : 1.0;

      /* If a TPV polynomial with radial terms is requested, add odd powers
         of the radial term. See
         https://fits.gsfc.nasa.gov/registry/tpvwcs/tpv.html */
      if(r_pow && deg>1 && deg%2) r_pow[deg] = r_pow[deg-2] * rsq;
    }
}





/* Compute and store the values of a row. The first argument is the
   starting position: this is useful for instance when the matrix is block
   diagonal, and hence the rows starts with several zeroes */
static void
fit_polynomial_row_2d(double *xi1_pow, double *xi2_pow, double *r_pow,
                      size_t degree, double *row)
{
  size_t j, k=0, deg;

  /* Column k is a combination of powers of the input values.  This will
     make it a polynomial. */
  for(deg=0; deg<=degree; deg++)
    {
      /* Use the previously computed powers to fill the matrix */
      for(j=deg+1; j-->0;)
        {
          /* For a check on the powers of the dimensions.
          if(i==0) printf("%s: %zu, %zu\n", __func__, j, deg-j);
          //*/

          row[k++] = xi1_pow[j] * xi2_pow[deg-j];
        }

      /* If a tpv polynomial is requested, add odd powers of the radial
         term. See https://fits.gsfc.nasa.gov/registry/tpvwcs/tpv.html */
      if(r_pow && deg%2) row[k++] = r_pow[deg];
    }
}





static void
fit_polynomial_matrix_fill_2d(gal_data_t *xin, int nconst,
                              double *xo, size_t degree,
                              uint8_t tpv, uint8_t tpvradial)
{
  size_t i;
  gal_data_t *xin1=xin, *xin2=xin->next;
  double *rp, *row, *xi1, *xi2, *xi1p, *xi2p;

  /* Allocate the intermediate arrays that keep the powers (hence the "p"
     suffix) of each coordiante and the radius. */
  xi1p=gal_pointer_allocate(GAL_TYPE_FLOAT64, degree+1, 1, __func__,
                            "xi1p");
  xi2p=gal_pointer_allocate(GAL_TYPE_FLOAT64, degree+1, 1, __func__,
                            "xi2p");
  rp = ( tpvradial
         ? gal_pointer_allocate(GAL_TYPE_FLOAT64, degree+1, 1,
                                __func__, "rp")
         : NULL );

  /* Initialize the array pointers and fill them. */
  xi1=xin1->array;
  xi2=xin2->array;
  for(i=0;i<xin1->size;i++)
    {
      /* Initialize and fill the matrix. */
      row=xo+i*nconst;
      fit_polynomial_row_init_2d(xi1[i], xi2[i], degree, xi1p, xi2p, rp);
      fit_polynomial_row_2d(xi1p, xi2p, rp, degree, row);

      /* For the TPV matrix, we actually need to fit two polynomials at the
         same time, so the matrix is double the size in each dimension:
         four times larger. The 'fit_polynomial_row_2d' function above only
         fills the first quarter so below we need to fill the last one (the
         other two quarters are empty). Furthermore, since the second
         polynomial is for the second coordinate, the 'xi1p' and 'xi2p' get
         reversed in the new call. */
      if(tpv)
        {
          /* Put the second coordinate in the second half of the matrix */
          row = xo + (xin1->size+i)*nconst + nconst/2;

          /* Fill the final quarter of the matrix. Note that xi1 and xi2
             are flipped because the last filled quarter of the matrix is
             for the second dimension in the TPV fit. */
          fit_polynomial_row_2d(xi2p, xi1p, rp, degree, row);
        }
    }

  /* For a check.
  {
    size_t j, checki=2;
    printf("Row %zu: ", checki);
    for(j=0;j<nconst;++j)
      printf("%.3f ", xo[ checki*nconst + j ]);
    printf("\n");
    exit(0);
  } //*/

  /* Clean up. */
  free(xi1p);
  free(xi2p);
  if(rp) free(rp);
}





static void
fit_polynomial_prepare(gal_data_t *xin,  gal_data_t *yin,
                       gal_data_t *ywht, int nconst,
                       gsl_matrix **x,   gsl_vector **c,
                       gsl_matrix **cov, gsl_vector *y,
                       gsl_vector *w,    size_t degree,
                       uint8_t matrixid)
{
  double *xo;

  /* Use GSL's own matrix allocation functions for the structures that need
     allocation and we can't use the same allocated space of the inputs. */
  *c   = gsl_vector_alloc(nconst);
  *cov = gsl_matrix_alloc(nconst, nconst);
  *x   = gsl_matrix_calloc(yin->size, nconst);

  /* Fill the design matrix. */
  xo=(*x)->data;
  switch(matrixid)
    {
    case GAL_FIT_MATRIX_POLYNOMIAL_1D:
      fit_polynomial_matrix_fill_1d(xin, nconst, xo);
      break;
    case GAL_FIT_MATRIX_POLYNOMIAL_2D:
      fit_polynomial_matrix_fill_2d(xin, nconst, xo, degree, 0, 0);
      break;
    case GAL_FIT_MATRIX_POLYNOMIAL_2D_TPV:
      fit_polynomial_matrix_fill_2d(xin, nconst, xo, degree, 1, 1);
      break;
    case GAL_FIT_MATRIX_POLYNOMIAL_2D_TPV_NO_RADIAL:
      fit_polynomial_matrix_fill_2d(xin, nconst, xo, degree, 1, 0);
      break;
    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at '%s' to "
            "fix the problem. '%d' isn't recognized as a design matrix ID",
            __func__, PACKAGE_BUGREPORT, matrixid);
    }

  /* Set the pointers of the 'y' and 'w' GSL vectors. */
  y->data=yin->array;
  if(ywht) w->data=ywht->array;
}





static void
fit_polynomial_sanity_check(gal_data_t *xin, gal_data_t *yin,
                            gal_data_t *ywht, uint8_t matrixid,
                            gal_data_t **xdata, gal_data_t **ydata,
                            gal_data_t **wdata, size_t degree,
                            uint8_t robustid, double tikhonovlambda)
{
  /* Check the dimensionality. */
  if(xin->next)
    {
      if(xin->next->next)
        error(EXIT_FAILURE, 0, "%s: only 1d and 2d polynomials are "
              "currently supported, but more than 2 columns have "
              "been given as the first independenat variable ('x'). "
              "This might happen if the fitting columns directly "
              "came from a reading of a table and the columns were "
              "not separated ('xin->next' or xin->next->next were "
              "are not set to NULL before calling 'gal_fit_polynomial')",
              __func__);
      else if(xin->next==yin)
        error(EXIT_FAILURE, 0, "%s: the second independent variable "
              "points to the measurement variable. Set 'x->next=NULL' "
              "before calling 'gal_fit_polynomial'", __func__);
      else if(matrixid < GAL_FIT_MATRIX_NUMBER_1D)
        error(EXIT_FAILURE, 0, "%s: 'matrixid=%s' describes a 1d design "
              "matrix, but 'xin' is a list of more than one datasets",
              __func__, fit_name_matrix_from_id(matrixid));
    }
  else if(matrixid > GAL_FIT_MATRIX_NUMBER_1D)
    error(EXIT_FAILURE, 0, "%s: 'matrixid=%s' describes a 2d design "
          "matrix, but 'xin' is a list of just one dataset",
          __func__, fit_name_matrix_from_id(matrixid));

  /* Tikhonov regularized regression is not yet implemented with weights,
     https://www.gnu.org/s/gsl/doc/html/lls.html#regularized-regression */
  if(!isnan(tikhonovlambda) && ywht)
    error(EXIT_FAILURE, 0, "%s: tikhonov regularized fitting with "
          "weights is still not implemented. Please use "
          "'gal_fit_polynomial()' or set 'ywht=NULL'", __func__);

  /* Regularized and robust fits cannot be called simultaneously. */
  if(!isnan(tikhonovlambda) && robustid!=GAL_FIT_ROBUST_INVALID)
    error(EXIT_FAILURE, 0, "%s: Robust (non-zero 'robustid') and "
          "regularized (non-NaN 'tikhonovlambda') polynomial fits "
          "cannot be called simultaneously. The values are: "
          "robustid=%u and tikhonovlambda=%lf)", __func__, robustid,
          tikhonovlambda);

  /* Make sure the types and lengths of each column are correct. */
  *xdata =        fit_sanity_check_col(xin,  xin, __func__);
  *ydata =        fit_sanity_check_col(yin,  yin, __func__);
  *wdata = ywht ? fit_sanity_check_col(ywht, yin, __func__) : NULL;
  (*xdata)->next = ( xin->next
                   ? fit_sanity_check_col(xin->next, xin, __func__)
                   : NULL );
}





static double
fit_polynomial_base_robust(gsl_matrix *x, gsl_vector *y, gsl_vector *c,
                           gsl_matrix *cov,
                           const gsl_multifit_robust_type *rtype)
{
  double out=NAN;
  gsl_multifit_robust_workspace *work_r;

  /* Initialize the worker and do the fit (depending on if a weight
     image was provided). */
  work_r=gsl_multifit_robust_alloc(rtype, x->size1, x->size2);
  gsl_multifit_robust(x, y, c, cov, work_r);

  /* Get the residual sum of squares, free the worker and return. */
  out=gsl_multifit_robust_statistics(work_r).sse;
  gsl_multifit_robust_free(work_r);
  return out;
}





static gal_data_t *
fit_polynomial_base(gal_data_t *xin, gal_data_t *yin,
                    gal_data_t *ywht, size_t degree,
                    uint8_t robustid, uint8_t matrixid,
                    double tikhonovlambda, double *redchisq)
{
  /* Low-level variable. */
  size_t nconst = fit_polynomial_nconst(degree, matrixid);

  /* Other variables */
  gsl_vector *c=NULL;
  double rnorm_t, snorm_t;
  double chisq=NAN, sse=NAN;
  gsl_matrix *x=NULL, *cov=NULL;
  size_t covsize[2]={nconst, nconst};
  gsl_multifit_linear_workspace *work_n;
  gal_data_t *xdata, *ydata, *wdata, *tmp, *out=NULL;

  /* For the 'y' and 'w' GSL vectors, we don't actually need to allocate
     any space, we can just use the allocated space within the
     'gal_data_t'. We can't set the pointers now because we aren't sure
     they have 'double' type yet. */
  gsl_vector yvec={yin->size, 1, NULL, NULL, 0}; /* Both have same size. */
  gsl_vector wvec={yin->size, 1, NULL, NULL, 0}; /* 'ywht' may be NULL!  */
  gsl_vector *y=&yvec, *w=&wvec; /* These have to be after the two above.*/

  /* Basic check of the inputs. */
  if(xin==NULL || xin->size==0 || yin==NULL || yin->size==0)
    error(EXIT_FAILURE, 0, "%s: the inputs are either NULL or "
          "do not contain any data (have a length of zero)", __func__);

  /* Fill all the GSL structures after a sanity check of th einput
     columns. */
  fit_polynomial_sanity_check(xin, yin, ywht, matrixid, &xdata,
                              &ydata, &wdata, degree, robustid,
                              tikhonovlambda);
  fit_polynomial_prepare(xdata, ydata, wdata, nconst,
                         &x, &c, &cov, y, w, degree, matrixid);

  /* Do the fit depending on the input arguments. */
  switch(robustid)
    {
    case GAL_FIT_ROBUST_BISQUARE:
      sse=fit_polynomial_base_robust(x, y, c, cov,
                                     gsl_multifit_robust_bisquare);
      break;
    case GAL_FIT_ROBUST_CAUCHY:
      sse=fit_polynomial_base_robust(x, y, c, cov,
                                     gsl_multifit_robust_cauchy);
      break;
    case GAL_FIT_ROBUST_FAIR:
      sse=fit_polynomial_base_robust(x, y, c, cov,
                                     gsl_multifit_robust_fair);
      break;
    case GAL_FIT_ROBUST_HUBER:
      sse=fit_polynomial_base_robust(x, y, c, cov,
                                     gsl_multifit_robust_huber);
      break;
    case GAL_FIT_ROBUST_OLS:
      sse=fit_polynomial_base_robust(x, y, c, cov,
                                     gsl_multifit_robust_ols);
      break;
    case GAL_FIT_ROBUST_WELSCH:
      sse=fit_polynomial_base_robust(x, y, c, cov,
                                     gsl_multifit_robust_welsch);
      break;
    case GAL_FIT_ROBUST_INVALID:
      if( !isnan(tikhonovlambda) ) /* Tikonov regularization. */
        {
          work_n = gsl_multifit_linear_alloc(xin->size, nconst);
          gsl_multifit_linear_svd(x, work_n);
          gsl_multifit_linear_solve(tikhonovlambda, x, y, c,
                                    &rnorm_t, &snorm_t, work_n);
        }
      else /* Linear fit. */
        {
          work_n = gsl_multifit_linear_alloc(xin->size, nconst);
          if(ywht) gsl_multifit_wlinear(x, w, y, c, cov, &chisq, work_n);
          else     gsl_multifit_linear( x,    y, c, cov, &sse,   work_n);
          gsl_multifit_linear_free(work_n);
        }
      break;
    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at '%s' to "
            "fix the problem. the 'robustid' value '%d' isn't recognize",
            __func__, PACKAGE_BUGREPORT, robustid);
    }

  /* For a check:
  {
    size_t i;
    double *ca=c->data;
    for(i=0;i<=nconst;++i) { printf("%f ", ca[i]); } printf("\n");
  } //*/

  /* Allocate the output dataset containing the fit results as first
     'gal_data_t'. */
  tmp=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &nconst, NULL, 0,
                     xin->minmapsize, xin->quietmmap, NULL, NULL,
                     NULL);
  memcpy(tmp->array, c->data, nconst*sizeof c->data);
  out=tmp;

  /* Allocate the second element of the output (the covariance matrix). */
  tmp=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 2, covsize, NULL, 0,
                     xin->minmapsize, xin->quietmmap, NULL, NULL,
                     NULL);
  memcpy(tmp->array, cov->data, nconst*nconst*sizeof cov->data);
  out->next=tmp;

  /* Calculate the reduced chi^2, see the description of same step in
     'fit_1d_linear_base'. */
  *redchisq = (isnan(chisq) ? sse : chisq) / (xdata->size-nconst);

  /* Clean up and return. */
  if(xdata->next && xdata->next!=xin->next) /* Must be before 'xdata'.*/
    { gal_data_free(xdata->next); xdata->next=NULL; }
  if(ywht && wdata!=ywht) gal_data_free(wdata);
  if(xdata!=xin) gal_data_free(xdata);
  if(ydata!=yin) gal_data_free(ydata);
  gsl_matrix_free(cov);
  gsl_matrix_free(x);
  gsl_vector_free(c);
  return out;
}





gal_data_t *
gal_fit_polynomial(gal_data_t *xin, gal_data_t *yin,
                   gal_data_t *ywht, size_t degree,
                   double *redchisq, uint8_t matrixid)
{
  return fit_polynomial_base(xin, yin, ywht, degree,
                             GAL_FIT_ROBUST_INVALID, matrixid,
                             NAN, redchisq);
}





gal_data_t *
gal_fit_polynomial_robust(gal_data_t *xin, gal_data_t *yin,
                          size_t degree, uint8_t robustid,
                          double *redchisq, uint8_t matrixid)
{
  /* Robust fitting doesn't use weights (the functions are effectively the
     weight). */
  return fit_polynomial_base(xin, yin, NULL, degree, robustid,
                             matrixid, NAN, redchisq);
}





gal_data_t *
gal_fit_polynomial_tikhonov(gal_data_t *xin, gal_data_t *yin,
                            size_t degree, double *redchisq,
                            uint8_t matrixid, double tikhonovlambda)
{
  return fit_polynomial_base(xin, yin, NULL, degree,
                             GAL_FIT_ROBUST_INVALID,
                             matrixid, tikhonovlambda, redchisq);
}





static void
fit_polynomial_estimate_2d(gal_data_t *fit, gal_data_t *xin,
                           size_t degree, gsl_vector *xvec,
                           gsl_vector *cvec, gsl_matrix *cmat,
                           gal_data_t *out)
{
  size_t i;
  double x1, x2, *xi1p, *xi2p, *xo=xvec->data;
  double *x1a, *x2a, *y=out->array, *yerr=out->next->array;

  /* Necessary allocations. */
  xi1p=gal_pointer_allocate(GAL_TYPE_FLOAT64, degree+1, 1, __func__,
                            "xi1p");
  xi2p=gal_pointer_allocate(GAL_TYPE_FLOAT64, degree+1, 1, __func__,
                            "xi2p");

  /* Fill every pixel of the input. Note that the coordinates are in FITS
     format that starts from one (as in 'gal_dimension_image_to_table' that
     created the inputs to the fit). Also, since we are filling the output
     'y' array by incrementing with 'i', we need to fill the fastest axis
     first, so the inner for-loop should be the first (FITS) coordinate. */
  if(xin->ndim==2)
    {
      i=0;
      for(x2=1; x2<xin->dsize[0]+1; ++x2)
        for(x1=1; x1<xin->dsize[1]+1; ++x1)
          {
            fit_polynomial_row_init_2d(x1, x2, degree, xi1p, xi2p, NULL);
            fit_polynomial_row_2d(xi1p, xi2p, NULL, degree, xo);
            gsl_multifit_linear_est(xvec, cvec, cmat, y+i, yerr+i);
            ++i;
          }
    }
  else
    {
      x1a=xin->array;
      x2a=xin->next->array;
      for(i=0;i<xin->size;++i)
        {
          fit_polynomial_row_init_2d(x1a[i], x2a[i], degree, xi1p,
                                     xi2p, NULL);
          fit_polynomial_row_2d(xi1p, xi2p, NULL, degree, xo);
          gsl_multifit_linear_est(xvec, cvec, cmat, y+i, yerr+i);
        }
    }

  /* Clean up. */
  free(xi1p);
  free(xi2p);
}





/* Estimate values from a polynomial fit. */
gal_data_t *
gal_fit_polynomial_estimate(gal_data_t *fit, gal_data_t *xin,
                            size_t degree)
{
  size_t i, ndim;
  size_t nconst=fit->size;
  double *y, *xi, *xo, *yerr;
  gal_data_t *xd=NULL, *out=NULL;

  /* We don't need to allocate space for the GSL vectors and matrices, we
     can just use the allocated space within the 'gal_data_t'. We can't set
     the pointers now because we aren't sure they have 'double' type
     yet. */
  gsl_vector xvec={nconst, 1, NULL, NULL, 0};
  gsl_vector cvec={nconst, 1, NULL, NULL, 0};
  gsl_matrix cmat={nconst, nconst, nconst, NULL, NULL, 0};

  /* Do the basic preparations. */
  out=fit_estimate_prepare(xin, fit, &xd, __func__);

  /* Set the pointers. */
  xi        = xd->array;
  cvec.data = fit->array;
  y         = out->array;
  yerr      = out->next->array;
  cmat.data = fit->next->array;
  xo = xvec.data = gal_pointer_allocate(GAL_TYPE_FLOAT64, nconst,
                                        0, __func__, "xvec.data");

  /* Do the estimation. */
  ndim = (xin->ndim==2 || xin->next) ? 2 : 1;
  switch(ndim)
    {

    /* 1D estimation. */
    case 1:
      for(i=0;i<xd->size;++i)
        {
          fit_polynomial_row_1d(xi, i, nconst, xo);
          gsl_multifit_linear_est(&xvec, &cvec, &cmat, y+i, yerr+i);
        }
      break;

    /* 2D estimation. */
    case 2:
      fit_polynomial_estimate_2d(fit, xin->array?xin:xd, degree, &xvec,
                                 &cvec, &cmat, out);
      break;

    /* Undefined. */
    default:
      error(EXIT_FAILURE, 0, "%s: only one or two dimensional data "
            "are currently supported", __func__);
    }

  /* Clean up and return. */
  if(xin->next && xd->next!=xin->next) gal_data_free(xd->next);
  if(xd!=xin) gal_data_free(xd); /*Must be after 'xin->next'.*/
  free(xvec.data);
  return out;
}
