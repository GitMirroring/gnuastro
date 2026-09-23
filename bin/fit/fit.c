/*********************************************************************
Fit - Regression analysis of the input to a certain model.
Fit is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2026-2026 Free Software Foundation, Inc.

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

#include <errno.h>
#include <error.h>

#include <gnuastro/fit.h>
#include <gnuastro/data.h>

#include <gnuastro-internal/checkset.h>

#include "main.h"

#include "fit.h"





static char *
fit_wht_nature(struct fitparams *p)
{
  switch(p->whtid)
    {
    case FIT_WHT_STD:    return "Standard deviation";
    case FIT_WHT_VAR:    return "Variance";
    case FIT_WHT_INVVAR: return "Inverse variance";
    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at '%s' "
            "to find and fix the problem. The value '%d' isn't a "
            "recognized weight type identifier", __func__,
            PACKAGE_BUGREPORT, p->whtid);
    }

  /* Control should not reach here. */
  error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at '%s' "
        "to find and fix the problem. Control should not reach the "
        "end of this function", __func__, PACKAGE_BUGREPORT);
  return NULL;
}





static void
fit_params_to_keys(struct fitparams *p, gal_data_t *fit, double redchisq)
{
  size_t i, j;
  char *kname, *kcomm;
  struct gal_fits_list_key_t *out=NULL;
  struct gal_options_common_params *cp=&p->cp;
  double *c=fit->array, *cov=fit->next?fit->next->array:NULL;

  /* Set the title and basic info (independent of the type of fit). */
  gal_fits_key_list_title_add(&out, "Regression analysis (fitting) "
                              "results", 0);
  gal_fits_key_list_add(&out, GAL_TYPE_STRING, "INPUT", 0,
                        p->inputname, 0,"Name of input file.",
                        0, NULL, 0);

  /* Add the Fitting results. */
  switch(p->model)
    {
    /* Linear with no constant */
    case FIT_MODEL_LINEAR_NO_CONSTANT:
      gal_fits_key_list_add(&out, GAL_TYPE_FLOAT64, "FITC1", 0,
                            c, 0, "C1: in y=C1*x ", 0, NULL, 0);
      gal_fits_key_list_add(&out, GAL_TYPE_FLOAT64, "FCOV11", 0,
                            c+1, 0, "Variance of C1 (only element of "
                            "cov. matrix).", 0, NULL, 0);
      break;

    /* Basic linear. */
    case FIT_MODEL_LINEAR:
      gal_fits_key_list_add(&out, GAL_TYPE_FLOAT64, "FITC0", 0,
                            c, 0, "C0: Constant in linear fit "
                            "(y=C0+C1*x).", 0, NULL, 0);
      gal_fits_key_list_add(&out, GAL_TYPE_FLOAT64, "FITC1", 0,
                            c+1, 0, "C1: Multiple of X in linear fit "
                            "(y=C0+C1*x).", 0, NULL, 0);
      gal_fits_key_list_add(&out, GAL_TYPE_FLOAT64, "FCOV11", 0,
                            c+2, 0, "Element (1,1) of covariance matrix.",
                            0, NULL, 0);
      gal_fits_key_list_add(&out, GAL_TYPE_FLOAT64, "FCOV12", 0,
                            c+3, 0, "Element (1,2)=(2,1) of covariance "
                            "matrix element.", 0, NULL, 0);
      gal_fits_key_list_add(&out, GAL_TYPE_FLOAT64, "FCOV22", 0,
                            c+4, 0, "Element (2,2) of covariance matrix.",
                            0, NULL, 0);
      break;

    /* Polynomial fit. */
    case FIT_MODEL_POLYNOMIAL:
      for(i=0;i<fit->size;++i)
        {
          if( asprintf(&kname, "FITC%zu", i)<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf in FITCxx name",
                  __func__);
          if( asprintf(&kcomm, "C%zu: polynomial constant (degree=%u).",
                       i, p->degree)<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf in FITCxx comment",
                  __func__);
          gal_fits_key_list_add(&out, GAL_TYPE_FLOAT64, kname, 1,
                                c+i, 0, kcomm, 1, NULL, 0);
        }
      for(i=0;i<fit->size;++i)
        for(j=0;j<fit->size;++j)
          {
            if( asprintf(&kname, "FCOV%zu%zu", i+1, j+1)<0 )
              error(EXIT_FAILURE, 0, "%s: asprintf in FCOVxx name",
                    __func__);
            if( asprintf(&kcomm, "Element (%zu,%zu) of covariance "
                         "matrix.", i+1, j+1)<0 )
              error(EXIT_FAILURE, 0, "%s: asprintf in FCOVxx comment",
                    __func__);
            gal_fits_key_list_add(&out, GAL_TYPE_FLOAT64, kname, 1,
                                  cov+(i*fit->size+j), 0, kcomm, 1,
                                  NULL, 0);
          }
      break;

    /* Unrecognized FIT ID. */
    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at '%s' to "
            "fix the problem. The code '%d' isn't recognized for 'fitid'",
            __func__, PACKAGE_BUGREPORT, p->model);
    }

  /* Add the reduced chi-squared. */
  gal_fits_key_list_add(&out, GAL_TYPE_FLOAT64, "FRDCHISQ", 0,
                        &redchisq, 0, "Reduced chi^2 of fit.",
                        0, NULL, 0);

  /* Reverse the last-in-first-out list to be in the same logical order we
     inserted the items here. */
  gal_fits_key_list_reverse(&out);

  /* Append this list to the end of the configuration keywords and then
     write it. */
  gal_fits_key_list_append(&cp->ckeys, out);
  gal_fits_key_write(cp->ckeys, cp->output, "0", "NONE", 1, 1);
}





/* Estimate the input based on the model. */
static gal_data_t *
fit_estimate_model(struct fitparams *p, gal_data_t *fit, gal_data_t *xin)
{
  switch(p->model)
    {
    case FIT_MODEL_LINEAR:
    case FIT_MODEL_LINEAR_NO_CONSTANT:
      return gal_fit_linear_estimate_1d(fit, xin);
      break;

    case FIT_MODEL_POLYNOMIAL:
      return gal_fit_polynomial_estimate(fit, xin, p->degree);
      break;

    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at '%s' "
            "to fix the problem. The code '%d' isn't recognized for "
            "'fitid'", __func__, PACKAGE_BUGREPORT, p->model);
      return NULL; /* Control never reaches here: to avoid warning. */
    }
}





static void
fit_correct_name(gal_data_t *data, char *name, char *unit)
{
  if(name)
    {
      if(data->name) free(data->name);
      gal_checkset_allocate_copy(name, &data->name);
    }
  if(unit)
    {
      if(data->unit) free(data->unit);
      gal_checkset_allocate_copy(unit, &data->unit);
    }
}





static void
fit_estimate_residual_nexttonull(gal_list_void_t **nexttonull,
                                 gal_data_t *input, gal_data_t *next)
{
  gal_data_t *tmp;
  if(input)
    {
      tmp=gal_list_data_last(input);
      gal_list_void_add(nexttonull, tmp);
      tmp->next=next;
    }
}





static void
fit_estimate_residual_write_1d(struct fitparams *p, gal_data_t *xin_res,
                               gal_data_t *y_est, gal_data_t *y_res)
{
  gal_list_void_t *tv, *ntn=NULL;
  gal_data_t *tmp, *oest=NULL, *ores=NULL;
  struct gal_options_common_params *cp=&p->cp;

  /* The input columns should be in the output. */
  if(p->outtablenoinput==0)
    {
      /* If a weight was given for Y, add it after Y. */
      if(p->ywht)
        fit_estimate_residual_nexttonull(&ntn, p->yin, p->ywht);

      /* In case the estimate and residual come from the same array, set
         them to be the same pointer and different otherwise. */
      if( xin_res == p->xin_est )
        {
          ores=oest=xin_res;
          fit_estimate_residual_nexttonull(&ntn, xin_res, p->yin);
        }
      else
        {
          /* For a residual, we will add the input Y here. */
          if(xin_res)
            {
              ores=xin_res;
              fit_estimate_residual_nexttonull(&ntn, xin_res, p->yin);
            }

          /* An independent estimation can have different rows from the
             input Y; so it is not necessary to add Y to 'xin_est'. */
          if(p->xin_est)
            {
              oest=p->xin_est;
              if(p->estisself)
                fit_estimate_residual_nexttonull(&ntn, p->xin_est, p->yin);
            }
        }
    }

  /* If we have a residual, that goes first. */
  if(y_res)
    {
      if(ores)
        {
          tmp=gal_list_data_last(ores);
          gal_list_void_add(&ntn, tmp); tmp->next=y_res;
        }
      else ores=y_res; /* No input columns. */
    }

  /* If we have an estimation, we should add it to the possibly existing
     output. */
  if(y_est)
    {
      if(oest)
        {
          tmp=gal_list_data_last(oest);
          gal_list_void_add(&ntn, tmp); tmp->next=y_est;
        }
      else /* No input columns. */
        {
          /* There is a residual column in the same table. */
          if( xin_res == p->xin_est )
            {
              tmp=gal_list_data_last(ores);
              gal_list_void_add(&ntn, tmp); tmp->next=y_est;
              oest=ores;
            }
          else oest=y_est;  /* No input or residual column. */
        }
    }

  /* Write the table(s). When we have one table, it is stored in 'ores'. */
  if(ores)
    gal_table_write(ores, NULL, NULL, cp->tableformat, cp->output,
                    ores==oest?"FIT":"FIT-RESIDUAL", 0, 1);
  if(oest && ores!=oest)
    gal_table_write(oest, NULL, NULL, cp->tableformat, cp->output,
                    "FIT-ESTIMATE", 0, 1);

  /* Reset all the pointers that were originally NULL so they do not
     interfere with the freeing of the various datasets. */
  for(tv=ntn; tv!=NULL; tv=tv->next)
    { tmp=(gal_data_t *)(tv->v); tmp->next=NULL; }
  gal_list_void_free(tv, 0);
}




static void
fit_estimate_residual_write_2d(struct fitparams *p, gal_data_t *y_est,
                               gal_data_t *y_res)
{
  struct gal_options_common_params *cp=&p->cp;
  struct wcsprm *wcsest=p->estisself?p->wcs_in:p->wcs_est;

  /* If a residual was requested. */
  if(y_res)
    {
      /* Residual and its error (if necessary: not same as estimate
         error) */
      y_res->wcs=p->wcs_in;
      gal_fits_img_write(y_res, cp->output, NULL, 0);
      y_res->wcs=NULL;
      if(y_res->next) /* No error when residual and estimate==self. */
        {
          y_res->next->wcs=p->wcs_in;
          gal_fits_img_write(y_res->next, cp->output, NULL, 0);
          y_res->next->wcs=NULL;
        }
    }

  /* If an estimate was requested. */
  if(y_est)
    {
      /* Estimated dataset. */
      y_est->wcs=wcsest;
      gal_fits_img_write(y_est, cp->output, NULL, 0);
      y_est->wcs=NULL;

      /* Error in estimation. */
      y_est->next->wcs=wcsest;
      gal_fits_img_write(y_est->next, cp->output, NULL, 0);
      y_est->next->wcs=NULL;
    }
}





static void
fit_estimate_residual(struct fitparams *p, gal_data_t *fit,
                      double redchisq)
{
  size_t ondim;
  double *y, *ye, *yf, *yr;
  gal_data_t *xin_res=NULL, *y_est=NULL, *y_res=NULL;

  /* When '--residual=self' and '--residual' are requested, we only need to
     do the estimation once (the residual is just the subtraction of the
     estimation and the input). Otherwise, the estimation needs to be done
     once for the input and once for the estimated input. */
  if( p->estisself && p->residual )
    {
      xin_res=p->xin_est;
      y_est=fit_estimate_model(p, fit, p->xin_est);
    }
  else /* Estimation and residual are not from the same dataset. */
    {
      /* If estimation was requested, do it. */
      if(p->xin_est)
        y_est=fit_estimate_model(p, fit, p->xin_est);

      /* If a residual was requested, do it. */
      if(p->residual)
        {
          xin_res=p->xin_r2d?p->xin_r2d:p->xin;
          y_res=fit_estimate_model(p, fit, xin_res);
        }
    }

  /* Set the metadata of the estimated dataset. */
  if(y_est)
    {
      fit_correct_name(y_est, "Y-ESTIMATE", p->input->unit);
      fit_correct_name(y_est->next, y_res?"Y-STD":"Y-ESTIMATE-STD",
                       p->input->unit);
    }

  /* Final residual dataset */
  if(p->residual)
    {
      /* When only one estimation was done (see above), we need to allocate
         another array to keep the residual. */
      y = p->xin_r2d ? p->input->array : p->yin->array;
      if(y_res==NULL)
        {
          /* Allocate and fill the array. */
          y_res=gal_data_alloc(NULL, y_est->type, y_est->ndim,
                               y_est->dsize, NULL, 0, y_est->minmapsize,
                               y_est->quietmmap, NULL, NULL, NULL);
          yr=y_res->array;
          yf=(ye=y_est->array)+y_est->size;
          do { *yr++=(*y-*ye)/(*y); y++; } while(++ye<yf);
        }
      else
        {
          /* Subtract the input from the estimate. */
          yf=(ye=y_res->array)+y_res->size;
          do { *ye=(*y-*ye)/(*y); y++; } while(++ye<yf);
        }

      /* Set the metadata. */
      fit_correct_name(y_res, "Y-RESIDUAL-FRAC", p->input->unit);
      if(y_res->next)
        fit_correct_name(y_res->next, "Y-RESIDUAL-ESTIMATE-STD",
                         p->input->unit);
    }

  /* Write the output: when there was only a single estimation. */
  ondim=y_est?y_est->ndim:y_res->ndim;
  fit_params_to_keys(p, fit, redchisq);
  switch(ondim)
    {
    case 1:
      fit_estimate_residual_write_1d(p, xin_res, y_est, y_res); break;
    case 2:
      fit_estimate_residual_write_2d(p, y_est, y_res); break;
    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at '%s' "
            "to find a fix the problem. The value %zu is not an "
            "acceptable value for 'ondim'", __func__, PACKAGE_BUGREPORT,
            ondim);
    }

  /* Inform the user if not in quiet mode. */
  if(p->cp.quiet==0) printf("Output: %s\n", p->cp.output);

  /* Clean up. */
  gal_list_data_free(y_est);
  gal_list_data_free(y_res);
}





static char *
fit_print_intro(struct fitparams *p)
{
  char *filename, *intro, *wcolstr=NULL;

  /* Set the full file name (for easy reading later!).*/
  filename=gal_fits_name_save_as_string(p->inputname, p->cp.hdu);

  /* Set the Weight column string(s). */
  if(p->ywht)
    {
      if( asprintf(&wcolstr, "\nWeight: %s [%s of Y in each row]",
                   p->weightname, fit_wht_nature(p))<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
    }

  /* Put everything into one string. */
  if( asprintf(&intro,
               "%s\n"
               "-------\n"
               "Fitting results (remove extra info with '--quiet' "
               "or '-q)\n"
               "Input file: %s with %zu elements."
               "%s",
               PROGRAM_STRING, filename, p->xin->size,
               wcolstr ? wcolstr : "")<0 )
    error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);

  /* Clean up and return. */
  free(filename);
  free(wcolstr);
  return intro;
}





static void
fit_print_linear(struct fitparams *p, gal_data_t *fit)
{
  char *intro, *funcvals;
  double redchisq=NAN, *f=fit->array;

  /* The prints depend on the user asking to be quiet or not. */
  if(p->cp.quiet)
    {
      /* The covariance matrix will only be present if a constant was
         present. After it, just print the residual (chi-squared or sum of
         squares of residuals). */
      switch(p->model)
        {
        case FIT_MODEL_LINEAR_NO_CONSTANT:
          printf("%+-.15e\n%+-.15e\n%+-.15e\n", f[0], f[1], f[2]);
          break;
        case FIT_MODEL_LINEAR:
          printf("%+-.15e %+-.15e\n"    /* The Two coefficients. */
                 "%+-20.15e %+-20.15e\n"/* First row of cov. matrix. */
                 "%+-20.15e %+-20.15e\n"/* Second row of cov. matrix.*/
                 "%+-.15e\n",           /* Residual. */
                 f[0], f[1], f[2], f[3], f[3], f[4], f[5]);
          break;
        default:
          error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at '%s' to "
                "fix the problem. The value '%u' is not recognized for "
                "'p->model' in the quiet condition", __func__,
                PACKAGE_BUGREPORT, p->model);
        }
    }

  /* Quiet was not given. */
  else
    {
      switch(p->model)
        {
        case FIT_MODEL_LINEAR_NO_CONSTANT:
          redchisq=f[2];
          if( asprintf(&funcvals,
                       "Fitting model: Y = c1 * X\n"
                       "  c1: %+-.15e\n\n"
                       "Variance of 'c1':\n"
                       "  %+-.15e\n", f[0], f[1])<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
        case FIT_MODEL_LINEAR:
          redchisq=f[5];
          if( asprintf(&funcvals,
                       "Fitting model: Y = c0 + (c1 * X)\n"
                       "  c0:  %+-.15e\n"
                       "  c1:  %+-.15e\n\n"
                       "Covariance matrix (off-diagonal are identical "
                       "same):\n"
                       "  %+-20.15e %+-20.15e\n"
                       "  %+-20.15e %+-20.15e\n", f[0], f[1], f[2], f[3],
                       f[3], f[4])<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
          break;
        default:
          error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at '%s' to "
                "fix the problem. The value '%u' is not recognized for "
                "'p->model' in the non-quiet condition", __func__,
                PACKAGE_BUGREPORT, p->model);
        }

      /* Program version, and input filenames. */
      intro=fit_print_intro(p);

      /* Final printed report.*/
      printf("%s\n\n%s\nReduced chi^2 of fit:\n  %+.15e\n", intro,
             funcvals, redchisq);

      /* Clean up. */
      free(intro);
      free(funcvals);
    }
}





static void
fit_print_polynomial(struct fitparams *p, gal_data_t *fit,
                     double redchisq)
{
  size_t i, j;
  char *intro;
  size_t nconst=fit->size;
  double *farr=fit->array, *carr=fit->next->array;

  /* Print fitted constants */
  if(p->cp.quiet)
    {
      if(p->xin_est==NULL)
        {
          for(i=0;i<nconst;++i) printf("%+-.15e ", farr[i]);
          printf("\b\n");
        }
    }
  else
    {
      /* Program version, and input filenames. */
      intro=fit_print_intro(p);

      /* Final printed report.*/
      printf("%s\n\n", intro);
      switch(p->ndim)
        {
        case 1:
          printf("Fitting model [1d]: Y = c0 + (c1 * X^1) + "
                 "(c2 * X^2) + ... (cN * X^N)\n"); break;
        case 2:
          printf("Fitting model [2d, Y=f(X1,X2)]: Y = c0 + c1.X1 "
                 "+ c2.X2 + c3.X1^2 + c4.X1.X2 + c5.X2^2 + c6.X1^3 "
                 "+ c7.X1^2.X2 + c8.X1.X2^2 + c9.X2^3 + ... "
                 "+ cn.X1^(j).X2^(d-j)\n"); break;
        }

      /* Notice for the (possible) robust function and degree of
         polynomial. */
      if(p->robustname)
        printf("  Robust function: %s\n", p->robustname);
      printf("  Degree:  %d\n", p->degree);

      /* Print the fitted values. */
      for(i=0;i<nconst;++i)
        printf("  c%zu: %s%+-.15e\n", i, i<10?" ":"", farr[i]);

      /* Print the information on the covariance matrix. */
      printf("\nCovariance matrix:\n");

      /* Clean up. */
      free(intro);
    }

  /* Print the covariance matrix (when no estimation is requested, this is
     the same for quiet or non-quiet mode). But when estimation is
     requested, they will be written in the FITS keywords of the output so
     there is no more need to have them on the command-line.*/
  if(p->cp.quiet==0 || p->xin_est==NULL)
    {
      for(i=0;i<nconst;++i)
        {
          if(p->cp.quiet==0) printf("  ");
          for(j=0;j<nconst;++j)
            printf("%+-20.15e ", carr[i*nconst+j]);
          printf("\b\n");
        }

      /* Print the chi^2. */
      if(p->cp.quiet==0)
        printf("\nReduced chi^2 of fit:\n");
      printf("%s%+-.15e\n", p->cp.quiet?"":"  ", redchisq);
    }
}





void
fit(struct fitparams *p)
{
  double redchisq=NAN;
  uint8_t matrixid, islinear=1;
  gal_data_t *residual=NULL, *fit=NULL;

  /* Do the fitting depending on the model. */
  switch(p->model)
    {
    case FIT_MODEL_LINEAR:
      fit=gal_fit_linear_1d(p->xin, p->yin, p->ywht, &redchisq);
      break;

    case FIT_MODEL_LINEAR_NO_CONSTANT:
      fit=gal_fit_linear_no_constant_1d(p->xin, p->yin, p->ywht,
                                        &redchisq);
      break;

    case FIT_MODEL_POLYNOMIAL:
      islinear=0;
      matrixid = ( p->ndim==1
                   ? GAL_FIT_MATRIX_POLYNOMIAL_1D
                   : GAL_FIT_MATRIX_POLYNOMIAL_2D );
      fit = ( p->robust
              ? gal_fit_polynomial_robust(p->xin, p->yin, p->degree,
                                          p->robust, &redchisq, matrixid)
              : gal_fit_polynomial(p->xin, p->yin, p->ywht, p->degree,
                                   &redchisq, matrixid) );
      break;

    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at '%s' to "
            "fix the problem. '%s' is not a recognized as a fit type",
            __func__, PACKAGE_BUGREPORT, p->modelname);
    }

  /* Estimate values (if requested), note that it involves writing the
     fitted parameters in the header. */
  if(p->cp.output) fit_estimate_residual(p, fit, redchisq);
  else
    {
      if(islinear) fit_print_linear(p, fit);
      else         fit_print_polynomial(p, fit, redchisq);
    }
  /* Clean up. */
  if(residual) gal_data_free(residual);
  gal_data_free(fit);
  return;
}
