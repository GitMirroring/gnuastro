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

#include <argp.h>
#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <string.h>

#include <gnuastro/fit.h>
#include <gnuastro/txt.h>
#include <gnuastro/wcs.h>
#include <gnuastro/fits.h>
#include <gnuastro/array.h>
#include <gnuastro/table.h>

#include <gnuastro-internal/timing.h>
#include <gnuastro-internal/options.h>
#include <gnuastro-internal/checkset.h>
#include <gnuastro-internal/tableintern.h>
#include <gnuastro-internal/fixedstringmacros.h>

#include "main.h"

#include "ui.h"
#include "fit.h"
#include "authors-cite.h"





/**************************************************************/
/*********      Argp necessary global entities     ************/
/**************************************************************/
/* Definition parameters for the Argp: */
const char *
argp_program_version = PROGRAM_STRING "\n"
                       GAL_STRINGS_COPYRIGHT
                       "\n\nWritten/developed by "PROGRAM_AUTHORS;

const char *
argp_program_bug_address = PACKAGE_BUGREPORT;

static char
args_doc[] = "ASTRdata";

const char
doc[] = GAL_STRINGS_TOP_HELP_INFO PROGRAM_NAME" is Gnuastro's program "
  "for regression analysis (also known as \"fitting\"). It takes a "
  "single input dataset (that can be a table or image; formatted in "
  "plain-text or FITS) and one of the pre-defined models (given to "
  "the '--model' option below) to do the analysis.\n"
  GAL_STRINGS_MORE_HELP_INFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;




















/**************************************************************/
/*********    Initialize & Parse command-line    **************/
/**************************************************************/
static void
ui_initialize_options(struct fitparams *p,
                      struct argp_option *program_options,
                      struct argp_option *gal_commonopts_options)
{
  size_t i;
  struct gal_options_common_params *cp=&p->cp;


  /* Set the necessary common parameters structure. */
  cp->program_struct     = p;
  cp->poptions           = program_options;
  cp->program_name       = PROGRAM_NAME;
  cp->program_exec       = PROGRAM_EXEC;
  cp->program_bibtex     = PROGRAM_BIBTEX;
  cp->program_authors    = PROGRAM_AUTHORS;
  cp->coptions           = gal_commonopts_options;

  /* Program-specific initializers. */
  p->degree              = GAL_BLANK_UINT8;

  /* Modify common options. */
  for(i=0; !gal_options_is_last(&cp->coptions[i]); ++i)
    {
      /* Select individually. */
      switch(cp->coptions[i].key)
        {
        case GAL_OPTIONS_KEY_SEARCHIN:
        case GAL_OPTIONS_KEY_MINMAPSIZE:
        case GAL_OPTIONS_KEY_TABLEFORMAT:
          cp->coptions[i].mandatory=GAL_OPTIONS_MANDATORY;
          break;

        case GAL_OPTIONS_KEY_LOG:
        case GAL_OPTIONS_KEY_TYPE:
          cp->coptions[i].flags=OPTION_HIDDEN;
          break;
        }

      /* Select by group. */
      switch(cp->coptions[i].group)
        {
        case GAL_OPTIONS_GROUP_TESSELLATION:
          cp->coptions[i].doc=NULL; /* Necessary to remove title. */
          cp->coptions[i].flags=OPTION_HIDDEN;
          break;
        }
    }
}





/* Parse a single option: */
error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
  struct fitparams *p = state->input;

  /* Pass 'gal_options_common_params' into the child parser. */
  state->child_inputs[0] = &p->cp;

  /* In case the user incorrectly uses the equal sign (for example
     with a short format or with space in the long format, then 'arg'
     start with (if the short version was called) or be (if the long
     version was called with a space) the equal sign. So, here we
     check if the first character of arg is the equal sign, then the
     user is warned and the program is stopped: */
  if(arg && arg[0]=='=')
    argp_error(state, "incorrect use of the equal sign ('='). For short "
               "options, '=' should not be used and for long options, "
               "there should be no space between the option, equal sign "
               "and value");

  /* Set the key to this option. */
  switch(key)
    {
    /* Read the non-option tokens (arguments): */
    case ARGP_KEY_ARG:

      /* Only one input file name is acceptable in this program. All other
         files should be given to options. */
      if(p->inputname)
        argp_error(state, "only one argument (input file) may be "
                   "given; the extra argument is '%s'", arg);
      else
        /* The user may give a shell variable that is empty! In that case
           'arg' will be an empty string! We don't want to account for such
           cases (and give a clear error that no input has been given). */
        if(arg[0]!='\0') p->inputname=arg;
      break;

    /* This is an option, set its value. */
    default:
      return gal_options_set_from_key(key, arg, p->cp.poptions, &p->cp);
    }

  /* This option/argument has been parsed successfully. */
  return 0;
}




















/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/

/* Read the desired parameters. */
static uint8_t
ui_model_id_from_string(char *name)
{
  if( !strcmp(name, "linear") )
    return FIT_MODEL_LINEAR;
  else if( !strcmp(name, "linear-no-constant") )
    return FIT_MODEL_LINEAR_NO_CONSTANT;
  else if( !strcmp(name, "polynomial") )
    return FIT_MODEL_POLYNOMIAL;
  else return FIT_MODEL_INVALID;

  /* If control reaches here, there was a bug! */
  error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at '%s' to "
        "find a fix it. Control should not have reached here",
        __func__, PACKAGE_BUGREPORT);
  return FIT_MODEL_INVALID;
}





/* Check ONLY the options. When arguments are involved, do the check
   in 'ui_check_options_and_arguments'. */
static void
ui_check_only_options(struct fitparams *p)
{
  /* Read and check the model name. */
  p->model=ui_model_id_from_string(p->modelname);
  switch(p->model)
    {
    case FIT_MODEL_LINEAR:
    case FIT_MODEL_LINEAR_NO_CONSTANT:
      if(p->degree!=GAL_BLANK_UINT8)
        error(EXIT_FAILURE, 0, "linear models are not compatible with "
              "--degree");
      if(p->robustname)
        error(EXIT_FAILURE, 0, "the linear models do not have any "
              "robust algorithm implementation. Use '--model=polynomial' "
              "with '--degree=1' instead: the linear models are highly "
              "optimized only for y = c0 + c1.X1");
      break;

    case FIT_MODEL_POLYNOMIAL:

      /* Degree is mandatory for polynomials. */
      if(p->degree==GAL_BLANK_UINT8)
        error(EXIT_FAILURE, 0, "'--degree' is necessary for polynomial "
              "model. This is the maximum power of the terms in the "
              "fitted polynomial");

      /* If the user asked for a robust fit (to remove outliers). */
      if( p->robustname )
        {
          p->robust=gal_fit_name_robust_to_id(p->robustname);
          if(p->robust==GAL_FIT_ROBUST_INVALID)
            error(EXIT_FAILURE, 0, "'%s' is not a recognized robust "
                  "algorithm name. Please see the description of "
                  "'--robust' in the manual (you can run `info %s`)",
                  p->robustname, PROGRAM_EXEC);
        }

      /* The weight is not compatible with robust methods. */
      if(p->robust && p->weightname)
        error(EXIT_FAILURE, 0, "robust methods do not take weights");
      break;

    /* The value given to '--model' could not be read. */
    case FIT_MODEL_INVALID:
      error(EXIT_FAILURE, 0, "'%s' is not a recognized model name. "
            "Please see the description of '--model' in the manual "
            "(you can run `info %s`)", p->modelname, PROGRAM_EXEC);

    /* Unexpected! */
    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at '%s' to "
            "fix the problem. The value '%d' is not recognized for "
            "'p->model'", __func__, PACKAGE_BUGREPORT, p->model);
    }

  /* If a weight is specified. */
  if(p->weightname)
    {
      /* Make sure a weight type is given. */
      if(p->weighttype==NULL)
        error(EXIT_FAILURE, 0, "no '--weight-type' given");

      /* Read the value given to '--weight-type'. */
      if(       !strcmp(p->weighttype, "std") ) p->whtid=FIT_WHT_STD;
      else if ( !strcmp(p->weighttype, "var") ) p->whtid=FIT_WHT_VAR;
      else if ( !strcmp(p->weighttype, "inv-var") )
        p->whtid=FIT_WHT_INVVAR;
      else
        error(EXIT_FAILURE, 0, "'%s' is not a recognized weight type! "
              "Please use either 'std' (standard deviation), 'var' "
              "(variance) or 'inv-var' (inverse variance)",
              p->weighttype);
    }

  /* If the estimation is to be done on the input. */
  if(p->estimatestr)
    p->estisself = p->estimatestr && !strcmp(p->estimatestr, "self");
}





static void
ui_sanity_file_check(char *name, char *hdu, char *hduoptionstr,
                     gal_list_str_t *columns, char *coloptionstr,
                     uint8_t txtisimg, uint8_t *isfits, uint8_t *isimg)
{
  int inhdutype;

  /* Check if it exists. */
  gal_checkset_check_file(name);

  /* If it is FITS, a HDU is also mandatory. */
  if( (*isfits=gal_fits_file_recognized(name)) )
    {
      /* A HDU is necessary. */
      if(hdu==NULL)
        error(EXIT_FAILURE, 0, "no HDU specified. When the input is a "
              "FITS file, a HDU must also be specified, you can use the "
              "'--hdu' ('-h') option and give it the HDU number "
              "(starting from zero), extension name, or anything "
              "acceptable by CFITSIO");

      /* Early checks purely based on the given HDU type. */
      inhdutype=gal_fits_hdu_format(name, hdu, hduoptionstr);
      *isimg=inhdutype==IMAGE_HDU;
    }

  /* Plain-text file: we rely on the user's input. */
  else *isimg=txtisimg;

  /* Now that we know if the input is an image or table (independent of
     format, we can do the necessary checks. */
  if(*isimg)
    {
      if(columns)
        error(EXIT_FAILURE, 0, "%s (hdu: %s): is a FITS image "
              "extension, but '--column' option is only relevant "
              "in tables", name, hdu);
    }
  else
    {
      /* In a table, the user should specify column(s). */
      if(columns==NULL)
        {
          if(*isfits)
            error(EXIT_FAILURE, 0, "%s (hdu: %s): is a table, but no "
                  "columns were given. Please use the '%s' option "
                  "to specify the column(s) to use", name, hdu,
                  coloptionstr);
          else
            error(EXIT_FAILURE, 0, "%s: is a table, but no columns "
                  "specified. Please use the '%s' option to specify "
                  "the column(s) to use. In case you want this "
                  "plain-text file to be read as an image, please run "
                  "with '--txt-is-img'", name, coloptionstr);
        }
    }
}





static gal_data_t *
ui_prepare_empty_2d(struct gal_options_common_params *cp,
                    gal_data_t *inimg, size_t *dsize)
{
  size_t tmp=1;
  gal_data_t *out;

  /* Allocate an empty dataset (that will be given to
     'gal_fit_polynomial_estimate'): we are giving a non-null first
     argument, so no 'array' is not allocated and the pointer is written in
     its place. Immediately afterwards, the 'array' pointer is set to
     NULL. If we do not do this trick (and pass a NULL as first argument),
     this function will allocate a space the size of the input (a big waste
     of RAM): the actual values inside it are not used later, only its size
     is important. */
  out=gal_data_alloc(&tmp, GAL_TYPE_FLOAT64, 2, inimg?inimg->dsize:dsize,
                     inimg?inimg->wcs:NULL, 0, cp->minmapsize,
                     cp->quietmmap, "FITESTIMATE", NULL, NULL);
  out->array=NULL;
  return out;
}





static gal_data_t *
ui_sanity_file_read(struct fitparams *p, char *name, char *hdu,
                    gal_list_str_t *lines, char *hduoptionstr,
                    gal_list_str_t *columns, uint8_t isfits,
                    uint8_t isimg)
{
  int nwcs;
  gal_data_t *out;
  size_t ncols, *dsize;
  struct gal_options_common_params *cp=&p->cp;

  /* Read the given file. In the following two scenarios, we need to
     read the input's data:
       - The main input:
         - If it was a file: identified by the pointer of its name in
           comparison with the one in 'p').
         - If the stdin 'lines' was given.
       - The estimate dataset when it is a table ('isimg==0'). */
  if( p->inputname==name || lines || isimg==0 )
    {
      out = ( isimg
              ? gal_array_read(name, hdu, lines,
                               cp->minmapsize, cp->quietmmap,
                               hduoptionstr)
              : gal_table_read(name, hdu, lines, columns,
                               cp->searchin, cp->ignorecase,
                               cp->numthreads, cp->minmapsize,
                               cp->quietmmap, NULL, hduoptionstr) );
      if(isimg)
        p->wcs_in = gal_wcs_read(name, hdu, cp->wcslinearmatrix,
                                 -1, -1, &nwcs, hduoptionstr);
    }

  /* This is the estimation input that is an image.*/
  else
    {
      /* We do not need to read the file, we just need its size. */
      gal_array_info_size(name, hdu, hduoptionstr, &dsize);
      out=ui_prepare_empty_2d(&p->cp, NULL, dsize);
      free(dsize);

      /* If the estimation dataset is an image, keep its WCS for the output
         (if it doesn't have any, this function will return NULL). */
      if(gal_fits_hdu_format(name, hdu, hduoptionstr)==IMAGE_HDU)
        p->wcs_est = gal_wcs_read(name, hdu, cp->wcslinearmatrix,
                                  -1, -1, &nwcs, hduoptionstr);
    }

  /* Correct any redundant dimensions from the input. For example some
     radio astronomy software make 2D image with four dimensions: the last
     two are just one element thick! */
  out->ndim=gal_dimension_remove_extra(out->ndim, out->dsize, p->wcs_in);

  /* Inputs larger than 2D are not supported. */
  if(out->ndim>2)
    error(EXIT_FAILURE, 0, "input has %zu dimensions; but Fit currently "
          "only supports 1D and 2D fits", out->ndim);

  /* When the input was a table, but in a vector column, we currently treat
     it as a 2D image. */
  if(out->next && out->ndim==2)
    error(EXIT_FAILURE, 0, "the first given column is a vector column "
          "which is currently treated as a 2D image (where the two "
          "independent variables are the positions of every value within "
          "it and the dependent variable is the value). Thefore, only a "
          "single column should be given. If you intended to fit vector "
          "columns in a different way, please contact us at '%s'",
          PACKAGE_BUGREPORT);

  /* In case the user's input does not have metadata, add them here. */
  if(out->ndim==1)
    {
      ncols=gal_list_data_number(out);
      switch(ncols)
        {
        case 2:
          if(out->name==NULL) free(out->name);
          gal_checkset_allocate_copy("X", &out->name);
          if(out->next->name==NULL) free(out->next);
          gal_checkset_allocate_copy("Y", &out->next->name);
          break;
        case 3:
          if(out->name==NULL) free(out->name);
          gal_checkset_allocate_copy("X1", &out->name);
          if(out->next->name==NULL) free(out->next->name);
          gal_checkset_allocate_copy("X2", &out->next->name);
          if(out->next->next->name==NULL) free(out->next->next->name);
          gal_checkset_allocate_copy("Y", &out->next->next->name);
          break;
        default:
          error(EXIT_FAILURE, 0, "there can only be 2 or 3 input columns "
                "to read, but there are %zu columns", ncols);
        }
    }

  /* Return the dataset. */
  return out;
}





static void
ui_check_options_and_arguments(struct fitparams *p)
{
  /* If an input file name was given and if it was a FITS file, that a HDU
     is also given. */
  if(p->inputname)
    ui_sanity_file_check(p->inputname, p->cp.hdu, "--hdu", p->columns,
                         "--column", p->txtisimg, &p->isfits, &p->isimg);

  /* When no input name was given (input is expected from stdin), the
     '--column' option specifies if we have a 1D or 2D input.*/
  else
    if(p->columns==NULL && p->txtisimg==0)
      error(EXIT_FAILURE, 0, "no '--column' given. Note that standard "
            "input (from pipes) is read as a plain-text table by "
            "default. If it should be read as a plain-text table, "
            "please use '--txt-is-img'");
}




















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
static void
ui_read_raw_input(struct fitparams *p)
{
  size_t ncols;
  gal_list_str_t *lines;
  struct gal_options_common_params *cp=&p->cp;

  /* Read and do basic checks on the input. */
  lines=gal_options_check_stdin(p->inputname, cp->stdintimeout, "input");
  p->input=ui_sanity_file_read(p, p->inputname, cp->hdu, lines, "--hdu",
                               p->columns, p->isfits, p->isimg);

  /* If we have a table input, we can either have two or three
     columns. Otherwise, the program should abort. This is done here (after
     reading the columns; not 'ui_check_only_options'), because the user
     may have used regular expressions to find their desired columns. */
  if(p->input->next)
    {
      ncols=gal_list_data_number(p->input);
      if(ncols==1)
        error(EXIT_FAILURE, 0, "only one column was read. For a 1D fit "
              "two columns are needed (for the X and Y) and for a 2D "
              "fit, three columns: X1, X2 and Y.");
      if(ncols>3)
        error(EXIT_FAILURE, 0, "the '--column' (or '-c') option should "
              "be used to give only one or two input columns, but %zu "
              "columns were produced (note that this option also "
              "accepts regular expressions for column selection)", ncols);
    }

  /* All the inputs should have the same dimentions (the user may
     mistakenly give a vector column: which is 2D). */
  if(p->input->next && p->input->ndim!=p->input->next->ndim)
    {
      error(EXIT_FAILURE, 0, "all the inputs should have the same "
            "number of dimensions, the first and second columns have "
            "'%zu' and '%zu' dimensions respectively", p->input->ndim,
            p->input->next->ndim);
      if(p->input->next->next
         && p->input->ndim!=p->input->next->next->ndim)
        error(EXIT_FAILURE, 0, "all the inputs should have the same "
              "number of dimensions, the first and third columns have "
              "'%zu' and '%zu' dimensions respectively", p->input->ndim,
              p->input->next->next->ndim);
    }

  /* The number of fitting dimensions can now be determined. If the input
     is 2D, then we have a 2D fit. Otherwise (when the input is a 1D
     column), the fit is 1D when there is only two columns and 2D when
     there are three columns. */
  p->ndim = ( p->input->ndim==2
              ? 2
              : p->input->next->next ? 2 : 1 );

  /* Check based on models. */
  if(p->ndim!=1
     && (    p->model==FIT_MODEL_LINEAR
          || p->model==FIT_MODEL_LINEAR_NO_CONSTANT) )
    error(EXIT_FAILURE, 0, "the '%s' model is only defined on 1D inputs. "
          "Instead, use the 'polynomial' model with '--degree=1'",
          p->model==FIT_MODEL_LINEAR?"linear":"linear-no-constant");
}





static void
ui_read_raw_weight(struct fitparams *p)
{
  char *wname;
  gal_list_str_t *wcol=NULL;
  struct gal_options_common_params *cp=&p->cp;

  /* Make sure weight column is given if we have a table input. */
  if(!p->isimg && p->weightcol==NULL)
    error(EXIT_FAILURE, 0, "please use '--weight-col' to specify the "
          "column name or number (counting from 1) that contains the "
          "weight values");

  /* The name of file to use for the weights. */
  wname=strcmp(p->weightname, "samefile") ? p->weightname : p->inputname;

  /* Read the weights dataset. */
  gal_list_str_add(&wcol, p->weightcol, 0);
  p->ywht = ( p->isimg
              ? gal_array_read(wname, p->weighthdu, NULL,
                               cp->minmapsize, cp->quietmmap,
                               "--weight-hdu")
              : gal_table_read(wname, p->weighthdu, NULL, wcol,
                               cp->searchin, cp->ignorecase,
                               cp->numthreads, cp->minmapsize,
                               cp->quietmmap, NULL,
                               "--weight-hdu") );
  gal_list_str_free(wcol, 0);

  /* Only a single dataset should have been read. */
  if(p->ywht->next)
    error(EXIT_FAILURE, 0, "only a single weight data set (image "
          "or column) should be given, not %zu data sets; the given "
          "value to '-weight-col' was '%s' in '%s' (note that this "
          "option also accepts regular expressions for column "
          "selection)" ,
          gal_list_data_number(p->ywht->next), p->weightcol, wname);

  /* It must have the same size as the input. */
  if( gal_dimension_is_different(p->input, p->ywht) )
    error(EXIT_FAILURE, 0, "the weight and input datasets must have the "
          "same dimensions and number of elements");

  /* If it is a 1D column and doesn't have any name give it a name . */
  if( p->ywht->ndim==1 && p->ywht->name==NULL )
    gal_checkset_allocate_copy("Y-WHT", &p->ywht->name);
}





static void
ui_prepare_input_wht_cols_clean(struct fitparams *p)
{
  size_t i=0;
  void *arr;
  int anyblank=0, hasblank[4]={0};
  gal_data_t *tmp, *col, *lastin=gal_list_data_last(p->input);

  /* Add the (potential) weight column into the list temporarily. Note that
     when no weight is given 'ywht==NULL', so it causes no problems.*/
  lastin->next=p->ywht;

  /* Check if we have blank values. We will parse all of them to also
     update all their flags. */
  for(col=p->input; col!=NULL; col=col->next)
    {
      /* Blank-related checks. */
      hasblank[i]=gal_blank_present(col, 1);
      anyblank+=hasblank[i];

      /* Convert this column to 'double'. */
      if(col->type!=GAL_TYPE_FLOAT64)
        {
          /* Do the conversion in a temporary dataset. */
          tmp=gal_data_copy_to_new_type(col, GAL_TYPE_FLOAT64);

          /* Preserve the 'array' pointer of 'dtmp' and replace it with
             'tmp->array' to free. */
          arr=tmp->array;
          tmp->type=col->type;
          tmp->array=col->array;
          gal_data_free(tmp);

          /* Correct the respctive elements in this input column. */
          col->array=arr;
          col->type=GAL_TYPE_FLOAT64;
        }
    }

  /* Remove all possible blank values. */
  if(anyblank)
    {
      /* If an estimation is requested with a value of 'self', then we
         should keep the input values here (before possibly removing any
         blanks later). */
      if( p->estisself )
        {
          /* If there are any blanks in independent variable column(s), we
             need to copy the input since the blanks will be removed in the
             next step.  */
          if( hasblank[0] || (p->ndim==2 && hasblank[1]) )
            {
              /* The copying also preserves the 'next' pointer, that we
                 need to set to NULL if it is 1D. */
              p->xin_est=gal_data_copy(p->input);
              p->xin_est->next = ( p->ndim==2
                                   ? gal_data_copy(p->input->next)
                                   : NULL );
            }
        }

      /* Remove all the rows with blank values. */
      gal_blank_remove_rows(p->input, NULL, 0);
    }

  /* Remove the weight dataset and set the other pointers. */
  lastin->next=NULL;
}





static void
ui_prepare_input_wht(struct fitparams *p)
{
  double *d, *df;

  switch(p->ndim)
    {
    /* 1D input. */
    case 1:
      ui_prepare_input_wht_cols_clean(p);
      p->xin=p->input;
      p->yin=p->xin->next;
      p->xin->next=NULL; /* Has to be after 'p->yin'. */
      if(p->estisself && p->xin_est==NULL) p->xin_est=p->xin;
      break;

    /* 2D input. */
    case 2:

      /* If the input is already prepared in 1D columns, then no processing
         is necessary and we can just set the pointers. */
      switch(p->input->ndim)
        {

        case 1: /* Input was from a table. */
          ui_prepare_input_wht_cols_clean(p);
          p->xin=p->input;
          p->yin=p->xin->next->next;
          p->xin->next->next=NULL; /* Has to be after 'p->yin'. */
          if(p->estisself && p->xin_est==NULL) p->xin_est=p->xin;
          break;

        case 2: /* Input was from an image. */

          /* We should extract the coordinates and values of both the image
             and the weight (if it was given). So we'll temporarily append
             the weight input to the inputs: note that if 'p->ywht==NULL',
             it has no effect. */
          p->input->next=p->ywht; /* no effect when 'p->ywht==NULL'. */
          p->xin=gal_dimension_image_to_table(p->input);
          p->input->next=NULL; /* Remove any potential next element. */

          /* Clean up and set the final pointers; order is important
             here. */
          gal_data_free(p->ywht); /* We have it as a table now. */
          p->yin=p->xin->next->next;
          p->xin->next->next=NULL;
          p->ywht=p->yin->next;
          p->yin->next=NULL;

          /* If the estimate is to be done on the same input dataset,
             create it. */
          if(p->estisself)
            p->xin_est=ui_prepare_empty_2d(&p->cp, p->input, NULL);

          /* When a residual is requested with a 2D image (this 'case'),
             then we will put an empty copy of the input in 'xin_r2d'. */
          if(p->residual)
            p->xin_r2d=ui_prepare_empty_2d(&p->cp, p->input, NULL);
          break;

        default: /* Unexpected. */
          error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at '%s' "
                "to fix the problem. The value '%zu' is not expected for "
                "'p->input->ndim'", __func__, PACKAGE_BUGREPORT,
                p->input->ndim);
        }
      break;

    /* Only 1D and 2D fits are available. */
    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at '%s' to "
            "fix the problem. The value '%zu' is not expected for "
            "'p->ndim'", __func__, PACKAGE_BUGREPORT, p->ndim);
    }

  /* In case there are no usable elements, inform the user and abort. */
  if(p->xin->size==0)
    error(EXIT_FAILURE, 0, "no usable information in the inputs (for "
          "example all values were NaN/blank)");

  /* The weight should be the inverse-variance. */
  if(p->ywht)
    {
      p->ywht=gal_data_copy_to_new_type_free(p->ywht, GAL_TYPE_FLOAT64);
      d=p->ywht->array;
      switch(p->whtid)
        {
        case FIT_WHT_STD:
          df=d+p->ywht->size; do *d=1/(*d * *d); while(++d<df); break;
        case FIT_WHT_VAR:
          df=d+p->ywht->size; do *d = 1 / *d;    while(++d<df); break;
        case FIT_WHT_INVVAR:
          /* This is the expected, no action necessary! */ break;
        default:
          error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at '%s' "
                "to fix the problem. The value '%u' is not recognized "
                "for 'p->weighttype'", __func__, PACKAGE_BUGREPORT,
                p->whtid);
        }
    }
}





static void
ui_read_estimate_numbers(struct fitparams *p, double *evals)
{
  size_t one=1;
  gal_data_t *tmp;

  /* First number. */
  tmp=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &one,
                     NULL, 0, -1, 1, NULL, NULL, NULL);
  ((double *)(tmp->array))[0]=evals[0];
  p->xin_est=tmp;

  /* Second number goes as a new node in the list (just like the
     inputs). */
  if( !isnan(evals[1]) )
    {
      /* If the fit is 1D, this is extra. */
      if(p->ndim==1)
        error(EXIT_FAILURE, 0, "a 1D fit was requested, but two "
              "numbers were given to '--estimate'");

      /* Write the number. */
      tmp=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &one,
                         NULL, 0, -1, 1, NULL, NULL, NULL);
      ((double *)(tmp->array))[0]=evals[1];
      p->xin_est->next=tmp;
    }
}





static void
ui_read_estimate_dataset(struct fitparams *p)
{
  gal_data_t *tmp;
  uint8_t isfits=0, isimg=0;
  char *hduoptionstr="--estimate-hdu";
  struct gal_options_common_params *cp=&p->cp;

  /* Basic checks on the given file and options. */
  ui_sanity_file_check(p->estimatestr, p->estimatehdu, hduoptionstr,
                       p->estimatecol, "--estimate-col", p->txtisimg,
                       &isfits, &isimg);
  p->xin_est=ui_sanity_file_read(p, p->estimatestr, p->estimatehdu,
                                 NULL, hduoptionstr, p->estimatecol,
                                 isfits, isimg);

  /* Checks and prepartions are based on dimensionality of the fit. */
  switch(p->ndim)
    {
    case 1: /* Sanity checks for a 1D fit. */
      if(p->xin_est->ndim!=1)
        error(EXIT_FAILURE, 0, "the estimation dataset has %zu "
              "dimensions, but a 1D fit was requested", p->xin_est->ndim);
      if(p->xin_est->next)
        error(EXIT_FAILURE, 0, "more than one estimation column is "
              "given, but a 1D fit was requested");
      break;

    case 2: /* Sanity checks for a 2D fit. */
      if(p->xin_est->next==NULL && p->xin_est->ndim!=2)
        error(EXIT_FAILURE, 0, "for a 2D fit estimation, the independent "
              "variables must either be a 2D image or two columns");

      /* Prepare an empty dataset for the estimation. */
      if(p->xin_est->ndim==2)
        {
          tmp=ui_prepare_empty_2d(cp, p->xin_est, NULL);
          gal_data_free(p->xin_est);
          p->xin_est=tmp;
        }
      break;

    default: /* Unexpected number! */
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at '%s' to "
            "fix the problem. The value '%zu' is not expected for "
            "'p->ndim'", __func__, PACKAGE_BUGREPORT, p->ndim);
    }
}




static void
ui_read_estimate(struct fitparams *p)
{
  size_t i;
  gal_list_str_t *stmp, *slist;
  double *dptr, evals[2]={NAN,NAN};

  /* First, we need to check if the given string can be read as a
     coma-seaprated list of doubles. */
  i=0;
  slist=gal_options_parse_csv_strings_to_list(p->estimatestr, NULL, 0);
  for(stmp=slist; stmp!=NULL; stmp=stmp->next)
    {
      /* Crash if we have more than two elements. */
      if(i>1)
        error(EXIT_FAILURE, 0, "at most two comma-separated values "
              "can be given to '--estimate'. The given value was '%s'",
              p->estimatestr);

      /* Set the pointer to read/write. */
      dptr=evals+i++;
      if( gal_type_from_string((void **)(&dptr), stmp->v,
                               GAL_TYPE_FLOAT64) )
        { *dptr=NAN; break; } /* !=0: could not be read as a number! */
    }

  /* If they were numbers, we need one function and another if not (it was
     a dataset. */
  if( isnan(evals[0]) ) ui_read_estimate_dataset(p);
  else                  ui_read_estimate_numbers(p, evals);

  /* If an estimate and residual are requested at the same time, they have
     to have the same number of dimensions. */
  if(p->residual && p->input->ndim!=p->xin_est->ndim)
    error(EXIT_FAILURE, 0, "the estimate does not have the same format "
          "as the input (for example the input is a 2D image, but the "
          "estimate is a number or table). This is important when both "
          "the estimate and residual are requested in one run");
}





static void
ui_output(struct fitparams *p)
{
  struct gal_options_common_params *cp=&p->cp;

  /* An output is only necessary when we have a residual or estimate.*/
  if(p->residual || p->xin_est)
    {
      /* Set the output file name. For a 1D fit, it is fine to have
         'output==NULL' because it will be printed to standrad output. If
         the user wants to save the printed output, they will give
         '--output'. */
      if(cp->output)
        {
          /* If a plain-text file is requested, and the user wants both the
             residual and an estimate, the output should be a FITS file. */
          if(   p->ndim==1
                && p->residual
                && p->estimatestr
                && p->estisself==0
                && gal_fits_file_recognized(cp->output)==0 )
            error(EXIT_FAILURE, 0, "both a residual and non-self "
                  "estimate have been requested on a 1D fit. Therefore "
                  "two output tables will be created. However, you have "
                  "requested a plain-text output format that does not "
                  "support multiple tables. The only possible output "
                  "for this scenario is FITS (with a '.fits' suffix). "
                  "You can later extract your desired HDU into a "
                  "plain-text file with Gnuastro's Table program "
                  "('asttable in.fits --hdu=1 --output=out.txt')");
        }
      else /* No output name given. */
        cp->output = gal_checkset_automatic_output(cp,
                                      p->inputname,
                                      ( (p->estimatestr && p->residual)
                                        ? "-fit.fits"
                                        : ( p->estimatestr
                                            ? "-fit-estimate.fits"
                                            : "-fit-residual.fits" ) ) );

      /* Make sure the location is writable and that the file does not
         already exist. */
      gal_checkset_writable_remove(cp->output, NULL, cp->keep,
                                   cp->dontdelete);
    }

  /* When no output should be built, free any existing output name. */
  else if(cp->output) { free(cp->output); cp->output=NULL; }
}





static void
ui_preparations(struct fitparams *p)
{
  /* Read the input and (possible) weight. In the two "raw" functions, we
     just read the raw files and corresponding sanity checks: no
     re-arrangement is done on them (for preparing the fit). */
  ui_read_raw_input(p);
  if(p->weightname) ui_read_raw_weight(p);

  /* Make the necessary corrections and set the pointers for the input and
     weight datasets. */
  ui_prepare_input_wht(p);

  /* In case an estimation was requested, and was not 'self' read it. If
     it was 'self', it was set 'ui_prepare_input_wht'. */
  if(p->estimatestr && p->xin_est==NULL)
    ui_read_estimate(p);

  /* Set the output file name. */
  ui_output(p);
}



















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/
void
ui_read_check_inputs_setup(int argc, char *argv[], struct fitparams *p)
{
  struct gal_options_common_params *cp=&p->cp;


  /* Include the parameters necessary for argp from this program ('args.h')
     and for the common options to all Gnuastro ('commonopts.h'). We want
     to directly put the pointers to the fields in 'p' and 'cp', so we are
     simply including the header here to not have to use long macros in
     those headers which make them hard to read and modify. This also helps
     in having a clean environment: everything in those headers is only
     available within the scope of this function. */
#include <gnuastro-internal/commonopts.h>
#include "args.h"


  /* Initialize the options and necessary information. */
  ui_initialize_options(p, program_options, gal_commonopts_options);


  /* Read the command-line options and arguments. */
  errno=0;
  if(argp_parse(&thisargp, argc, argv, 0, 0, p))
    error(EXIT_FAILURE, errno, "parsing arguments");


  /* Read the configuration files and set the common values. */
  gal_options_read_config_set(&p->cp);


  /* Sanity check only on options. */
  ui_check_only_options(p);


  /* Print the option values if asked. Note that this needs to be done
     after the option checks so un-sane values are not printed in the
     output state. */
  gal_options_print_state(&p->cp);


  /* Prepare all the options as FITS keywords to write in output later. */
  gal_options_as_fits_keywords(&p->cp);


  /* Check that the options and arguments fit well with each other. Note
     that arguments don't go in a configuration file. So this test should
     be done after (possibly) printing the option values. */
  ui_check_options_and_arguments(p);


  /* Read/allocate all the necessary starting arrays. */
  ui_preparations(p);
}




















/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
ui_free_report(struct fitparams *p, struct timeval *t1)
{
  /* Free the allocated arrays. */
  free(p->cp.hdu);
  free(p->cp.output);
  gal_wcs_free(p->wcs_in);
  gal_wcs_free(p->wcs_est);
  gal_list_data_free(p->yin);
  gal_list_data_free(p->ywht);
  gal_list_data_free(p->input);
  gal_list_str_free(p->columns, 1);

  /* These may not be allocated. */
  if(p->xin!=p->input) gal_list_data_free(p->xin);
  if(p->xin!=p->xin_est) gal_list_data_free(p->xin_est);
  if(p->estimatecol!=p->columns) gal_list_str_free(p->estimatecol, 1);

  /* Print the final message. */
  if(!p->cp.quiet)
    gal_timing_report(t1, PROGRAM_NAME" finished in: ", 0);
}
