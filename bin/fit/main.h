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
#ifndef MAIN_H
#define MAIN_H

/* Include necessary headers */
#include <wcslib/wcs.h>

#include <gnuastro/data.h>

#include <gnuastro-internal/options.h>

/* Progarm names.  */
#define PROGRAM_NAME   "Fit"    /* Program full name.       */
#define PROGRAM_EXEC   "astfit" /* Program executable name. */
#define PROGRAM_STRING PROGRAM_NAME" (" PACKAGE_NAME ") " PACKAGE_VERSION







/* Main program parameters structure */
struct fitparams
{
  /* From command-line */
  struct gal_options_common_params     cp; /* Common parameters.         */
  uint8_t            txtisimg;  /* Input text file is an image.          */
  char             *inputname;  /* Input filename.                       */
  gal_list_str_t     *columns;  /* Col. name or num. when input is table.*/
  char            *weightname;  /* File containing weight dataset.       */
  char             *weighthdu;  /* HDU containing weight dataset.        */
  char             *weightcol;  /* Column containing weight dataset.     */
  char            *weighttype;  /* Type of weight (e.g., STD or VAR).    */
  char           *estimatestr;  /* Name or number to estimate the fit.   */
  char           *estimatehdu;  /* HDU containing to-estimate values.    */
  gal_list_str_t *estimatecol;  /* Column containing to-estimate values. */
  char             *modelname;  /* Model name.                           */
  uint8_t              degree;  /* Degree of model to use.               */
  char            *robustname;  /* Robust function name to use.          */
  uint8_t            residual;  /* Calculate the residual.               */
  uint8_t     outtablenoinput;  /* Do not include input X,Y in 1D output.*/

  /* Internal variables. */
  uint8_t              isfits;  /* If the input is a FITS file.          */
  uint8_t               isimg;  /* If input is an image or table.        */
  gal_data_t           *input;  /* Raw input (from file).                */
  gal_data_t             *xin;  /* Independent variable(s) as columns.   */
  gal_data_t             *yin;  /* Dependent variable as a column.       */
  gal_data_t            *ywht;  /* Weight of dependent variable.         */
  gal_data_t         *xin_est;  /* Estimation independent var.(s).       */
  gal_data_t         *xin_r2d;  /* Residual on 2D independent var.(s).   */
  uint8_t           estisself;  /* The estimation is on the same input.  */
  uint8_t               whtid;  /* Code for the nature of the weight col.*/
  uint8_t               model;  /* ID of desired model to fit.           */
  size_t                 ndim;  /* Number of dimensions in the fit.      */
  uint8_t              robust;  /* ID of robust fit type.                */
  struct wcsprm       *wcs_in;  /* WCS of input (to write in output).    */
  struct wcsprm      *wcs_est;  /* WCS of estimation input.              */

  /* Output: */
  time_t              rawtime;  /* Starting time of the program.         */
};

#endif
