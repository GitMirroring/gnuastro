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
#ifndef ARGS_H
#define ARGS_H






/* Array of acceptable options. */
struct argp_option program_options[] =
  {
    {
      "columns",
      UI_KEY_COLUMN,
      "STR",
      0,
      "Column name(s) or number(s; counting from 1).",
      GAL_OPTIONS_GROUP_INPUT,
      &p->columns,
      GAL_TYPE_STRLL,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_parse_csv_strings_append
    },
    {
      "weight",
      UI_KEY_WEIGHT,
      "STR/FLT",
      0,
      "Weight ('samefile' if same file as in input).",
      GAL_OPTIONS_GROUP_INPUT,
      &p->weightname,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "weight-hdu",
      UI_KEY_WEIGHTHDU,
      "STR/INT",
      0,
      "HDU of the weights (in '--weight').",
      GAL_OPTIONS_GROUP_INPUT,
      &p->weighthdu,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "weight-col",
      UI_KEY_WEIGHTCOL,
      "STR/INT",
      0,
      "Column of the weights (only for 1D inputs).",
      GAL_OPTIONS_GROUP_INPUT,
      &p->weightcol,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "weight-type",
      UI_KEY_WEIGHTTYPE,
      "STR",
      0,
      "Weight type: 'std', 'var' or 'inv-var'.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->weighttype,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "estimate",
      UI_KEY_ESTIMATE,
      "STR/FLT",
      0,
      "Estimate fit (column number, file or 'self').",
      GAL_OPTIONS_GROUP_INPUT,
      &p->estimatestr,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "estimate-hdu",
      UI_KEY_ESTIMATEHDU,
      "STR/INT",
      0,
      "HDU containing the --fitestimate values.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->estimatehdu,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "estimate-col",
      UI_KEY_ESTIMATECOL,
      "STR/INT",
      0,
      "Estimate column(s) name(s) or number(s).",
      GAL_OPTIONS_GROUP_INPUT,
      &p->estimatecol,
      GAL_TYPE_STRLL,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET,
      gal_options_parse_csv_strings_append
    },
    {
      "txt-is-img",
      UI_KEY_TXTISIMG,
      0,
      0,
      "Read plain-text input as image, not table.",
      GAL_OPTIONS_GROUP_INPUT,
      &p->txtisimg,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },




    /* Output options. */
    {
      "residual",
      UI_KEY_RESIDUAL,
      0,
      0,
      "Residual of the fit in the output.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->residual,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "out-table-no-input",
      UI_KEY_OUTTABLENOINPUT,
      0,
      0,
      "No inputs in table outputs.",
      GAL_OPTIONS_GROUP_OUTPUT,
      &p->outtablenoinput,
      GAL_OPTIONS_NO_ARG_TYPE,
      GAL_OPTIONS_RANGE_0_OR_1,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },



    /* Fitting parameters. */
    {
      0, 0, 0, 0,
      "Customizing the fit",
      UI_GROUP_PARAMS
    },
    {
      "model",
      UI_KEY_MODEL,
      "STR",
      0,
      "'polynomial', 'linear', 'linear-no-constant'.",
      UI_GROUP_PARAMS,
      &p->modelname,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "degree",
      UI_KEY_DEGREE,
      "INT",
      0,
      "Degree of the model (e.g., in a polynomial).",
      UI_GROUP_PARAMS,
      &p->degree,
      GAL_TYPE_UINT8,
      GAL_OPTIONS_RANGE_GE_0,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },
    {
      "robust",
      UI_KEY_ROBUST,
      "STR",
      0,
      "The robust function name to use.",
      UI_GROUP_PARAMS,
      &p->robustname,
      GAL_TYPE_STRING,
      GAL_OPTIONS_RANGE_ANY,
      GAL_OPTIONS_NOT_MANDATORY,
      GAL_OPTIONS_NOT_SET
    },



    {0}
  };





/* Define the child argp structure
   -------------------------------

   NOTE: these parts can be left untouched.*/
struct argp
gal_options_common_child = {gal_commonopts_options,
                            gal_options_common_argp_parse,
                            NULL, NULL, NULL, NULL, NULL};

/* Use the child argp structure in list of children (only one for now). */
struct argp_child
children[]=
{
  {&gal_options_common_child, 0, NULL, 0},
  {0, 0, 0, 0}
};

/* Set all the necessary argp parameters. */
struct argp
thisargp = {program_options, parse_opt, args_doc, doc, children, NULL, NULL};
#endif
