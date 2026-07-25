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
#ifndef FIT_H
#define FIT_H

/* Definitions. */
enum fit_models
  {
    FIT_MODEL_INVALID,   /* Invalid (=0 by C standard). */
    FIT_MODEL_LINEAR,
    FIT_MODEL_LINEAR_NO_CONSTANT,
    FIT_MODEL_POLYNOMIAL,
  };

enum fit_weight_types
  {
    FIT_WHT_INVALID,   /* Invalid (=0 by C standard). */
    FIT_WHT_STD,
    FIT_WHT_VAR,
    FIT_WHT_INVVAR,
  };

void
fit(struct fitparams *p);

#endif
