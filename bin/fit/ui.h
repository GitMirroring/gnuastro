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
#ifndef UI_H
#define UI_H

/* For common options groups. */
#include <gnuastro-internal/options.h>





/* Option groups particular to this program. */
enum program_args_groups
{
  UI_GROUP_PARAMS = GAL_OPTIONS_GROUP_AFTER_COMMON,
};





/* Available letters for short options:

   a b e f g i j k l n s t p u v x y z
   A B C E G H J L O Q W X Y
*/
enum option_keys_enum
{
  /* With short-option version. */
  UI_KEY_COLUMN             = 'c',
  UI_KEY_MODEL              = 'm',
  UI_KEY_DEGREE             = 'd',
  UI_KEY_ROBUST             = 'r',
  UI_KEY_ESTIMATE           = 'e',
  UI_KEY_WEIGHT             = 'w',
  UI_KEY_ESTIMATEHDU        = 'E',
  UI_KEY_RESIDUAL           = 'R',

  /* Only with long version (start with a value 1000, the rest will be set
     automatically). */
  UI_KEY_TXTISIMG      = 1000,
  UI_KEY_WEIGHTHDU,
  UI_KEY_WEIGHTCOL,
  UI_KEY_WEIGHTTYPE,
  UI_KEY_ESTIMATECOL,
  UI_KEY_OUTTABLENOINPUT,
};





void
ui_read_check_inputs_setup(int argc, char *argv[], struct fitparams *p);

void
ui_free_report(struct fitparams *p, struct timeval *t1);

#endif
