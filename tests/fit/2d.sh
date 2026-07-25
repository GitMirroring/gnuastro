# 2D fitting
#
# See the Tests subsection of the manual for a complete explanation
# (in the Installing gnuastro section).
#
# Original author:
#     Mohammad Akhlaghi <mohammad@akhlaghi.org>
# Contributing author(s):
# Copyright (C) 2026-2026 Free Software Foundation, Inc.
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.





# Preliminaries
# =============
#
# Set the variables (The executable is in the build tree). Do the
# basic checks to see if the executable is made or if the defaults
# file exists (basicchecks.sh is in the source tree).
prog=fit
execname=../bin/$prog/ast$prog
execarith=../bin/arithmetic/astarithmetic




# Skip?
# =====
#
# If the dependencies of the test don't exist, then skip it. There are two
# types of dependencies:
#
#   - The executable was not made (for example due to a configure option),
#
#   - The input data was not made (for example the test that created the
#     data file failed).
if [ ! -f $execname  ]; then echo "$execname not created."; exit 77; fi
if [ ! -f $execarith ]; then echo "$execarith not created.";  exit 77; fi





# Actual test script
# ==================
#
# 'check_with_program' can be something like Valgrind or an empty
# string. Such programs will execute the command if present and help in
# debugging when the developer doesn't have access to the user's system.
output=fit1d.fits
input=fit1d-input.fits
export GSL_RNG_SEED=1788653609
$execarith 100 100 2 makenew indexonly set-i \
           i 100 % 1 + set-X1 \
           i 100 / 1 + set-X2 \
           5 10 X2 x + X1 X2 x + f64 set-Yraw \
           Yraw sqrt set-Ystd \
           Yraw Ystd mknoise-sigma set-Ynoised \
           i 50 constant 100 mknoise-uniform set-rand \
           Ynoised rand 60 gt nan where set-Y \
           Ystd Y --writeall --envseed --output=$input
$check_with_program $execname $input -cX,Y --model=polynomial --degree=2 \
                    --weight=samefile --weight-col=Y-STD \
                    --estimate=self --residual --output=$output
rm $input
