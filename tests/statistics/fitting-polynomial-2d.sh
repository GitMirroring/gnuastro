# Do a 2D fitting
#
# See the Tests subsection of the manual for a complete explanation
# (in the Installing gnuastro section).
#
# Original author:
#     Mohammad Akhlaghi <mohammad@akhlaghi.org>
# Contributing author(s):
#     Giacomo Lorenzetti <glorenzetti@cefca.es>
# Copyright (C) 2025-2025 Free Software Foundation, Inc.
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
prog=statistics
prog2=arithmetic
execname=../bin/$prog/ast$prog
execarith=../bin/$prog2/ast$prog2





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
if [ ! -f $execname ];  then echo "$execname not created.";  exit 77; fi
if [ ! -f $execarith ]; then echo "$execarith not created."; exit 77; fi





# Actual test script
# ==================
#
# 'check_with_program' can be something like Valgrind or an empty
# string. Such programs will execute the command if present and help in
# debugging when the developer doesn't have access to the user's system.
fitin=fitting-polynomial-2d-input.fits
$execarith 100 100 2 makenew indexonly set-i \
           i 100 / 1 + set-Y \
           i 100 % 1 + set-X \
           5 X Y x + 10 Y x + f64 set-p \
           i 0 constant 100 mknoise-sigma set-n \
           p n 100 lt nan where --output=$fitin

$check_with_program $execname $fitin --fit=polynomial-robust \
                    --fitmaxpower=2 --output=fitting-estimate.txt
rm $fitin
