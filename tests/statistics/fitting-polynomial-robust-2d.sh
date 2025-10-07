# Do a robust polynomial fit to the given data
#
# See the Tests subsection of the manual for a complete explanation
# (in the Installing gnuastro section).
#
# Original author:
#     Giacomo Lorenzetti <glorenzetti@cefca.es>
# Contributing author(s):
# Copyright (C) 2025-2025 Free Software Foundation, Inc.
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.





# The test is not yet implemented, as these features are not yet
# implemented in aststatistics.
# The lines below help setting up a testing environment: they create a
# table with points in a 1000x1000 image with values along '5+xy'.
noisethresh=50
astarithmetic 1000 1000 2 makenew indexonly set-i i 1000 / \
               set-Y i 1000 % set-X 5 X Y x + \
               f64 set-p i 0 constant 100 mknoise-sigma set-n p \
               n $noisethresh lt nan where set-V V to-1d Y to-1d X to-1d \
               --writeall --output=raw-table.fits
asttable raw-table.fits --noblankend=_all -oinput.fits

# While setting up a testing program, remember to chain the first two
# columns, while isolating the third one with a col->next->next=NULL before
# invoking gal_fit_polynomial.
