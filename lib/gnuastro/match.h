/*********************************************************************
match -- Functions to match catalogs and WCS.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
     Sachin Kumar Singh <sachinkumarsingh092@gmail.com>
Copyright (C) 2017-2025 Free Software Foundation, Inc.

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
#ifndef __GAL_MATCH_H__
#define __GAL_MATCH_H__

/* Include other headers if necessary here. Note that other header files
   must be included before the C++ preparations below */
#include <gnuastro/data.h>


/* C++ Preparations */
#undef __BEGIN_C_DECLS
#undef __END_C_DECLS
#ifdef __cplusplus
# define __BEGIN_C_DECLS extern "C" {
# define __END_C_DECLS }
#else
# define __BEGIN_C_DECLS                /* empty */
# define __END_C_DECLS                  /* empty */
#endif
/* End of C++ preparations */





/* The names are inspired by SQL's joining features:
   https://www.w3schools.com/sql/sql_join.asp . We do not need two OUTER
   modes because the order of inputs is important here (always using the
   first one as reference):

   - Inner: The output will contain matched rows in both catalogs. For
     example when you have two catalogs of galaxies from different
     telescopes and want to find the same galaxies in both.

   - Outer: The output will have the same number of rows as the second
     input, and the nearest row of the first input (within the given
     aperture) will be allocated for it. This is generally useful in
     classification scenarios: the first catalog will contain the
     "coordinate" columns and a certain value for that "coordinate". A
     second catalog's rows then need to be given the value of the nearest
     "coordinate" in the first.

   - Full (outer): Matching and non-matching rows of both. For example, you
     want to build a combined catalog from two separate input catalogs
     where some items match and some don't. */
enum gal_match_arrange{
  GAL_MATCH_ARRANGE_INVALID,   /* ==0 by default. */
  GAL_MATCH_ARRANGE_FULL,
  GAL_MATCH_ARRANGE_INNER,
  GAL_MATCH_ARRANGE_OUTER,
  GAL_MATCH_ARRANGE_OUTERWITHINAPERTURE,
};





gal_data_t *
gal_match_sort_based(gal_data_t *coord1, gal_data_t *coord2,
                      double *aperture, int sorted_by_first,
                      int inplace, size_t minmapsize, int quietmmap,
                      size_t *nummatched);

gal_data_t *
gal_match_kdtree(gal_data_t *coord1, gal_data_t *coord2,
                 gal_data_t *coord1_kdtree, size_t kdtree_root,
                 uint8_t arrange, double *aperture, size_t numthreads,
                 size_t minmapsize, int quietmmap, size_t *nummatched);





__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_TABLE_H__ */
