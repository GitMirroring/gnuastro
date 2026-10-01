/*********************************************************************
Arithmetic - Do arithmetic operations on images.
Arithmetic is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2017-2026 Free Software Foundation, Inc.

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
#ifndef AUTHORS_CITE_H
#define AUTHORS_CITE_H

/* When any specific citation is necessary, please add its BibTeX (from ADS
   hopefully) to this variable along with a title decribing what this
   paper/book does for the progarm in a short line. In the following line
   put a row of '-' with the same length and then put the BibTeX.

   This macro will be used in 'gal_options_print_citation' function of
   'lib/options.c' (from the top Gnuastro source code directory). */

#define PROGRAM_BIBTEX ""                                               \
  "Paper describing '*-maskfilled' or 'collapse-*clip-fill-*' "         \
  "operators\n"                                                         \
  "-----------------------------------------------------------"         \
  "---------\n"                                                         \
  "@ARTICLE{2026RNAAS..10..271A,\n"                                     \
  "       author = {{Akhlaghi}, Mohammad "                              \
                    "and {Vives-Arias}, H{\'e}ctor "                    \
                    "and {Renard}, Pablo "                              \
                    "and {V{\'a}zquez Rami{\'o}}, H{\'e}ctor "          \
                    "and {Infante-Sainz}, Ra{\'u}l},\n"                 \
  "        title = \"{Gnuastro: Removing Extended/Diffuse Outliers "    \
                   "while Coadding}\",\n"                               \
  "      journal = {Research Notes of the American Astronomical "       \
                   "Society},\n"                                        \
  "     keywords = {Astronomy data reduction, "                         \
                   "Astronomy image processing, "                       \
                   "Astronomy software, 1861, 2306, 1855, "             \
                   "Instrumentation and Methods for Astrophysics},\n"   \
  "         year = 2026,\n"                                             \
  "        month = sep,\n"                                              \
  "       volume = {10},\n"                                             \
  "       number = {9},\n"                                              \
  "          eid = {271},\n"                                            \
  "        pages = {271},\n"                                            \
  "          doi = {10.3847/2515-5172/aea78e},\n"                       \
  "archivePrefix = {arXiv},\n"                                          \
  "       eprint = {2609.15529},\n"                                     \
  " primaryClass = {astro-ph.IM},\n"                                    \
  "       adsurl = {https://ui.adsabs.harvard.edu/abs/2026RNAAS..10..271A},\n" \
  "      adsnote = {Provided by the SAO/NASA Astrophysics Data System}\n" \
  "}"




#define PROGRAM_AUTHORS "Mohammad Akhlaghi"

#endif
