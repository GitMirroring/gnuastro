/*********************************************************************
Measure the gradient of pixel values around each pixel.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Corresponding author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
     Sepideh Eskandarlou <sepideh.eskandarlou@gmail.com>
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

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>

#include <gnuastro/data.h>
#include <gnuastro/gradient.h>
#include <gnuastro/dimension.h>





/* Low-level function that will return the magnitude and direction
   depending on the arguments. But it is not the main interface. */
static gal_data_t *
gradient_lowlevel(gal_data_t *input, uint8_t m0d1, uint8_t allngb,
                  uint8_t both, const char *caller)
{
  gal_data_t *out, *indbl=NULL;
  size_t i, ndim=input->ndim, *dsize=input->dsize;
  double *v, *o1, *o2, dx, dy, px, nx, py, ny, dir=NAN;
  size_t w=dsize[1], *dinc=gal_dimension_increment(ndim, dsize);

  /* Basic sanity checks. */
  if(m0d1!=0 && m0d1!=1)
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at '%s' to "
          "address the problem. The value to 'm0d1' should either be "
          "0 or 1, but it is %u", __func__, PACKAGE_BUGREPORT, m0d1);
  if(ndim!=2)
    error(EXIT_FAILURE, 0, "%s: this function is currently designed "
          "only for 2D inputs but the input has %zu dimensions. "
          "Please get in touch with us at '%s' if you need this "
          "feature", caller, ndim, PACKAGE_BUGREPORT);

  /* If the input is not double-precision floating point, convert it before
     continuing and allocate the output. */
  indbl = ( input->type==GAL_TYPE_FLOAT64
            ? input
            : gal_data_copy_to_new_type(input, GAL_TYPE_FLOAT64) );
  out=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, ndim, input->dsize,
                     input->wcs, 0, input->minmapsize, input->quietmmap,
                     NULL, NULL, NULL);
  if(both)
    gal_list_data_add_alloc(&out, NULL, GAL_TYPE_FLOAT64, ndim, dsize,
                            input->wcs, 0, input->minmapsize,
                            input->quietmmap, NULL, NULL, NULL);

  /* Parse all the pixels. */
  o1=out->array;
  v=indbl->array;
  o2=out->next?out->next->array:NULL;
  for(i=0;i<input->size;++i)
    {
      /* Initialize the prev/next elements and fill them up based on the
         neighbors. The initialization is important for the edges of the
         image (where the prev/next may not exist and we do not want to use
         values from the previous pixel). */
      px=py=nx=ny=NAN;
      GAL_DIMENSION_NEIGHBOR_OP(i, ndim, dsize, (int)(allngb?ndim:1), dinc,
          {
            /* First check if this neighbor is a 4-connected one (that is
               always needed. */
            if     (nind==i+1) nx=v[nind];
            else if(nind==i-1) px=v[nind];
            else if(nind==i+w) ny=v[nind];
            else if(nind==i-w) py=v[nind];

            /* In case the user asked for all neighbors, check for them too
               (and take the mean with the 4-connected neighbors). Note
               that because of NaN pixels in the images, it can happen that
               'nx', 'ny', 'px' and 'py' are not initialized yet. */
            else if(allngb)
              {
                if(     nind==i-w-1) py=isnan(py)?v[nind]:(py+v[nind])/2;
                else if(nind==i+w+1) ny=isnan(ny)?v[nind]:(ny+v[nind])/2;
                else if(nind==i-w+1) nx=isnan(nx)?v[nind]:(nx+v[nind])/2;
                else if(nind==i+w-1) px=isnan(px)?v[nind]:(px+v[nind])/2;
                else
                  error(EXIT_FAILURE, 0, "%s: a bug! Please contact us "
                        "at '%s' to solve the problem. The value of "
                        "'nind' (%zu) is not a recognized 8-connected "
                        "neighbor of pixel %zu", __func__,
                        PACKAGE_BUGREPORT, nind, i);
              }

            /* Should not happen. */
            else
              error(EXIT_FAILURE, 0, "%s: a bug! Please contact us "
                    "at '%s' to solve the problem. The value of "
                    "'nind' (%zu) is not a recognized 4-connected "
                    "neighbor of pixel %zu", __func__,
                    PACKAGE_BUGREPORT, nind, i);
          });

      /* Calculate the deltas and measure the gradient's magnitude or
         direction. For the direction, libc's 'atan2' will return values
         from -180 to +180, but that is not easy for human interpretation,
         so we change the range of angels from 0 to 360 by adding 360 to
         the negative directions. */
      dx = nx - px;
      dy = ny - py;
      if(m0d1) dir = atan2(dy, dx)*180.0f/M_PI;
      if(both)
        {
          o1[i] = dir<0 ? dir+360.0f : dir;
          o2[i] = sqrt( dx*dx + dy*dy );
        }
      else
        o1[i] = m0d1 ? (dir<0 ? dir+360.0f : dir) : sqrt( dx*dx + dy*dy );
    }

  /* Clean up (if necessary) and return the output. */
  if(indbl!=input) gal_data_free(indbl);
  return out;
}





/* High-level function to get the gradient magnitdue. */
gal_data_t *
gal_gradient_magnitude(gal_data_t *input, uint8_t allngb)
{
  return gradient_lowlevel(input, 0, allngb, 0, __func__);
}





/* High-level function to get the gradient direction. */
gal_data_t *
gal_gradient_direction(gal_data_t *input, uint8_t allngb)
{
  return gradient_lowlevel(input, 1, allngb, 0, __func__);
}





/* High-level function to get the gradient direction and magnitude in one
   run. */
gal_data_t *
gal_gradient_direction_magnitude(gal_data_t *input, uint8_t allngb)
{
  return gradient_lowlevel(input, 1, allngb, 1, __func__);
}
