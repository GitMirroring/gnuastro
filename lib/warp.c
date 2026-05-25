/*********************************************************************
Warp -- Warp pixels of one dataset to another pixel grid.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Corresponding author:
     Pedram Ashofteh-Ardakani <pedramardakani@pm.me>
Contributing author(s):
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Copyright (C) 2022-2026 Free Software Foundation, Inc.

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

#include <errno.h>
#include <error.h>
#include <stdio.h>

#include <gnuastro/wcs.h>
#include <gnuastro/tile.h>
#include <gnuastro/type.h>
#include <gnuastro/warp.h>
#include <gnuastro/blank.h>
#include <gnuastro/pointer.h>
#include <gnuastro/polygon.h>
#include <gnuastro/threads.h>
#include <gnuastro/dimension.h>

#include <gnuastro-internal/timing.h>





/* Macros. */
#define WARP_NEXT_ODD(D) \
  ( (size_t)( ceil((D)) ) + ( (size_t)( ceil((D)) )%2==0 ? 1 : 0 ) )


#define WARP_WCSALIGN_H(IND,ES,IS1)     \
  (size_t)( ( (IND)%(IS1) ) * ( (ES)+1 ) \
            + ( (IND)/(IS1) ) * (1+(IS1)*( (ES)+1 )) )


#define WARP_WCSALIGN_V0(ES,IS0,IS1) \
  (size_t)(  1+(IS0)+(IS1)*( (IS0)+1 )*( (ES)+1) )


#define WARP_WCSALIGN_V(IND,ES,V0,IS1) \
  (size_t)( (V0)+(ES)*( (IND)+(IND)/(IS1) ) )





/* Per-thread parameters.*/
typedef struct warp_params_on_thread_t
{
  gal_warp_wcsalign_t *wa;  /* Pointer to main structure                 */

  gal_data_t        *tile;  /* Pointer to tile to be processed           */
  size_t     *start_coord;  /* start coord of tile                       */

  struct wcsprm     *iwcs;  /* Thread-safe copy of input wcs.            */
  struct wcsprm     *owcs;  /* Thread-safe copy of output wcs.           */

  /* Internal variables (allocated and freed internally)  */
  size_t               v0;  /* The first vertical corner in each pixel.  */
  size_t             gcrn;  /* Gap between corners of each row.          */
  gal_data_t    *vertices;  /* Stores all vertice coords of tile.        */
  double            *ocrn;  /* Temporary coords of pixel                 */
} warp_params_on_thread_t;





/* Generate the points on the outer boundary of a dsize[0] x dsize[1]
   matrix and return the array. */
static gal_data_t *
warp_alloc_perimeter(gal_data_t *input)
{
  /* Low level variables. */
  size_t ind, i;
  gal_data_t *pcrn;
  double *x=NULL, *y=NULL;
  size_t is0=input->dsize[0];
  size_t is1=input->dsize[1];
  int quietmmap=input->quietmmap;
  size_t minmapsize=input->minmapsize;

  /* High level variables. */
  size_t npcrn=2*(is0+is1);
  pcrn=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &npcrn, NULL, 0,
                      minmapsize, quietmmap, NULL, NULL, NULL);
  pcrn->next=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &npcrn, NULL, 0,
                            minmapsize, quietmmap, NULL, NULL, NULL);

  /* Find outermost pixel coordinates of the input image. Cover two
     corners at once to shorten the loop. */
  x=pcrn->array; y=pcrn->next->array;

  /* Top and bottom. */
  ind=0;
  for(i=is1+1; i--;)
    {
      x[ind]=i+0.5f;     y[ind]=0.5f; ind++;
      x[ind]=i+0.5f; y[ind]=is0+0.5f; ind++;
    }

  /* Left and right. */
  for(i=is0-1; i--;)
    {
      x[ind]=0.5f;     y[ind]=1.5f+i; ind++;
      x[ind]=0.5f+is1; y[ind]=1.5f+i; ind++;
    }

  /* Sanity check: let's make sure we have correctly covered the input
     perimeter. */
  if(ind!=npcrn)
    error(EXIT_FAILURE, 0, "%s: the input img perimeter of size <%zu> "
          "is not covered correctly. Currently on ind <%zu>.", __func__,
          npcrn, ind);

  return pcrn;
}





/* Create a base image with WCS consisting of the basic geometry
   keywords. */
static void
warp_wcsalign_init_output_from_params(gal_warp_wcsalign_t *wa)
{
  /* Low level variables. */
  size_t i, nkcoords, *osize;
  gal_data_t *input=wa->input;
  int quietmmap=input->quietmmap;
  size_t minmapsize=input->minmapsize;
  struct wcsprm *bwcs=NULL, *rwcs=NULL;
  gal_data_t *kcoords=NULL, *pcrn=NULL, *converted=NULL;
  double ocrpix[2], *xkcoords, *ykcoords, *x=NULL, *y=NULL;
  double pmin[2]={DBL_MAX, DBL_MAX}, pmax[2]={-DBL_MAX, -DBL_MAX}, tmp;

  /* Base WCS default parameters. */
  double pc[4]={-1, 0, 0, 1};
  double rcrpix[2]={1.0, 1.0};
  char **ctype=wa->ctype->array;
  struct wcsprm *iwcs=input->wcs;
  double *cdelt=wa->cdelt->array;
  double *center=wa->center->array;

  /* High level variables. */
  char  *cunit[2]={iwcs->cunit[0], iwcs->cunit[1]};

  /* Determine the output image size. */
  size_t iminr=GAL_BLANK_SIZE_T, imaxr=GAL_BLANK_SIZE_T; /* indexs of  */
  size_t imind=GAL_BLANK_SIZE_T, imaxd=GAL_BLANK_SIZE_T; /* extreme-um */
  double pminr=DBL_MAX, pmind=DBL_MAX, pmaxr=-DBL_MAX, pmaxd=-DBL_MAX;

  /* Create the reference WCS. */
  rwcs=gal_wcs_create(rcrpix, center, cdelt, pc, cunit, ctype, 2,
                      GAL_WCS_LINEAR_MATRIX_PC);

  /* Calculate the outer boundary of the input. */
  pcrn=warp_alloc_perimeter(input);
  converted=gal_wcs_img_to_world(pcrn, iwcs, 0);

  /* Get the minimum/maximum of the outer boundary. */
  x=converted->array; y=converted->next->array;
  for(i=converted->size; i--;)
    {
      if(x[i]<pminr) { pminr=x[i]; iminr=i; }
      if(y[i]<pmind) { pmind=y[i]; imind=i; }
      if(x[i]>pmaxr) { pmaxr=x[i]; imaxr=i; }
      if(y[i]>pmaxd) { pmaxd=y[i]; imaxd=i; }
    }

  /* Prepare the key world coorinates and change to image coordinates
     later. We are doing this to determine the CRPIX and NAXISi size for
     the final image. */
  nkcoords=5;
  kcoords=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &nkcoords, NULL, 0,
                         minmapsize, quietmmap, NULL, NULL, NULL);
  kcoords->next=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &nkcoords, NULL,
                               0, minmapsize, quietmmap, NULL, NULL, NULL);
  xkcoords=kcoords->array; ykcoords=kcoords->next->array;
  xkcoords[0]=x[iminr];   ykcoords[0]=y[iminr];   /* min RA       */
  xkcoords[1]=x[imaxr];   ykcoords[1]=y[imaxr];   /* max RA       */
  xkcoords[2]=x[imind];   ykcoords[2]=y[imind];   /* min Dec      */
  xkcoords[3]=x[imaxd];   ykcoords[3]=y[imaxd];   /* max Dec      */
  xkcoords[4]=center[0];  ykcoords[4]=center[1];  /* Image center */

  /* Convert to pixel coords. */
  gal_wcs_world_to_img(kcoords, rwcs, 1);

  /* Determine output image size. */
  if( wa->widthinpix )
    osize=wa->widthinpix->array;
  else
    {
      /* Automatic: the first four coordinates are the extreme-um RA/Dec. */
      for(i=4; i--;)
        {
          pmin[0] = xkcoords[i] < pmin[0] ? xkcoords[i] : pmin[0];
          pmin[1] = ykcoords[i] < pmin[1] ? ykcoords[i] : pmin[1];
          pmax[0] = xkcoords[i] > pmax[0] ? xkcoords[i] : pmax[0];
          pmax[1] = ykcoords[i] > pmax[1] ? ykcoords[i] : pmax[1];
        }

      /* Size must be odd so the image would have a center value. Also, the
         indices are swapped since number of columns defines the horizontal
         part of the center and vice versa. To calculate the output image
         size, measure the difference between center and outermost edges of
         the input image (in pixels). Since this is the distance from
         center to the furthest edge of the image, the value must be
         multiplied by two. */
      osize=gal_pointer_allocate(GAL_TYPE_SIZE_T, 2, 0,
                                 __func__, "osize");
      tmp=2*fmax( fabs(ykcoords[4]-pmin[1]),
                  fabs(ykcoords[4]-pmax[1]) );
      osize[0] = WARP_NEXT_ODD(tmp);
      tmp=2*fmax( fabs(xkcoords[4]-pmin[0]),
                  fabs(xkcoords[4]-pmax[0]) );
      osize[1] = WARP_NEXT_ODD(tmp);
    }

  /* Set the CRPIX value

     Note: os1 is number of columns, so we use it to define CRPIX in the
     horizontal axis, and vice versa. */
  ocrpix[0]= 1.5f + osize[1]/2.0f - xkcoords[4];
  ocrpix[1]= 1.5f + osize[0]/2.0f - ykcoords[4];

  /* Create the base WCS. */
  bwcs=gal_wcs_create(ocrpix, center, cdelt, pc, cunit,
                      ctype, 2, GAL_WCS_LINEAR_MATRIX_PC);

  /* Make sure that the size is reasonable (i.e., less than 100000 pixels
     on a side). This can happen when a wrong central coordinate is
     requested. */
  if(osize[0]>100000 || osize[1]>100000)
    error(EXIT_SUCCESS, 0, "%s: the output image size (%zu x %zu pixels) "
          "is unreasonably large. This may be due to a mistake in the "
          "given central coordinate compared to the input image (the "
          "given center is too far from the image)", __func__, osize[1],
          osize[0]);

  /* Create the output image dataset with the base WCS. */
  wa->output=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 2, osize, bwcs, 0,
                            minmapsize, quietmmap,
                            GAL_WARP_OUTPUT_NAME_WARPED, NULL, NULL);

  /* Clean up. */
  wcsfree(bwcs);
  wcsfree(rwcs);
  gal_list_data_free(pcrn);
  gal_list_data_free(kcoords);
  gal_list_data_free(converted);
  if(wa->widthinpix==NULL) free(osize);

  /* Must be freed after wcsfree() is called! */
  free(bwcs);
  free(rwcs);
}





/* Initialize the vertices of each output pixel (accounting for any
   edge-sampling). */
static void
warp_wcsalign_init_vertices(warp_params_on_thread_t *wt)
{
  size_t es=wt->wa->edgesampling;

  size_t ind, i, j, ix, iy;
  double gap=1.0f/(es+1.0f);
  double *x=NULL, *y=NULL, x0, y0, row, col;
  uint8_t quietmmap=wt->wa->input->quietmmap;
  size_t minmapsize=wt->wa->input->minmapsize;
  size_t os0= wt->tile ? wt->tile->dsize[0] : wt->wa->output->dsize[0];
  size_t os1= wt->tile ? wt->tile->dsize[1] : wt->wa->output->dsize[1];

  size_t nvcrn=es*os0;
  size_t nhcrn=es*os1+os1+1;
  size_t v0=WARP_WCSALIGN_V0(es,os0,os1);
  size_t nvertices=nvcrn*(os1+1)+nhcrn*(os0+1);

  /* Allocate the space for all the vertice coordinates. */
  wt->vertices=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &nvertices,
                               NULL, 0, minmapsize, quietmmap, NULL,
                               NULL, NULL);
  wt->vertices->next=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &nvertices,
                                    NULL, 0, minmapsize, quietmmap,
                                    NULL, NULL, NULL);

  /* Parse all pixels. */
  x=wt->vertices->array;
  y=wt->vertices->next->array;

  /* Account for the tile shift */
  x0=wt->start_coord[1];
  y0=wt->start_coord[0];

  for(ind=os0*os1; ind--;)
    {
      /* For easy reading. */
      row=ind%os1;
      col=floor(ind/os1);
      ix=WARP_WCSALIGN_H(ind,es,os1);
      iy=WARP_WCSALIGN_V(ind,es,v0,os1);

      /* Bottom left edge of pixel. */
      x[ix]= x0+0.5f+row;
      y[ix]= y0+0.5f+col;

      for(i=es; i--;)
        {
          /* Horizontal */
          j=ix+i+1;
          x[j]=x0+0.5f+row+gap+i*gap;
          y[j]=y0+0.5f+col;

          /* Vertical */
          j=iy+i;
          x[j]=x0+0.5f+row;
          y[j]=y0+0.5f+col+gap+i*gap;
        }
    }

  /* Top */
  for(i=nhcrn; i--;)
    {
      j=v0-nhcrn+i;
      x[j]=x0+0.5f+gap*i;
      y[j]=y0+0.5f+os0;
    }

  /* Right */
  for(ind=os1-1; ind<os1*os0; ind+=os1)
    {
      /* For easy reading. */
      row=(size_t)(ind%os1);
      col=(size_t)(ind/os1);
      ix=WARP_WCSALIGN_H(ind,es,os1);
      iy=WARP_WCSALIGN_V(ind,es,v0,os1);

      /* Bottom right */
      j=ix+es+1;
      x[j]=x0+0.5f+os1;
      y[j]=y0+0.5f+col;

      /* Right vertice */
      for(i=es; i--;)
        {
          j=iy+es+i;
          x[j]=x0+0.5f+os1;
          y[j]=y0+0.5f+col+gap+i*gap;
        }
    }

  /* For a check:
  for(i=0;i<nvertices;++i)
    printf("%zu: %f, %f\n", i, x[i], y[i]);
  printf("%s: GOOD\n", __func__); exit(0);
  */
}





static void
warp_pixel_perimeter_ccw(warp_params_on_thread_t *wt, size_t ind)
{
  /* Low-level variables. */
  size_t i, j, hor, ver, ic;
  double *xcrn=NULL, *ycrn=NULL, *ocrn=wt->ocrn;

  /* High-level variables. */
  size_t v0=wt->v0;
  size_t gcrn=wt->gcrn;
  size_t es=wt->wa->edgesampling;
  size_t os1= wt->tile ? wt->tile->dsize[1] : wt->wa->output->dsize[1];

  /* Set ocrn, the corners of each output pixel. */
  xcrn=wt->vertices->array;
  ycrn=wt->vertices->next->array;

  /* Index of surrounding vertices for this pixel. */
  hor=WARP_WCSALIGN_H(ind, es, os1);
  ver=WARP_WCSALIGN_V(ind, es, v0, os1);

  /* All four corners

     WARNING: this block of code highly depends on the ordering, take
     extra care while refactoring

     ocrn: all output edges transformed into the input image pixel
     coordiantes

     ic: index (position in array) of the current pixel edges

     io: index of the output (reference) pixel edges


       left edge -> +---------+ <- top edge
                    |         |
                    |         |
                    |         |
     bottom edge -> +---------+ <- right edge */


  /* bottom left */
  i=0;
  j=hor;
  ocrn[2*i]=xcrn[j]; ocrn[2*i+1]=ycrn[j];

  /* bottom right */
  i=es+1;
  j=hor+es+1;
  ocrn[2*i]=xcrn[j]; ocrn[2*i+1]=ycrn[j];

  /* top right */
  i=2*(es+1);
  j=hor+es+1+gcrn;
  ocrn[2*i]=xcrn[j]; ocrn[2*i+1]=ycrn[j];

  /* top left */
  i=3*(es+1);
  j=hor+gcrn;
  ocrn[2*i]=xcrn[j]; ocrn[2*i+1]=ycrn[j];

  /* Sampling corners of the output pixel on the input image. */
  for(i=es; i--;)
    {
      /* Bottom vertice: 0*(es+1)+(i+1) */
      ic=i+1;
      j=hor+i+1;
      ocrn[2*ic]=xcrn[j]; ocrn[2*ic+1]=ycrn[j];

      /* Right vertice:  1*(es+1)+(i+1) */
      ic=i+2+es;
      j=ver+es+i;
      ocrn[2*ic]=xcrn[j]; ocrn[2*ic+1]=ycrn[j];

      /* Top vertice:    2*(es+1)+(i+1) */
      ic=i+3+2*es;
      j=hor+es+gcrn-i;
      ocrn[2*ic]=xcrn[j]; ocrn[2*ic+1]=ycrn[j];

      /* Left vertice:   3*(es+1)+(i+1) */
      ic=i+4+3*es;
      j=ver+es-i-1;
      ocrn[2*ic]=xcrn[j]; ocrn[2*ic+1]=ycrn[j];
    }
}





static void
warp_wcsalign_check_2d(gal_data_t *in, uint8_t type, const char *func,
                        char *name, char *comment)
{
  if(!in)
    error(EXIT_FAILURE, 0, "%s: no '%s' specified. %s", func,
          name, comment);
  if(in->size!=2)
    error(EXIT_FAILURE, 0, "%s: '%s' takes exactly 2 values, currently "
          "detected %zu values", func, name, in->size);
  if(in->type!=type)
    error(EXIT_FAILURE, 0, "%s: '%s' must have a type of '%s' but has "
          "type '%s'", func, name, gal_type_name(type, 1),
          gal_type_name(in->type, 1));
}





/* Create the output image using the WCS struct from the given 'gridfile'
   and 'gridhdu'. */
static void
warp_wcsalign_init_output_from_wcs(gal_warp_wcsalign_t *wa,
                                   const char *func)
{
  gal_data_t *output=NULL;
  int quietmmap=wa->input->quietmmap;
  size_t *dsize=wa->widthinpix->array, minmapsize=wa->input->minmapsize;

  /* Create the output image dataset with the target WCS given. */
  output=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 2, dsize, wa->twcs, 0,
                        minmapsize, quietmmap, GAL_WARP_OUTPUT_NAME_WARPED,
                        NULL, NULL);

  /* Write to wcsalign data type for later use. */
  wa->output=output;
}





/* Check parameters for aligning and pixelarea. */
static void
warp_check_basic_params(gal_warp_wcsalign_t *wa, const char *func)
{
  size_t es=wa->edgesampling;

  /* Check if input and 'wa' are not NULL! */
  if(wa==NULL) error(EXIT_FAILURE, 0, "%s: 'wa' structure is NULL", func);
  if(wa->input==NULL) error(EXIT_FAILURE, 0, "%s: input is NULL", func);

  /* This function assumes the input is double precision. */
  if(wa->input->type != GAL_TYPE_FLOAT64)
    error(EXIT_FAILURE, 0, "%s: input must have a double precision "
          "floating point type, but its type is '%s', you can use "
          "'gal_data_copy_to_new_type' or "
          "'gal_data_copy_to_new_type_free' for the conversion", func,
          gal_type_name(wa->input->type, 1));

  /* Check 'edgesampling', can't compare to '0' since it has meaning, can't
     check if negative since it is an unsigned type. */
  if(wa->edgesampling==GAL_BLANK_SIZE_T)
    error(EXIT_FAILURE, 0, "%s: no 'edgesampling' specified. This is the "
          "Order of samplings along each pixel edge", func);
  if(wa->edgesampling>999)
    error(EXIT_FAILURE, 0, "%s: edgesampling takes zero OR a positive "
          "integer value of type 'size_t', <%zu> is too big which might "
          "be a bad cast", func, wa->edgesampling);

  /* If 'numthreads' is 0, use the number of threads available to the
     system. */
  if(wa->numthreads==GAL_BLANK_SIZE_T || wa->numthreads==0)
    wa->numthreads=gal_threads_number();

  /* Initialize the internal parameters. */
  wa->ncrn=4*es+4;
  wa->nhor=GAL_BLANK_SIZE_T;
  wa->gcrn=GAL_BLANK_SIZE_T;
}





static void
warp_wcsalign_init_internals(warp_params_on_thread_t *wt)
{
  size_t es  = wt->wa->edgesampling;
  size_t os0 = wt->tile ? wt->tile->dsize[0] : wt->wa->output->dsize[0];
  size_t os1 = wt->tile ? wt->tile->dsize[1] : wt->wa->output->dsize[1];

  /* Common variables to simplify next functions. */
  wt->gcrn=1+os1*(es+1);
  wt->v0=WARP_WCSALIGN_V0(es, os0, os1);
}





static void *
warp_wcsalign_init_params(gal_warp_wcsalign_t *wa, const char *func)
{
  gal_data_t *output=NULL;
  int quietmmap=wa->input->quietmmap;
  size_t minmapsize=wa->input->minmapsize, *dsize=NULL, *tmp=NULL;

  warp_check_basic_params(wa, func);

  /* Check 'coveredfrac'. */
  if(isnan( wa->coveredfrac ))
    error(EXIT_FAILURE, 0, "%s: no 'coveredfrac' specified. This is the "
          "acceptable fraction of output covered", func);
  if(wa->coveredfrac<0.0f || wa->coveredfrac>1.0f)
    error(EXIT_FAILURE, 0, "%s: coveredfrac takes exactly on positive "
          "value less than or equal to 1.0, but it is given a value "
          "of %f", func, wa->coveredfrac);

  /* If a target WCS is given ignore other variables and initialize the
     output image. */
  if(wa->twcs)
    {
      /* Check the widthinpix element. */
      warp_wcsalign_check_2d(wa->widthinpix, GAL_TYPE_SIZE_T, func,
                             "widthinpix", "This is the output image "
                             "size in pixels");
      warp_wcsalign_init_output_from_wcs(wa, func);

      /* Warp will ignore the following parameters, warn the user if
         detected any. */
      if(wa->cdelt || wa->center || wa->ctype)
        error(EXIT_SUCCESS, 0, "%s: WARNING: target WCS is already "
              "defined with 'gridfile' and 'gridhdu', ignoring extra "
              "non-linear parameter(s) given", func);
      return NULL;
    }

  /* No 'gridfile' given, Warp must create the output WCS using given
     parameters. Proceed with checking the 2D input parameters. */
  warp_wcsalign_check_2d(wa->ctype, GAL_TYPE_STRING, func, "ctype",
                         "Any pair of valid WCSLIB ctype is "
                         "allowed, e.g. 'RA---TAN, DEC--TAN'");
  warp_wcsalign_check_2d(wa->cdelt, GAL_TYPE_FLOAT64, func, "cdelt",
                         "This is the pixel scale in degrees");
  warp_wcsalign_check_2d(wa->center, GAL_TYPE_FLOAT64, func, "center",
                         "This is the output image center in degrees");

  /* Check 'widthinpix', it can be null for automatic detection. */
  if(wa->widthinpix)
    {
      warp_wcsalign_check_2d(wa->widthinpix, GAL_TYPE_SIZE_T, func,
                             "widthinpix", "This is the output "
                             "image size");
      tmp=wa->widthinpix->array;
      if(tmp[0]%2==0 || tmp[1]%2==0)
        error(EXIT_FAILURE, 0, "%s: 'widthinpix' takes exactly 2 ODD "
              "values, detected EVEN value in %zux%zu", func, tmp[0],
              tmp[1]);
    }

  /* Initialize the output image for further processing. */
  warp_wcsalign_init_output_from_params(wa);

  /* Create the check maximum fraction covered dataset if asked to, then
     return. */
  output=wa->output;
  dsize=output->dsize;
  if(wa->checkmaxfrac)
    output->next=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 2, dsize,
                                wa->twcs, 0, minmapsize, quietmmap,
                                GAL_WARP_OUTPUT_NAME_MAXFRAC,
                                NULL, NULL);
  return NULL;
}





/* Convert the necessary vertice coordinates. */
static void *
warp_wcsalign_init_convert(warp_params_on_thread_t *wt)
{
  gal_data_t *vertices=NULL;
  size_t size=wt->vertices->size;
  double *xarr=wt->vertices->array;
  int quietmmap=wt->vertices->quietmmap;
  double *yarr=wt->vertices->next->array;
  size_t minmapsize=wt->vertices->minmapsize;

  /* Allocate the non-allocated vertices table for this thread. */
  gal_list_data_add_alloc(&vertices, xarr, GAL_TYPE_FLOAT64,
                          1, &size, NULL, 0, minmapsize, quietmmap,
                          NULL, NULL, NULL);
  gal_list_data_add_alloc(&vertices, yarr, GAL_TYPE_FLOAT64,
                          1, &size, NULL, 0, minmapsize, quietmmap,
                          NULL, NULL, NULL);
  gal_list_data_reverse(&vertices); /* '_add' is last-in-first-out. */

  /* Convert the coordinates. */
  gal_wcs_img_to_world(vertices, wt->owcs, 1);
  gal_wcs_world_to_img(vertices, wt->iwcs, 1);

  /* Clean up: since the 'array' pointer is within a larger allocated
     array, we shouldn't free it when freeing the table, so we'll set it to
     NULL. */
  vertices->array=vertices->next->array=NULL;
  gal_list_data_free(vertices);

  return NULL;
}





/* Determine the final image size and allocate the output array
   accordingly.

   'is0' indicates number of rows available in the input fits image,
   while 'is1' indicates nummber of columns. The same goes for 'os0'
   and 'os1'. See the following figure:

                       +------------------------+
                    /  |                        |
                    |  |             N          |
                    |  |             ^          |
                    |  |             |          |
               is0 <   |             |          |
                    |  |     E <-----+          |
                    |  |                        |
                    |  |      input image       |
                    \  |                        |
                       +------------------------+
                        \__________ ___________/
                                   v
                                  is1

   NOTE: please keep in mind that dsize[1] is NAXIS1 and dsize[0] is
   NAXIS2 in FITS format.

   In preparations, 'pcrn' is of type 'gal_data_t' linked lists. The
   first list holds RA coords, and the next is Dec coords. This
   variable is filled with the outer-most pixel coordinates. The
   purpose of this variable is to hold the min and max RA and Dec
   coordinates 'temporarily', so we can determine the output image size
   later.

   nhcrn and nvcrn:

   'nhcrn' stands for number of horizontal corners, and similarly
   'nvcrn' stands for number of vertical corners. Note that number of
   'corners' is different from number of 'pixels'. Each pixel has many
   corners.

   Also bear in mind, to keep from counting repeated corners on the
   image edges, we let the horizontal corners devour the first and
   last vertical corners. See the following figure:

                       hc6   hc7   hc8   hc9   hc10
                       +-----+-----+-----+-----+
                    /  |                       |
                    |  |       img: 4x3        |
                    |  x vc2                   x vc4
            is0=3   |  |                  N    |
                    |  |                  ^    |
                    |  x vc1              |    x vc3
                    |  |            E <---+    |
                    \  |                       |
                       +-----+-----+-----+-----+
                       hc1   hc2   hc3   hc4   hc5

                        \_____________________/
                                 is1=4

   In the example above, for an image of 4x3, there will be 5
   horizontal and 2 vertical corners in each axis, hence we would have
   10 horizontal and 4 vertical corners: 'nhcrn=2*(is1+1)' and
   'nvcrn=2*(is0-1)'. The total number of corners will be:
   'nhcrn+nvcrn=2*(is0-1)+2*(is1+1)=2*(is0+is1)=2*(4+3)=14'.


   After finding out the min and max RA and Decs, the 'pcrn' is
   projected back to pixel coordinates.
*/
void
warp_wcsalign_init(warp_params_on_thread_t *wt)
{
  /* Set up the output image corners in pixel coords. */
  warp_wcsalign_init_vertices(wt);

  /* Project the output image corners to the input image pixel coords. We
     only want one job per thread, so the number of jobs and the number of
     threads are the same. */
  warp_wcsalign_init_convert(wt);

  /* Now that the output image is ready, initialize the helper internal
     variables for future processing. */
  warp_wcsalign_init_internals(wt);
}





static void
warp_wcsalign_onpix(warp_params_on_thread_t *wt, size_t ind, size_t oind)
{
  size_t ic, temp, numinput=0;
  double xmin, xmax, ymin, ymax;
  gal_data_t *input=wt->wa->input;
  gal_data_t *output=wt->wa->output;
  long xstart, ystart, xend, yend, x, y; /* Might be negative */
  double filledarea, v, *ocrn=wt->ocrn, pcrn[8], opixarea;

  size_t numcrn=0;
  size_t ncrn=wt->wa->ncrn;
  size_t is0=input->dsize[0];
  size_t is1=input->dsize[1];
  double *inputarr=input->array;
  double *outputarr=output->array;
  double ccrn[2*GAL_POLYGON_MAX_CORNERS], area;
  double *maxfrac=output->next ? output->next->array : NULL;

  /* Initialize if asked for each pixel's maximum coverage fraction. */
  if(maxfrac) maxfrac[oind]=-DBL_MAX;

  /* Initialize the output pixel value: */
  outputarr[oind] = filledarea = 0.0f;

  /* get coordinates of pixel's vertices */
  warp_pixel_perimeter_ccw(wt, ind);

  /* Find overlapping pixels. */
  xmin =  DBL_MAX; ymin =  DBL_MAX;
  xmax = -DBL_MAX; ymax = -DBL_MAX;
  for(ic=ncrn; ic--;)
    {
      temp=ic*2;
      if(xmin > ocrn[ temp   ]) { xmin = ocrn[ temp   ]; }
      if(xmax < ocrn[ temp   ]) { xmax = ocrn[ temp   ]; }
      if(ymin > ocrn[ temp+1 ]) { ymin = ocrn[ temp+1 ]; }
      if(ymax < ocrn[ temp+1 ]) { ymax = ocrn[ temp+1 ]; }
    }

  /* Start and end in both dimensions. */
  xstart = GAL_DIMENSION_NEARESTINT_HALFHIGHER( xmin );
  ystart = GAL_DIMENSION_NEARESTINT_HALFHIGHER( ymin );
  xend   = GAL_DIMENSION_NEARESTINT_HALFLOWER(  xmax ) + 1;
  yend   = GAL_DIMENSION_NEARESTINT_HALFLOWER(  ymax ) + 1;

  /* Check which input pixels we are covering. */
  for(y=ystart;y<yend;++y)
    {
      /* If the pixel isn't in the image (note that the pixel
         coordinates start from 1), skip this pixel. */
      if( y<1 || y>is0 ) continue;

      /* Y of base pixel vertices, in pixel coords. */
      pcrn[1]=y-0.5f; pcrn[3]=y-0.5f;
      pcrn[5]=y+0.5f; pcrn[7]=y+0.5f;

      for(x=xstart;x<xend;++x)
        {
          if( x<1 || x>is1 ) continue;

          /* X of base pixel vertices, in pixel coords. */
          pcrn[0]=x-0.5f; pcrn[2]=x+0.5f;
          pcrn[4]=x+0.5f; pcrn[6]=x-0.5f;

          /* Read the value of the input pixel. */
          v=inputarr[(y-1)*is1+x-1];

          /* Find the overlapping (clipped) polygon and its area.

             In theory, instead of 'gal_polygon_area_flat', we should be
             using the spherical polygon area ('gal_polygon_area_flat').
             But the area that is calculate here is not used in an absolute
             sense: it is only relative for comparison with other input
             pixels that overlap with this output pixel. Therefore, because
             the flat area calculation is faster, we'll suffice to that
             unless we discover there is any problem with it. */
          numcrn=0; /* initialize it. */
          gal_polygon_clip(ocrn, ncrn, pcrn, 4, ccrn, &numcrn);
          area=gal_polygon_area_flat(ccrn, numcrn);

          /* Write each pixel's maximum coverage fraction if asked. */
          if( maxfrac ) maxfrac[oind] = fmax(area, maxfrac[oind]);

          /* Add the fractional value of this pixel. If this output
             pixel covers a NaN pixel in the input grid, then
             calculate the area of this NaN pixel to account for it
             later. */
          if( !isnan(v) )
            {
              numinput+=1;
              filledarea+=area;
              outputarr[oind]+=v*area;

              /* Check
                 printf("Check: numinput %zu filledarea %f "
                 "outputarr[%zu]=%f\n",
                 filledarea, ind, outputarr[ind]);
              */
            }
        }
    }

  /* Replace untouched pixels with NAN in the 'maxfrac' array. */
  if( maxfrac && maxfrac[oind]==-DBL_MAX ) maxfrac[oind]=NAN;

  /* See if the pixel value should be set to NaN or not (because of not
     enough coverage). For a description of why we are not using
     'gal_polygon_area_sky', see the comment above the previous call to
     'gal_polygon_area_flat' above. */
  opixarea=gal_polygon_area_flat(ocrn, ncrn);
  if( numinput && filledarea/opixarea < wt->wa->coveredfrac-1e-5)
    numinput=0;

  /* Write the final value and return. */
  if( numinput==0 ) outputarr[oind]=NAN;
}





/* Clean up the internally allocated variables from the 'wt' struct. */
static void
warp_wcsalign_free(warp_params_on_thread_t *wt)
{
  free(wt->ocrn);
  free(wt->start_coord);
  if(wt->iwcs) gal_wcs_free(wt->iwcs);
  if(wt->owcs) gal_wcs_free(wt->owcs);
  if(wt->vertices) gal_list_data_free(wt->vertices);

  wt->vertices=NULL;
  wt->start_coord=NULL;
}





static void *
warp_wcsalign_onthread(void *inparam)
{
  struct gal_threads_params *tprm=(struct gal_threads_params *)inparam;
  gal_warp_wcsalign_t *wa=(gal_warp_wcsalign_t *)tprm->params;

  size_t *coord;
  size_t ncrn=wa->ncrn;
  size_t i, ind, oind, tile;
  warp_params_on_thread_t wt;

  /* Basic initializations. */
  wt.wa=wa;
  coord=malloc(wa->tiles[0].ndim * sizeof(*coord));
  wt.start_coord=gal_pointer_allocate(GAL_TYPE_SIZE_T, wa->tiles[0].ndim,
                                      0, __func__, "start_coord");
  wt.ocrn=gal_pointer_allocate(GAL_TYPE_FLOAT64, (2*ncrn), 0, __func__,
                               "ocrn");
  wt.iwcs=gal_wcs_copy(wa->input->wcs);
  wt.owcs=gal_wcs_copy(wa->output->wcs);

  /* Go over every tile. */
  for(i=0; tprm->indexs[i] != GAL_BLANK_SIZE_T; ++i)
    {
      tile=tprm->indexs[i];
      wt.tile=&wa->tiles[tile];

      /* Get the starting coordinate of the processed tile */
      gal_tile_start_coord(wt.tile, wt.start_coord);

      /* Calculate and allocate the output image size and WCS. */
      warp_wcsalign_init(&wt);

      /* Loop over pixels of the processed tile. */
      for(ind=0; ind<wt.tile->size; ++ind)
        {
          /* Get also the 'global index' of the current pixel */
          coord[0]=wt.start_coord[0]+ind/wt.tile->dsize[1];
          coord[1]=wt.start_coord[1]+ind%wt.tile->dsize[1];
          oind=gal_dimension_coord_to_index(wa->output->ndim,
                                            wa->output->dsize, coord);

          warp_wcsalign_onpix(&wt, ind, oind);
        }

      /* Vertices are initialized by 'warp_wcsalign_init'. We free after
         each processed tile, as tile size may vary a little and so the
         number of vertices. */
      gal_list_data_free(wt.vertices);
      wt.vertices=NULL;
    }

  /* Clean up the internally allocated variables. */
  free(coord);
  warp_wcsalign_free(&wt);

  /* Wait for all the other threads to finish, then return. */
  if(tprm->b) { pthread_barrier_wait(tprm->b); }
  return NULL;
}





/* Helper function that returns an empty set of the wcsalign data structure
   to prevent using uninitialized variables without warnings. Please note
   if you are not using this template to set 'gal_warp_wcsalign_t' values,
   you MUST pass NULL to unused pointers at least. */
gal_warp_wcsalign_t
gal_warp_wcsalign_init()
{
  gal_warp_wcsalign_t wa;

  /* Initialize pointers with NULL. */
  wa.twcs=NULL;
  wa.cdelt=NULL;
  wa.ctype=NULL;
  wa.input=NULL;
  wa.tiles=NULL;
  wa.center=NULL;
  wa.output=NULL;
  wa.widthinpix=NULL;

  /* Initialize values. */
  wa.checkmaxfrac=0;
  wa.gcrn=GAL_BLANK_SIZE_T;
  wa.ncrn=GAL_BLANK_SIZE_T;
  wa.nhor=GAL_BLANK_SIZE_T;
  wa.numthreads=GAL_BLANK_SIZE_T;
  wa.coveredfrac=GAL_BLANK_FLOAT64;
  wa.edgesampling=GAL_BLANK_SIZE_T;

  /* Retrun the initialized array. */
  return wa;
}





/* Finalize the output 'gal_data_t' image in 'wa->output'. */
void
gal_warp_wcsalign(gal_warp_wcsalign_t *wa)
{
  /* Fill the output image. */
  size_t i, ntiles=1, *ntilesperdim, *junk;

  /* Run a sanity check on the input parameters and initialize the output
     image. */
  warp_wcsalign_init_params(wa, __func__);
  ntilesperdim=gal_tile_full(wa->output, wa->tilesize, wa->remainderfrac,
                             &wa->tiles, 1, &junk);

  /* Get the final number of tiles and spin-off the threads. */
  for(i=0;i<wa->output->ndim;i++) ntiles*=ntilesperdim[i];
  gal_threads_spin_off(warp_wcsalign_onthread, wa,
                       ntiles, wa->numthreads,
                       wa->input->minmapsize, wa->input->quietmmap);

  /* Free the tile structure. */
  gal_data_array_free(wa->tiles, ntiles, 0);
}





static void
warp_pixelarea_init_convert(warp_params_on_thread_t *wt)
{
  /* Change the CRVALs to be on 0 and 180 (this helps to avoid issues with
     an image passing the RA=0.0, or high declination problems. */
  double crval[2]={180.0f,0.0f};

  /* Initialize the input/output WCSs and convert all the vertices. */
  wt->owcs=NULL;
  wt->iwcs=gal_wcs_copy_new_crval(wt->wa->input->wcs, crval);
  gal_wcs_img_to_world(wt->vertices, wt->iwcs, 1);

  /* For a check (after conversion).
  size_t j;
  double *ra=wa->vertices->array;
  double *dec=wa->vertices->next->array;
  for(j=0;j<wa->vertices->size;++j)
    printf("%-20.15f %-20.15f\n", ra[j], dec[j]);
  */

  /* Clean up and reset the input WCS to NULL. */
  gal_wcs_free(wt->iwcs);
  wt->iwcs=NULL;
}





static void
warp_pixelarea_init(warp_params_on_thread_t *wt)
{
  /* Initialize vertices coordinates */
  warp_wcsalign_init_vertices(wt);

  /*Convert coordinates via wcs */
  warp_pixelarea_init_convert(wt);

  /* Initialize relevant internal values */
  warp_wcsalign_init_internals(wt);
}





static void *
warp_pixelarea_onthread(void *inparam)
{
  /* Thread variables. */
  struct gal_threads_params *tprm=(struct gal_threads_params *)inparam;
  gal_warp_wcsalign_t *wa=(gal_warp_wcsalign_t *)tprm->params;

  /* Low-level variables. */
  size_t *coord;
  size_t ncrn=wa->ncrn;
  size_t i, ind, oind, tile;
  warp_params_on_thread_t wt;
  double area, *outputarr=wa->output->array;

  /* Necessary allocations. */
  wt.wa=wa;
  coord=malloc(wa->tiles[0].ndim * sizeof(*coord));
  wt.start_coord=gal_pointer_allocate(GAL_TYPE_SIZE_T, wa->input->ndim,
                                      0, __func__, "start_coord");
  wt.ocrn=gal_pointer_allocate(GAL_TYPE_FLOAT64, (2*ncrn), 0, __func__,
                               "ocrn");

  /* Loop over pixels given from the 'warp' function. */
  for(i=0; tprm->indexs[i] != GAL_BLANK_SIZE_T; ++i)
    {
      tile=tprm->indexs[i];
      wt.tile=&wa->tiles[tile];

      /* Get the starting coordinate of the processed tile */
      gal_tile_start_coord(wt.tile, wt.start_coord);

      /* Calculate and allocate the output image size and WCS. */
      warp_pixelarea_init(&wt);

      /* Loop over pixels given from the 'warp' function. */
      for(ind=0; ind<wt.tile->size; ++ind)
        {
          /* Get also the 'global index' of the current pixel */
          coord[0]=wt.start_coord[0]+ind/wt.tile->dsize[1];
          coord[1]=wt.start_coord[1]+ind%wt.tile->dsize[1];
          oind=gal_dimension_coord_to_index(wa->output->ndim,
                                            wa->output->dsize, coord);

          /* Get the coordinates of the pixel vertices (wt.ocrn) */
          warp_pixel_perimeter_ccw(&wt, ind);

          /* Calculate the area. */
          area=gal_polygon_area_sky(wt.ocrn, wa->ncrn);
          outputarr[oind]=area;
        }

        /* Vertices are initialized by warp_wcsalign_init().
           We free after each processed tile, as tile size may vary a little
           and so the number of vertices. */
        gal_list_data_free(wt.vertices);
        wt.vertices=NULL;
    }

  /* Clean up the internally allocated variables. */
  free(coord);
  warp_wcsalign_free(&wt);

  /* Wait for all the other threads to finish, then return. */
  if(tprm->b) { pthread_barrier_wait(tprm->b); }
  return NULL;
}





/* Calculate input pixel area covering the sky and write it to output. This
   function has a debugging nature and is not used through the aligning
   process. */
void
gal_warp_pixelarea(gal_warp_wcsalign_t *wa)
{
  gal_data_t *input=wa->input;
  size_t i, ntiles=1, *ntilesperdim, *junk;

  /* Basic sanity check. */
  warp_check_basic_params(wa, __func__);

  /* Create the output dataset. */
  wa->output=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, input->ndim,
                            input->dsize, input->wcs, 0,
                            input->minmapsize, input->quietmmap,
                            "PIX-AREA", NULL, NULL);

  /* Create the tile structure and set the number of tiles. */
  ntilesperdim=gal_tile_full(wa->output, wa->tilesize,
                             wa->remainderfrac, &wa->tiles, 1, &junk);
  for(i=0;i<wa->output->ndim;i++) ntiles*=ntilesperdim[i];

  /* Calculate pixel area on WCS and write to output. */
  gal_threads_spin_off(warp_pixelarea_onthread, wa, ntiles,
                       wa->numthreads, input->minmapsize,
                       input->quietmmap);

  /* Clean up. */
  gal_data_array_free(wa->tiles, ntiles, 0);
}
