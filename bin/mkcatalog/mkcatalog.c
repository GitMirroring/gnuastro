/*********************************************************************
MakeCatalog - Make a catalog from an input and labeled image.
MakeCatalog is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2015-2025 Free Software Foundation, Inc.

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
#include <float.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>

#include <gnuastro/git.h>
#include <gnuastro/wcs.h>
#include <gnuastro/data.h>
#include <gnuastro/fits.h>
#include <gnuastro/match.h>
#include <gnuastro/units.h>
#include <gnuastro/kdtree.h>
#include <gnuastro/threads.h>
#include <gnuastro/pointer.h>
#include <gnuastro/dimension.h>
#include <gnuastro/statistics.h>
#include <gnuastro/permutation.h>

#include <gnuastro-internal/timing.h>
#include <gnuastro-internal/checkset.h>

#include "main.h"
#include "mkcatalog.h"

#include "ui.h"
#include "parse.h"
#include "columns.h"
#include "outkeys.h"
#include "upperlimit.h"










/*********************************************************************/
/**************       Manage a single object       *******************/
/*********************************************************************/
static void
mkcatalog_clump_starting_index(struct mkcatalog_passparams *pp)
{
  struct mkcatalogparams *p=pp->p;

  /* Lock the mutex if we are working on more than one thread. NOTE: it is
     very important to keep the number of operations within the mutex to a
     minimum so other threads don't get delayed. */
  if(p->cp.numthreads>1)
    pthread_mutex_lock(&p->mutex);

  /* Put the current total number of rows filled into the output, then
     increment the total number by the number of clumps. */
  pp->clumpstartindex = p->clumprowsfilled;
  p->clumprowsfilled += pp->clumpsinobj;

  /* Unlock the mutex (if it was locked). */
  if(p->cp.numthreads>1)
    pthread_mutex_unlock(&p->mutex);
}




/* Vector allocation (short name is used since it is repeated a lot). */
static void
mkcatalog_vec_alloc(struct mkcatalog_passparams *pp, size_t onum,
                    size_t vnum, uint8_t type)
{
  if( pp->p->oiflag[ onum ] )
    gal_data_initialize(&pp->vector[vnum], 0, type, 1,
                        &(pp->p->objects->dsize[0]), NULL, 1,
                        pp->p->cp.minmapsize, pp->p->cp.quietmmap,
                        NULL, NULL, NULL);
}





/* Allocate all the necessary space. */
static void
mkcatalog_single_object_init(struct mkcatalogparams *p,
                             struct mkcatalog_passparams *pp)
{
  uint8_t *oif=p->oiflag;
  size_t ndim=p->objects->ndim;
  uint8_t i32=GAL_TYPE_INT32, f64=GAL_TYPE_FLOAT64; /* For short lines.*/

  /* Initialize the mkcatalog_passparams elements. */
  pp->p               = p;
  pp->clumpstartindex = 0;
  pp->rng             = p->rng ? gsl_rng_clone(p->rng) : NULL;
  pp->oi              = gal_pointer_allocate(GAL_TYPE_FLOAT64,
                                             OCOL_NUMCOLS, 0, __func__,
                                             "pp->oi");

  /* If we have second order measurements, allocate the array keeping the
     temporary shift values for each object of this thread. Note that the
     clumps catalog (if requested), will have the same measurements, so its
     just enough to check the objects. */
  pp->shift = ( (    oif[ OCOL_GXX ]
                  || oif[ OCOL_GYY ]
                  || oif[ OCOL_GXY ]
                  || oif[ OCOL_VXX ]
                  || oif[ OCOL_VYY ]
                  || oif[ OCOL_VXY ] )
                ? gal_pointer_allocate(GAL_TYPE_SIZE_T, ndim, 0, __func__,
                                       "pp->shift")
                : NULL );

  /* If we have upper-limit mode, then allocate the container to keep the
     values to calculate the standard deviation. */
  if(p->upperlimit)
    {
      /* Allocate the space to keep the upper-limit values. */
      pp->up_vals = gal_data_alloc(NULL, GAL_TYPE_FLOAT32, 1, &p->upnum,
                                  NULL, 0, p->cp.minmapsize,
                                  p->cp.quietmmap, NULL, NULL, NULL);

      /* Set the blank checked flag to 1. By definition, this dataset won't
         have any blank values. Also 'flag' is initialized to '0'. So we
         just have to set the checked flag ('GAL_DATA_FLAG_BLANK_CH') to
         one to inform later steps that there are no blank values. */
      pp->up_vals->flag |= GAL_DATA_FLAG_BLANK_CH;
    }
  else
    pp->up_vals=NULL;

  /* If any vector measurements are necessary, do the necessary
     allocations: first the general array, to keep all that are necessary,
     then the individual ones. */
  pp->vector = ( (    oif[ OCOL_NUMINSLICE         ]
                   || oif[ OCOL_SUMINSLICE         ]
                   || oif[ OCOL_NUMALLINSLICE      ]
                   || oif[ OCOL_SUMVARINSLICE      ]
                   || oif[ OCOL_SUMPROJINSLICE     ]
                   || oif[ OCOL_NUMPROJINSLICE     ]
                   || oif[ OCOL_SUMPROJVARINSLICE  ]
                   || oif[ OCOL_NUMOTHERINSLICE    ]
                   || oif[ OCOL_SUMOTHERINSLICE    ]
                   || oif[ OCOL_SUMOTHERVARINSLICE ]
                   || oif[ OCOL_NUMALLOTHERINSLICE ] )
                 ? gal_data_array_calloc(VEC_NUM)
                 : NULL );
  mkcatalog_vec_alloc(pp, OCOL_NUMINSLICE,      VEC_NUMINSLICE,     i32);
  mkcatalog_vec_alloc(pp, OCOL_NUMALLINSLICE,   VEC_NUMALLINSLICE,  i32);
  mkcatalog_vec_alloc(pp, OCOL_NUMPROJINSLICE,  VEC_NUMPROJINSLICE, i32);
  mkcatalog_vec_alloc(pp, OCOL_NUMOTHERINSLICE, VEC_NUMOTHERINSLICE,i32);
  mkcatalog_vec_alloc(pp, OCOL_SUMINSLICE,      VEC_SUMINSLICE,     f64);
  mkcatalog_vec_alloc(pp, OCOL_SUMVARINSLICE,   VEC_SUMVARINSLICE,  f64);
  mkcatalog_vec_alloc(pp, OCOL_SUMPROJINSLICE,  VEC_SUMPROJINSLICE, f64);
  mkcatalog_vec_alloc(pp, OCOL_SUMOTHERINSLICE, VEC_SUMOTHERINSLICE,f64);
  mkcatalog_vec_alloc(pp, OCOL_SUMOTHERVARINSLICE, VEC_SUMOTHERVARINSLICE,
                      f64);
  mkcatalog_vec_alloc(pp, OCOL_NUMALLOTHERINSLICE, VEC_NUMALLOTHERINSLICE,
                      i32);
  mkcatalog_vec_alloc(pp, OCOL_SUMPROJVARINSLICE, VEC_SUMPROJVARINSLICE,
                      f64);
}





/* Fill column values with blanks in case the label does not actually exist
   in the input table (can happen with '--inbetweenints'). This is based on
   the 'columns_fill' function. */
static void
mkcatalog_single_object_not_exist(struct mkcatalog_passparams *pp)
{
  void *vstart;
  size_t i, oind;
  gal_data_t *column;
  struct mkcatalogparams *p=pp->p;

  /* If the two metadata arrays are allocated, set them for this column. */
  if(p->numclumps_c)
    p->numclumps_c[ p->obj_to_tile
                    ? p->obj_to_tile[pp->object]
                    : pp->object-1 ] = 0;

  /* Set this object's rows in all columns to blank. */
  for(column=p->objectcols; column!=NULL; column=column->next)
    {
      /* This object's row index. */
      oind=pp->object-1;

      /* Write single-valued columns. */
      switch(column->ndim)
        {
        case 1: /* Single value columns. */
          gal_blank_write(gal_pointer_increment(column->array, oind,
                                                column->type),
                          column->type);
          break;

        case 2: /* Vector columns. */
          vstart=gal_pointer_increment(column->array,
                                      oind*column->dsize[1],
                                      column->type);
          for(i=0;i<column->dsize[1];++i)
            gal_blank_write(gal_pointer_increment(vstart, i,
                                                  column->type),
                            column->type);
          break;

        default:
          error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at "
                "'%s' to find and fix the problem. 'col->ndim' "
                "should be either 1 or 2, but has a value of %zu",
                __func__, PACKAGE_BUGREPORT, column->ndim);
        }
    }
}





/* Each thread will call this function once. It will go over all the
   objects that are assigned to it. */
static void *
mkcatalog_single_object(void *in_prm)
{
  struct gal_threads_params *tprm=(struct gal_threads_params *)in_prm;
  struct mkcatalogparams *p=(struct mkcatalogparams *)(tprm->params);

  size_t i;
  struct mkcatalog_passparams pp;

  /* Initialize and allocate all the necessary values. */
  mkcatalog_single_object_init(p, &pp);

  /* Fill the desired columns for all the objects given to this thread. */
  for(i=0; tprm->indexs[i]!=GAL_BLANK_SIZE_T; ++i)
    {
      /* For easy reading. Note that the object IDs start from one while
         the tile indexs (that we parallelize over) start from 0. */
      pp.ci = NULL;
      pp.tile = &p->tiles[ tprm->indexs[i] ];
      pp.object = ( p->obj_from_tile
                    ? p->obj_from_tile[ tprm->indexs[i] ]
                    : tprm->indexs[i] + 1 );

      /* When the 'ndim' of the tile is zero, this label did not exist in
         the image (only present because '--inbetweenints' was called). In
         this case, simply set all the column values to blank and go to the
         next object. */
      if(pp.tile->ndim==0)
        { mkcatalog_single_object_not_exist(&pp); continue; }

      /* Initialize the parameters for this object/tile. */
      parse_initialize(&pp);

      /* Get the first pass information. */
      parse_objects(&pp);

      /* Currently the second pass is only necessary when there is a clumps
         image. */
      if(p->clumps)
        {
          /* Allocate space for the properties of each clump. */
          pp.ci = gal_pointer_allocate(GAL_TYPE_FLOAT64,
                                       pp.clumpsinobj * CCOL_NUMCOLS, 1,
                                       __func__, "pp.ci");

          /* Get the starting row of this object's clumps in the final
             catalog. This index is also necessary for the unique random
             number generator seeds of each clump. */
          mkcatalog_clump_starting_index(&pp);

          /* Get the second pass information. */
          parse_clumps(&pp);
        }

      /* If an order-based calculation is requested, another pass is
         necessary. */
      if(    p->oiflag[ OCOL_MEDIAN        ]
          || p->oiflag[ OCOL_MAXIMUM       ]
          || p->oiflag[ OCOL_HALFMAXSUM    ]
          || p->oiflag[ OCOL_HALFMAXNUM    ]
          || p->oiflag[ OCOL_HALFSUMNUM    ]
          || p->oiflag[ OCOL_SIGCLIPNUM    ]
          || p->oiflag[ OCOL_SIGCLIPSTD    ]
          || p->oiflag[ OCOL_SIGCLIPMEAN   ]
          || p->oiflag[ OCOL_FRACMAX1NUM   ]
          || p->oiflag[ OCOL_FRACMAX2NUM   ]
          || p->oiflag[ OCOL_SIGCLIPMEDIAN ])
        parse_order_based(&pp);

      /* Calculate the upper limit magnitude (if necessary). */
      if(p->upperlimit) upperlimit_calculate(&pp);

      /* Write the pass information into the columns. */
      columns_fill(&pp);

      /* Clean up for this object. */
      if(pp.ci) free(pp.ci);
    }

  /* Clean up. */
  free(pp.oi);
  free(pp.shift);
  gal_data_free(pp.up_vals);
  if(pp.rng) gsl_rng_free(pp.rng);
  gal_data_array_free(pp.vector, VEC_NUM, 1);

  /* Wait until all the threads finish and return. */
  if(tprm->b) pthread_barrier_wait(tprm->b);
  return NULL;
}




















/*********************************************************************/
/********         Processing after threads finish        *************/
/*********************************************************************/
/* Convert internal image coordinates to WCS for table.

   Note that from the beginning (during the passing steps), we saved FITS
   coordinates. Also note that we are doing the conversion in place. */
static void
mkcatalog_wcs_conversion(struct mkcatalogparams *p)
{
  gal_data_t *c;
  gal_data_t *column;

  /* Flux weighted center positions for clumps and objects. */
  if(p->wcs_vo)
    {
      gal_wcs_img_to_world(p->wcs_vo, p->objects->wcs, 1);
      if(p->wcs_vc)
        gal_wcs_img_to_world(p->wcs_vc, p->objects->wcs, 1);
    }


  /* Geometric center positions for clumps and objects. */
  if(p->wcs_go)
    {
      gal_wcs_img_to_world(p->wcs_go, p->objects->wcs, 1);
      if(p->wcs_gc)
        gal_wcs_img_to_world(p->wcs_gc, p->objects->wcs, 1);
    }


  /* All clumps flux weighted center. */
  if(p->wcs_vcc)
    gal_wcs_img_to_world(p->wcs_vcc, p->objects->wcs, 1);


  /* All clumps geometric center. */
  if(p->wcs_gcc)
    gal_wcs_img_to_world(p->wcs_gcc, p->objects->wcs, 1);


  /* Go over all the object columns and fill in the values. */
  for(column=p->objectcols; column!=NULL; column=column->next)
    {
      /* Definitions */
      c=NULL;

      /* Set 'c' for the columns that must be corrected. Note that this
         'switch' statement doesn't need any 'default', because there are
         probably columns that don't need any correction. */
      switch(column->status)
        {
        case UI_KEY_W1:           c=p->wcs_vo;                break;
        case UI_KEY_W2:           c=p->wcs_vo->next;          break;
        case UI_KEY_W3:           c=p->wcs_vo->next->next;    break;
        case UI_KEY_GEOW1:        c=p->wcs_go;                break;
        case UI_KEY_GEOW2:        c=p->wcs_go->next;          break;
        case UI_KEY_GEOW3:        c=p->wcs_go->next->next;    break;
        case UI_KEY_CLUMPSW1:     c=p->wcs_vcc;               break;
        case UI_KEY_CLUMPSW2:     c=p->wcs_vcc->next;         break;
        case UI_KEY_CLUMPSW3:     c=p->wcs_vcc->next->next;   break;
        case UI_KEY_CLUMPSGEOW1:  c=p->wcs_gcc;               break;
        case UI_KEY_CLUMPSGEOW2:  c=p->wcs_gcc->next;         break;
        case UI_KEY_CLUMPSGEOW3:  c=p->wcs_gcc->next->next;   break;
        }

      /* Copy the elements into the output column. */
      if(c)
        memcpy(column->array, c->array,
               column->size*gal_type_sizeof(c->type));
    }


  /* Go over all the clump columns and fill in the values. */
  for(column=p->clumpcols; column!=NULL; column=column->next)
    {
      /* Definitions */
      c=NULL;

      /* Set 'c' for the columns that must be corrected. Note that this
         'switch' statement doesn't need any 'default', because there are
         probably columns that don't need any correction. */
      switch(column->status)
        {
        case UI_KEY_W1:           c=p->wcs_vc;                break;
        case UI_KEY_W2:           c=p->wcs_vc->next;          break;
        case UI_KEY_W3:           c=p->wcs_vc->next->next;    break;
        case UI_KEY_GEOW1:        c=p->wcs_gc;                break;
        case UI_KEY_GEOW2:        c=p->wcs_gc->next;          break;
        case UI_KEY_GEOW3:        c=p->wcs_gc->next->next;    break;
        }

      /* Copy the elements into the output column. */
      if(c)
        memcpy(column->array, c->array,
               column->size*gal_type_sizeof(c->type));
    }
}





/* Since all the measurements were done in parallel (and we didn't know the
   number of clumps per object a-priori), the clumps informtion is just
   written in as they are measured. Here, we'll sort the clump columns by
   object ID (when it is requested). */
static void
sort_clumps_by_objid(struct mkcatalogparams *p)
{
  gal_data_t *col;
  size_t i, j, tind, *permute, *rowstart;
  size_t ptwc=GAL_BLANK_SIZE_T; /* previous tile with clump. */

  /* Make sure everything is fine. */
  if(p->hosttind_c==NULL || p->numclumps_c==NULL)
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix "
          "the problem. 'p->hostobjid_c' and 'p->numclumps_c' must not "
          "be NULL.", __func__, PACKAGE_BUGREPORT);

  /* Allocate the necessary arrays. */
  rowstart=gal_pointer_allocate(GAL_TYPE_SIZE_T, p->numtiles, 0,
                                __func__, "rowstart");
  permute=gal_pointer_allocate(GAL_TYPE_SIZE_T, p->numclumps, 0,
                               __func__, "permute");

  /* The tile/object array is already sorted by tile ID. So we should just
     add up the number of clumps in each object to find the row where each
     object's clumps should start from in the final sorted clumps
     catalog. */
  for(i=0;i<p->numtiles;++i)
    {
      /* This tile has clumps, use the previous tile's total clumps to find
         the starting row of this tile's clumps. */
      if(p->numclumps_c[i])
        {
          rowstart[i] = ( ptwc==GAL_BLANK_SIZE_T
                          ? 0
                          : rowstart[ptwc]+p->numclumps_c[ptwc] );
          ptwc=i;

          /* For a check.
          printf("%s: %zu, %zu, %zu\n", __func__, i,
                 p->numclumps_c[i], rowstart[i]);
          */
        }

      /* There are no clumps on this tile, so we'll just put a blank value
         to create a segmentation fault if not used properly (better than a
         logical error). */
      else rowstart[i] = GAL_BLANK_SIZE_T;
    }

  /* Fill the permutation array to rearrange the clump rows. Note that WE
     KNOW that all the clumps within for one objects are already
     immediately after each other. */
  i=0;
  while(i<p->numclumps)
    {
      /* Find the tile-id of this clump and depending on how many clumps
         there are in this tile, increment the clump-index ('i'). */
      tind = p->hosttind_c[i];
      if(tind==GAL_BLANK_SIZE_T)
        error(EXIT_FAILURE, 0, "the number of clumps in the image does "
              "not match the 'NUMLABS' keyword of the input clumps HDU. "
              "This can happen when the clumps have been derived from "
              "Segment's output over one image while the object labels "
              "were derived from another image. As a result, there are "
              "multiple separate clumps with the same label but not with "
              "different object labels. If this the case, the safest "
              "solution is to re-run Segment with the '--noobjects' "
              "option for generating the clump labels independent of any "
              "objects (which will also be faster!). If you can no "
              "longer run Segment, you can re-label the clumps with this "
              "command: 'astarithmetic file-with-clumps.fits -hCLUMPS 0 "
              "gt 2 connected-components --output=out.fits' and replace "
              "the output with the clumps image you gave in this run");
      for(j=0; j<p->numclumps_c[tind]; ++j)
        {
          permute[i++] = rowstart[tind] + j;

          /* For a check.
             printf("%s: clump-%zu --> obj-%zu: %zu\n", __func__, i-1,
                    tind, permute[i-1]);
          */
        }
    }

  /* Permute all the clump columns. */
  for(col=p->clumpcols; col!=NULL; col=col->next)
      gal_permutation_apply_inverse(col, permute);

  /* Clean up */
  free(permute);
  free(rowstart);
}





/* Write the produced columns into the output */
static void
mkcatalog_write_outputs(struct mkcatalogparams *p)
{
  char *out=p->cp.output;
  int tf=p->cp.tableformat;
  gal_fits_list_key_t *tmp;
  gal_list_str_t *comments=NULL;
  char *ohdu="OBJECTS", *chdu="CLUMPS", *cuhdu="CHECK-UPPERLIMIT";

  /* If a catalog is to be generated. */
  if(p->objectcols)
    {
      /* Reverse the comments list (so it is printed in the same order
         here), write the objects catalog and free the comments. */
      gal_list_str_reverse(&comments);
      gal_table_write(p->objectcols, NULL, NULL, tf, out, ohdu, 0, 1);
      gal_list_str_free(comments, 1);

      /* CLUMPS catalog */
      if(p->clumps)
        {
          /* Reverse the comments list (so it is printed in the same order
             here), write the objects catalog and free the comments. */
          gal_list_str_reverse(&comments);
          gal_table_write(p->clumpcols, NULL, comments, tf, out, chdu,
                          0, 1);
          gal_list_str_free(comments, 1);
        }

      /* Upper limit check table. */
      if(p->upcheck)
        gal_table_write(p->upcheck, NULL, comments, tf, out, cuhdu, 0, 1);
    }

  /* Configuration information. */
  gal_fits_key_write_filename("input", p->objectsfile, &p->cp.ckeys,
                              1, p->cp.quiet);

  /* MakeCatalog's analysis/measurement headers after the configuration
     (input options) headers. */
  for(tmp=p->cp.ckeys; tmp!=NULL; tmp=tmp->next)
    if(tmp->next==NULL) { tmp->next=outkeys_write(p); break; }

  /* Write all the keywords in the 0-th HDU of the output. This function
     also frees the whole key list. */
  gal_fits_key_write(p->cp.ckeys, out, "0", "NONE", 1, 0);

  /* Inform the user */
  if(!p->cp.quiet)
    {
      printf("  - Output written to (HDUs listed below): %s\n"
             "    - OBJECTS: measurements on the object labels.\n",
             out);
      if(p->clumps)
        printf("    - CLUMPS: measurements on the clump labels.\n");
      if(p->upcheck)
        printf("    - CHECK-UPPERLIMIT: measurements on the "
               "clump labels.\n");
    }
}




















/*********************************************************************/
/*****************       Top-level function        *******************/
/*********************************************************************/
void
mkcatalog(struct mkcatalogparams *p)
{
  /* When more than one thread is to be used, initialize the mutex: we need
     it to assign a column to the clumps in the final catalog. */
  if( p->cp.numthreads > 1 ) pthread_mutex_init(&p->mutex, NULL);

  /* Do the processing on each thread. */
  gal_threads_spin_off(mkcatalog_single_object, p, p->numtiles,
                       p->cp.numthreads, p->cp.minmapsize,
                       p->cp.quietmmap);

  /* Post-thread processing, for example to convert image coordinates to RA
     and Dec. */
  mkcatalog_wcs_conversion(p);

  /* If the columns need to be sorted (by object ID), then some adjustments
     need to be made to the clumps catalog.. */
  if(p->hosttind_c) sort_clumps_by_objid(p);

  /* Write the filled columns into the output. */
  mkcatalog_write_outputs(p);

  /* Destroy the mutex. */
  if( p->cp.numthreads>1 ) pthread_mutex_destroy(&p->mutex);
}
