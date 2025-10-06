/*********************************************************************
MakeCatalog - Make a catalog from an input and labeled image.
MakeCatalog is part of GNU Astronomy Utilities (Gnuastro) package.

Authors:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
     Giacomo Lorenzetti <glorenzetti@cefca.es>
Copyright (C) 2025-2025 Free Software Foundation, Inc.

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
#include <string.h>

#include <gnuastro/wcs.h>
#include <gnuastro/fits.h>
#include <gnuastro/units.h>
#include <gnuastro/match.h>
#include <gnuastro/kdtree.h>
#include <gnuastro/pointer.h>
#include <gnuastro/statistics.h>

#include <gnuastro-internal/checkset.h>

#include <main.h>

#include <ui.h>
#include <mkcatalog.h>
#include <upperlimit.h>




void
outkeys_numeric(gal_fits_list_key_t **keylist, void *number,
                uint8_t type, char *namei, char *commenti, char *uniti)
{
  void *value;
  char *name, *unit, *comment;

  /* Allocate and copy the value. */
  value=gal_pointer_allocate(type, 1, 0, __func__, "value");
  memcpy(value, number, gal_type_sizeof(type));

  /* Allocate and copy the strings. */
  gal_checkset_allocate_copy(namei, &name);
  gal_checkset_allocate_copy(uniti, &unit);
  gal_checkset_allocate_copy(commenti, &comment);

  /* Add this keyword to the list. */
  gal_fits_key_list_add_end(keylist, type, name, 1, value, 1, comment, 1,
                            unit, 0);
}





static void
outkeys_sbl_nml(struct mkcatalogparams *p,
                gal_fits_list_key_t **keylist, int m0s1)
{
  char *rname = m0s1 ? "SBL" : "NML";
  float fvalue, pixarea=p->pixelarcsecsq;
  float area = m0s1 ? p->sblarea : p->nmlarea;
  float sigma = m0s1 ? p->sblsigma : p->nmlsigma;
  char *rcom = ( m0s1
                 ? "Surface brightness limit."
                 : "Noise-based magnitude limit." );

  /* The title keyword. */

  if( !isnan(p->medstd) && !isnan(sigma) )
    {
      /* Only print magnitudes if a zeropoint is given. */
      if( !isnan(p->zeropoint) )
        {
          /* Only print the SBL in fixed area if a WCS is present and a
             pixel area could be deduced. */
          if( !isnan(pixarea) )
            {
              fvalue=gal_units_counts_to_mag(
                        m0s1
                        ? ( sigma * p->medstd / sqrt( area * pixarea) )
                        : ( sigma * p->medstd * area / pixarea ),
                        p->zeropoint);
              outkeys_numeric(keylist, &fvalue, GAL_TYPE_FLOAT32, rname,
                              rcom, m0s1?"mag/arcsec^2":"mag");
            }
          else
            gal_fits_key_list_fullcomment_add_end(keylist, "Can not "
                   "report SBL or NML because the input doesn't "
                   "have a world coordinate system (WCS), or the "
                   "first two coordinates of the WCS weren't angular "
                   "positions in units of degrees.", 0);
        }
      else
        gal_fits_key_list_fullcomment_add_end(keylist, "Can not write "
               "SBL or NML because no '--zeropoint' has been given.", 0);
    }
  else
    {
      gal_fits_key_list_fullcomment_add_end(keylist, "Can not report "
             "SBL or NML because none of the requested columns needed "
             "the standard deviation (so it was not read!). Please "
             "call with '--meta-measures'.", 0);
      gal_fits_key_list_fullcomment_add_end(keylist, "", 0);
    }
}





static void
outkeys_metameasure_noisebased(struct mkcatalogparams *p,
                               gal_fits_list_key_t **keylist)
{
  gal_fits_key_list_title_add_end(keylist,
                                  "Noise-based metameasures", 0);
  outkeys_sbl_nml(p, keylist, 1);
  outkeys_sbl_nml(p, keylist, 0);
}





/* Write the check table for the confusion limit. */
static void
outkeys_confusion_limit_dist(struct mkcatalogparams *p, gal_data_t *wcscrd,
                             gal_data_t *matched)
{
  float *dp, *dpf;
  size_t i, *ind=matched->array;
  double ps, *mr=NULL, *md=NULL, *dd=NULL;
  double *r=wcscrd->array, *d=wcscrd->next->array;
  gal_data_t *dista=NULL, *distp=NULL, *chtab=NULL, *dists=matched->next;

  /* Allocate the arrays and pointers in the matched-ra and matched-dec
     columns. */
  if(p->cnlcheck)
    {
      /* Prepare the list. */
      gal_list_data_add_alloc(&chtab, NULL, GAL_TYPE_FLOAT64, 1,
                              wcscrd->dsize, NULL, 0, p->cp.minmapsize,
                              p->cp.quietmmap, "DEC-NEAREST", "deg",
                              "Declination of nearest label.");
      gal_list_data_add_alloc(&chtab, NULL, GAL_TYPE_FLOAT64, 1,
                              wcscrd->dsize, NULL, 0, p->cp.minmapsize,
                              p->cp.quietmmap, "RA-NEAREST", "deg",
                              "Right ascension of nearest label.");

      /* To keep distances in pixels. */
      distp=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, 1, dists->dsize, NULL,
                           0, dists->minmapsize, dists->quietmmap,
                           "DISTANCE_IMG", "pixels", "Same as DISTANCE_WCS "
                           "but in units of pixels.");

      /* Set the pointers for easy reading. */
      mr=chtab->array;
      md=chtab->next->array;
    }

  /* Fill the matched columns with distance on sphere. */
  dd=dists->array;
  for(i=0;i<wcscrd->size;++i)
    {
      if(p->cnlcheck) { mr[i]=r[ind[i]]; md[i]=d[ind[i]]; }
      dd[i]=gal_wcs_angular_distance_deg(r[i],      d[i],
                                         r[ind[i]], d[ind[i]])*3600;

      /* For a check:
      if(i==0)
        printf("%s: (%f,%f) --> (%f,%f) with %f\n", __func__,
               r[i], d[i], r[ind[i]], d[ind[i]], dd[i]);
      */
    }

  /* Put the new columns at the end of the existing coordinates and write
     it out. */
  if(p->cnlcheck)
    {
      /* Copy the distances array into float32 (the extra precision is just
         floating-point errors on the arcsec/pixel level): on the level of
         1-millionth of a pixel/arcsec! */
      dista=gal_data_copy_to_new_type(dists, GAL_TYPE_FLOAT32);

      /* Correct the distf metadata. */
      free(dista->name);
      free(dista->unit);
      free(dista->comment);
      gal_checkset_allocate_copy("arcsec", &dista->unit);
      gal_checkset_allocate_copy("DISTANCE_WCS", &dista->name);
      gal_checkset_allocate_copy("Distance on a sphere to nearest point.",
                                 &dista->comment);

      /* Fill up the column with distances in pixels. */
      ps=sqrt(p->pixelarcsecsq);
      dpf=(dp=distp->array)+distp->size;
      do { *dp=*dd/ps; ++dd; } while (++dp<dpf);

      /* Put the floating-point and pixel distances at the end of the check
         table. Then put the whole check table after the input coordinates
         for writing. */
      dista->next=distp;
      chtab->next->next=dista;
      wcscrd->next->next=chtab;

      /* Write the table. */
      gal_table_write(wcscrd, NULL, NULL, p->cp.tableformat,
                      p->cp.output, "CONFUSION-LIMIT", 0, 0);

      /* Remove the extra pointers for the test table. */
      wcscrd->next->next=NULL;

      /* Clean up the allocated tables for this check. */
      gal_list_data_free(chtab);
    }
}





#define OUTKEYS_CNL_NUMBER 5
static void
outkeys_confusion_limit_write(struct mkcatalogparams *p,
                              gal_fits_list_key_t **keylist,
                              gal_data_t *dists)
{
  size_t i;
  float tmp, *fval;
  gal_data_t *stat;
  char *name, *comment;
  double ps=sqrt(p->pixelarcsecsq);
  uint8_t percentiles[OUTKEYS_CNL_NUMBER]={5, 25, 50, 75, 95};

  /* Necessary allocations. */
  fval=gal_pointer_allocate(GAL_TYPE_FLOAT32, OUTKEYS_CNL_NUMBER, 0,
                            __func__, "fval");
  name=gal_pointer_allocate(GAL_TYPE_UINT8, 9, 0, __func__, "fval");
  comment=gal_pointer_allocate(GAL_TYPE_UINT8, FLEN_VALUE+1, 0,
                               __func__, "fval");

  /* Fill the values. We are claculating the values here, but not writing
     them, so the 'CNL' can be the first keyword. */
  for(i=0; i<OUTKEYS_CNL_NUMBER; ++i)
    {
      /* Calculate the statistic (in pixels). */
      stat=gal_statistics_quantile(dists, percentiles[i]/100.0f, 1);
      fval[i]=((double *)(stat->array))[0]/ps;
    }

  /* Write the final CNL in units of arcsec. */
  tmp = fval[3] - fval[1];
  outkeys_numeric(keylist, &tmp, GAL_TYPE_FLOAT32, "CNL",
                  "Confusion limit: CNLP75-CNLP25.", "pixels");
  printf("%s: %f, %f\n", __func__, fval[1], fval[3]);

  /* Write the values. */
  for(i=0; i<OUTKEYS_CNL_NUMBER; ++i)
    {
      /* Set the strings to write. */
      sprintf(name, "CNLP%02u", percentiles[i]);
      sprintf(comment, "Quant. %.2f of. dist. to nearest lab.",
              percentiles[i]/100.0f);

      /* Write the keyword. */
      outkeys_numeric(keylist, &fval[i], GAL_TYPE_FLOAT32, name,
                      comment, "pixels");
    }

  /* Clean up. */
  free(fval);
  free(name);
  free(comment);
}




/* Confusion limit: create a k-d tree using the image coordinates of the
   clumps, then find the nearest neighbour of each clump and return a struct
   with the mad-clipped stats. */
static void
outkeys_confusion_limit(struct mkcatalogparams *p, gal_data_t *wcscrd,
                        gal_fits_list_key_t **keylist)
{
  size_t root;
  size_t nummatched;
  uint8_t *flag=NULL;
  gal_data_t *kdtree, *matched;

  /* Make sure the inputs are good (we plan to later take this function
     outside of MakeCatalog). */
  if(wcscrd->next==NULL || wcscrd->next->next!=NULL)
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at '%s' "
          "to fix the problem. The 'coords' argument should only "
          "have two columns, but it has %zu columns", __func__,
          PACKAGE_BUGREPORT, gal_list_data_number(wcscrd));
  if( wcscrd->type!=GAL_TYPE_FLOAT64
      && wcscrd->next->type!=GAL_TYPE_FLOAT64)
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at '%s' "
          "to fix the problem. The 'coords' argument should all be "
          "float64, but they are '%s' and '%s' respectively", __func__,
          PACKAGE_BUGREPORT, gal_type_name(wcscrd->type, 1),
          gal_type_name(wcscrd->next->type, 1));

  /* Initialize the k-d tree. */
  kdtree=gal_kdtree_create(wcscrd, &root);

  /* Find the nearest neighbour of each clump */
  matched=gal_match_kdtree(wcscrd, wcscrd, kdtree, root,
                           GAL_MATCH_ARRANGE_OUTER, NULL,
                           p->cp.numthreads, p->cp.minmapsize,
                           p->cp.quietmmap, &nummatched, &flag, 1);

  /* This should be tested later. */
  if(flag)
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at '%s' "
          "to fix the problem. The match for finding the confusion "
          "limit should not have produced any flags", __func__,
          PACKAGE_BUGREPORT);

  /* Calculate spherical distances and write the check table in the output
     if requested */
  outkeys_confusion_limit_dist(p, wcscrd, matched);

  /* Calculate the various quantiles. 'matched->next' is the array to use
     (was corrected in 'outkeys_confusion_limit_dist'). */
  outkeys_confusion_limit_write(p, keylist, matched->next);

  /* Clean up and return. */
  gal_list_data_free(matched);
  gal_list_data_free(kdtree);
}





/* Confusion limit highlevel function. It extract the info from the clumps
   catalog and, if possible, writes the headers related to confusion limit */
static void
outkeys_confusion(struct mkcatalogparams *p, gal_fits_list_key_t **keylist)
{
  int cleancrds=0;
  gal_data_t *rn=NULL, *dn=NULL, *clim=NULL, *cat=NULL;
  gal_data_t *c1, *c2, *tmp, *x=NULL, *y=NULL, *r=NULL, *d=NULL;

  /* Title keyword (always gets printed). */
  gal_fits_key_list_title_add_end(keylist, "Confusion limit "
                                  "(CNL)", 0);

  /* In case the input is not two dimensional, print a comment and return
     (other dimensions need to be added). */
  if(p->objects->ndim!=2)
    {
      gal_fits_key_list_fullcomment_add_end(keylist, "Can't write "
               "confusion limit values because it is currently only "
               "implemented for 2D datasets. Please inform us at "
               PACKAGE_BUGREPORT" to add this feature", 0);
      return;
    }

  /* Set the main catalog for measuring the confusion limit and return if
     we need clumps and no clumps data set was present. */
  cat = p->cnlwithobjects ? p->objectcols : p->clumpcols;
  if(cat==NULL)
    {
      gal_fits_key_list_fullcomment_add_end(keylist, "No confusion "
               "cimit calculations because a clumps catalog was not "
               "requested. If you used Gnuastro's Segment to generate "
               "your input labels, re-run mkcatalog with --clumpscat "
               "to solve the issue. Otherwise, use '--cnl-with-objects'",
               0);
      return;
    }

  /* Check available coordinates. */
  for(tmp=cat; tmp!=NULL; tmp=tmp->next)
    {
      if(tmp->status==UI_KEY_X)  x=tmp;
      if(tmp->status==UI_KEY_Y)  y=tmp;
      if(tmp->status==UI_KEY_W1) r=tmp;
      if(tmp->status==UI_KEY_W2) d=tmp;
    }

  /* Set the columns to use for the distance measure in an independent
     table. RA,Dec will be used for the confusion limit. But the user may
     have only asked for X,Y! In that case, we need to create the
     independent table by converting the X,Y to RA,Dec.*/
  if ( r && d )
    { c1=r; c2=d; rn=r->next; dn=d->next; c1->next=c2; c2->next=NULL; }
  else if ( x && y )
    {
      /* Only continue if the input had WCS. */
      if(!p->objects->wcs)
        {
          gal_fits_key_list_fullcomment_add_end(keylist, "Can't write "
                   "the confusion limit (for example CNL) because the "
                   "input doesn't have a world coordinate system (WCS)",
                   0);
          return;
        }

      /* Convert the X,Y to float64 (by default, they are float32) and into
         an independent table. */
      c1=gal_data_copy_to_new_type(x, GAL_TYPE_FLOAT64);
      c2=gal_data_copy_to_new_type(y, GAL_TYPE_FLOAT64);
      c1->next=c2; c2->next=NULL;
      cleancrds=1;

      /* Convert the coordinates. It is done in-place, so we are just
         writing the returned array back into the 'c1' to avoid possible
         compiler warning ('c1' will not change!). */
      c1=gal_wcs_img_to_world(c1, p->objects->wcs, 1);
    }
  else
    {
      gal_fits_key_list_fullcomment_add_end(keylist, "Can't write "
               "Confusion Limit values because no coordinates were "
               "requested in the output catalogs (in image or WCS). "
               "Rerun MakeCatalog with '--ra --dec' or '--x --y'.", 0);
      return;
    }

  /* Compute the confusion limit and extract the necessary components. Note
     that it can crash and in that case, we will just return from this
     function (since all next steps depend on it). */
  outkeys_confusion_limit(p, c1, keylist);

  /* If the 'next' pointers of the input table were modified, set them to
     their original values. */
  if(rn) { r->next=rn; d->next=dn; }

  /* Clean up. */
  if(cleancrds) gal_list_data_free(c1);
  if(clim) gal_data_free(clim);
}





void
outkeys_infiles(struct mkcatalogparams *p, gal_fits_list_key_t **keylist)
{
  int quiet=p->cp.quiet;
  char *stdname, *stdhdu, *stdvalcom;

  /* Title for this group of keywords. */
  gal_fits_key_list_title_add_end(keylist,
                                  "Input file(s) and HDUs", 0);

  /* Object labels. */
  gal_fits_key_write_filename("INLAB", p->objectsfile, keylist, 0,
                              quiet);
  gal_fits_key_write_filename("INLABHDU", p->cp.hdu, keylist, 0,
                              quiet);

  /* Clump labels. */
  if(p->clumps)
    {
      gal_fits_key_write_filename("INCLU", p->usedclumpsfile, keylist, 0,
                                  quiet);
      gal_fits_key_write_filename("INCLUHDU", p->clumpshdu, keylist, 0,
                                  quiet);
    }

  /* Values image. */
  if(p->values)
    {
      gal_fits_key_write_filename("INVAL", p->usedvaluesfile, keylist, 0,
                                  quiet);
      gal_fits_key_write_filename("INVALHDU", p->valueshdu, keylist, 0,
                                  quiet);
    }

  /* Sky image/value. */
  if(p->sky)
    {
      if(p->sky->size==1)
        outkeys_numeric(keylist, p->sky->array, p->sky->type, "INSKYVAL",
                        "Value of Sky used (a single number).", NULL);
      else
        {
          gal_fits_key_write_filename("INSKY", p->usedskyfile, keylist, 0,
                                      quiet);
          gal_fits_key_write_filename("INSKYHDU", p->skyhdu, keylist, 0,
                                      quiet);
        }
    }

  /* Standard deviation (or variance) image. */
  if(p->variance)
    {
      stdname="INVAR"; stdhdu="INVARHDU";
      stdvalcom="Value of Sky variance (a single number).";
    }
  else
    {
      stdname="INSTD"; stdhdu="INSTDHDU";
      stdvalcom="Value of Sky STD (a single number).";
    }
  if(p->std)
    {
      if(p->std->size==1)
        outkeys_numeric(keylist, p->std->array, p->std->type, stdname,
                        stdvalcom, NULL);
      else
        {
          gal_fits_key_write_filename(stdname, p->usedstdfile, keylist, 0,
                                      quiet);
          gal_fits_key_write_filename(stdhdu, p->stdhdu, keylist, 0,
                                      quiet);
        }
    }

  /* Upper limit mask. */
  if(p->upmaskfile)
    {
      gal_fits_key_write_filename("INUPM", p->upmaskfile, keylist, 0,
                                  quiet);
      gal_fits_key_write_filename("INUPMHDU", p->upmaskhdu, keylist, 0,
                                  quiet);
    }
}





void
outkeys_pixinfo(struct mkcatalogparams *p, gal_fits_list_key_t **keylist)
{
  float fval, pixarea=p->pixelarcsecsq;

  /* Title for this group of keywords. */
  gal_fits_key_list_title_add_end(keylist,
           "Input pixel grid and value properties", 0);

  /* Pixel area. */
  if(p->objects->wcs)
    {
      if( isnan(p->pixelarcsecsq)==0 )
        {
          fval=sqrt(p->pixelarcsecsq);
          outkeys_numeric(keylist, &fval, GAL_TYPE_FLOAT32, "PIXWIDTH",
                          "Pixel width of input image.", "arcsec");
          outkeys_numeric(keylist, &pixarea, GAL_TYPE_FLOAT32, "PIXAREA",
                          "Pixel area of input image.", "arcsec^2");
        }
    }

  /* Zeropoint magnitude. */
  if( !isnan(p->zeropoint) )
    outkeys_numeric(keylist, &p->zeropoint, GAL_TYPE_FLOAT32, "ZEROPNT",
                    "Zero point (photometric calibration).", "mag");

  /* Standard devaition for noise-based measurements. */
  if( !isnan(p->medstd) )
    outkeys_numeric(keylist, &p->medstd, GAL_TYPE_FLOAT32, "MMSTD",
                    "Pixel STD for noise-based meta-measures.", NULL);

  /* The count-per-second correction. */
  if(p->cpscorr>1.0f)
    outkeys_numeric(keylist, &p->cpscorr, GAL_TYPE_FLOAT32, "CPSCORR",
                    "Counts-per-second correction.", NULL);
}





/* It is necessary to write the upperlimit parameters into the output
   tables. The same set of information will thus be necessary both in the
   upperlimit check table and also the final output. This function will do
   the job in both cases.

   Note that in the check output, the sigma-clipping information is not
   used/necessary, so to avoid confusion, we won't write it.
*/
void
outkeys_upperlimit(struct mkcatalogparams *p,
                   gal_fits_list_key_t **keylist, int withsigclip)
{
  /* Write a title for  */
  gal_fits_key_list_title_add_end(keylist, "Upper-limit (UP) parameters",
                                  0);

  /* Basic settings. */
  gal_fits_key_list_add_end(keylist, GAL_TYPE_FLOAT32, "UPSIGMA", 0,
                            &p->upnsigma, 0,
                            "Multiple of sigma to measure upper-limit.",
                            0, NULL, 0);
  gal_fits_key_list_add_end(keylist, GAL_TYPE_SIZE_T, "UPNUMBER", 0,
                            &p->upnum, 0,
                            "Number of usable random samples.", 0,
                            "counter", 0);
  gal_fits_key_list_add_end(keylist, GAL_TYPE_STRING, "UPRNGNAM", 0,
                            (void *)(p->rng_name), 0,
                            "Random number generator name.",
                            0, NULL, 0);
  outkeys_numeric(keylist, &p->rng_seed, GAL_TYPE_ULONG, "UPRNGSEE",
                  "Random number generator seed.", NULL);

  /* Range of upper-limit values. */
  if(p->uprange)
    {
      gal_fits_key_list_add_end(keylist, GAL_TYPE_SIZE_T, "UPRANGE1", 0,
                                &p->uprange[p->objects->ndim-1], 0,
                                "Range about target in axis 1.", 0,
                                "pixels", 0);
      gal_fits_key_list_add_end(keylist, GAL_TYPE_STRING, "UPRANGE2", 0,
                                &p->uprange[p->objects->ndim==2 ? 0 : 1], 0,
                                "Range about target in axis 2.", 0,
                                "pixels", 0);
      if(p->objects->ndim==3)
        gal_fits_key_list_add_end(keylist, GAL_TYPE_STRING, "UPRANGE3", 0,
                                  &p->uprange[0], 0,
                                  "Range about target in axis 3.", 0,
                                  "pixels", 0);
    }

  /* If the upper-limit measurement included sigma-clipping. */
  if(withsigclip)
    {
      gal_fits_key_list_add_end(keylist, GAL_TYPE_FLOAT64, "UPSCMLTP", 0,
                                &p->upsigmaclip[0], 0,
                                "Multiple of STD used for sigma-clipping.",
                                0, NULL, 0);
      if(p->upsigmaclip[1]>=1.0f)
        gal_fits_key_list_add_end(keylist, GAL_TYPE_FLOAT64, "UPSCNUM", 0,
                                  &p->upsigmaclip[1], 0,
                                  "Number of clips for sigma-clipping.", 0,
                                  NULL, 0);
      else
        gal_fits_key_list_add_end(keylist, GAL_TYPE_FLOAT64, "UPSCTOL", 0,
                                  &p->upsigmaclip[1], 0,
                                  "Tolerance level to sigma-clipping.", 0,
                                  NULL, 0);

    }

  /* If an upper limit check table was requested, write the object/clump
     ID(s) that the check was for.. */
  if(p->upcheck)
    {
      outkeys_numeric(keylist, &p->checkuplim[0], GAL_TYPE_INT32,
                      "UPCHKOBJ", "Object label for upper-limit "
                      "check target.", NULL);
      if( p->checkuplim[1]!=GAL_BLANK_INT32 )
        outkeys_numeric(keylist, &p->checkuplim[1], GAL_TYPE_INT32,
                        "UPCHKCLU", "Clump label for upper-limit check "
                        "target.", NULL);
    }
}





/* Write the output keywords. */
gal_fits_list_key_t *
outkeys_write(struct mkcatalogparams *p)
{
  gal_fits_list_key_t *keylist=NULL;

  /* First, add the file names and pixel grid information.*/
  outkeys_infiles(p, &keylist);
  outkeys_pixinfo(p, &keylist);

  /* Print upper-limit parameters (if there actually was any upper-limit
     request). */
  if(p->upperlimit) outkeys_upperlimit(p, &keylist, 1);

  /* Meta measurements. */
  if(p->metameasures)
    {
      /* Noise-based metameasures. */
      outkeys_metameasure_noisebased(p, &keylist);

      /* Catalog-based metameasures. */
      outkeys_confusion(p, &keylist);
    }

  /* Return the list of keywords. */
  return keylist;
}
