/*********************************************************************
Units -- Convert data from one unit to other.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Kartik Ohri <kartikohri13@gmail.com>
Contributing author(s):
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
     Pedram Ashofteh-Ardakani <pedramardakani@pm.me>
Copyright (C) 2020-2026 Free Software Foundation, Inc.

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

#include <time.h>
#include <math.h>
#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gnuastro/type.h>
#include <gnuastro/pointer.h>

#include <gnuastro-internal/checkset.h>





/**********************************************************************/
/****************      Functions to parse strings     *****************/
/**********************************************************************/
/* Parse the input string consisting of numbers separated by given
   delimiter into an array. */
int
gal_units_extract_decimal(char *convert, const char *delimiter,
                          double *args, size_t n)
{
  size_t i = 0;
  char *copy, *token, *end;

  /* Create a copy of the string to be parsed and parse it. This is because
     it will be modified during the parsing. */
  copy=strdup(convert);
  do
    {
      /* Check if the required number of arguments are passed. */
      if(i==n+1)
        {
          free(copy);
          error(0, 0, "%s: input '%s' exceeds maximum number of arguments "
                "(%zu)", __func__, convert, n);
          return 0;
        }

      /* Extract the substring till the next delimiter. */
      token=strtok(i==0?copy:NULL, delimiter);
      if(token)
        {
          /* Parse extracted string as a number, and check if it worked. */
          args[i++] = strtod (token, &end);
          if (*end && *end != *delimiter)
            {
              /* In case a warning is necessary
              error(0, 0, "%s: unable to parse element %zu in '%s'\n",
                    __func__, i, convert);
              */
              free(copy);
              return 0;
            }
        }
    }
  while(token && *token);
  free (copy);

  /* Check if the number of elements parsed. */
  if (i != n)
    {
      /* In case a warning is necessary
      error(0, 0, "%s: input '%s' must contain %lu numbers, but has "
            "%lu numbers\n", __func__, convert, n, i);
      */
      return 0;
    }

  /* Numbers are written, return successfully. */
  return 1;
}


















/**********************************************************************/
/****************      Convert string to decimal      *****************/
/**********************************************************************/

/* Parse the right ascension input as a string in form of hh:mm:ss to a
 * single decimal value calculated by (hh + mm / 60 + ss / 3600 ) * 15. */
double
gal_units_ra_to_degree(char *convert)
{
  double val[3];
  double decimal=0.0;

  /* Check whether the string is successfully parsed. */
  if(gal_units_extract_decimal(convert, ":hms", val, 3))
    {
      /* Check whether the first value is in within limits, and add it. We
         are using 'signbit(val[0])' instead of 'val[0]<0.0f' because
         'val[0]<0.0f' can't distinguish negative zero (-0.0) from an
         unsigned zero (in other words, '-0.0' will be interpretted to be
         positive). For the declinations it is possible (see the comments
         in 'gal_units_dec_to_degree'), so a user may mistakenly give that
         format in Right Ascension. */
      if(signbit(val[0]) || val[0]>24.0) return NAN;
      decimal += val[0];

      /* Check whether value of minutes is within limits, and add it. */
      if(signbit(val[1]) || val[1]>60.0) return NAN;
      decimal += val[1] / 60;

      /* Check whether value of seconds is in within limits, and add it. */
      if(signbit(val[2]) || val[2]>60.0) return NAN;
      decimal += val[2] / 3600;

      /* Convert value to degrees and return. */
      decimal *= 15.0;
      return decimal;
    }
  else return NAN;

  /* Control shouldn't reach this point. If it does, its a bug! */
  error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix the "
        "problem. Control should not reach the end of this function",
        __func__, PACKAGE_BUGREPORT);
  return NAN;
}





/* Parse the declination input as a string in form of dd:mm:ss to a decimal
 * calculated by (dd + mm / 60 + ss / 3600 ). */
double
gal_units_dec_to_degree(char *convert)
{
  int sign;
  double val[3], decimal=0.0;

  /* Parse the values in the input string. */
  if(gal_units_extract_decimal(convert, ":dms", val, 3))
    {
      /* Check whether the first value is in within limits. */
      if(val[0]<-90.0 || val[0]>90.0) return NAN;

      /* If declination is negative, the first value in the array will be
         negative and all other values will be positive. In that case, we
         set sign equal to -1. Therefore, we multiply the first value by
         sign to make it positive. The final answer is again multiplied by
         sign to make its sign same as original.

        We are using 'signbit(val[0])' instead of 'val[0]<0.0f' because
        'val[0]<0.0f' can't distinguish negative zero (-0.0) from an
        unsigned zero (in other words, '-0.0' will be interpretted to be
        positive). In the case of declination, this can happen just below
        the equator (where the declination is less than one degree), for
        example '-00d:12:34'. */
      sign = signbit(val[0]) ? -1 : 1;
      decimal += val[0] * sign;

      /* Check whether value of arc-minutes is in within limits. */
      if(signbit(val[1]) || val[1]>60.0) return NAN;
      decimal += val[1] / 60;

      /* Check whether value of arc-seconds is in within limits. */
      if (signbit(val[2]) || val[2] > 60.0) return NAN;
      decimal += val[2] / 3600;

      /* Make the sign of the decimal value same as input and return. */
      decimal *= sign;
      return decimal;
    }
  else return NAN;

  /* Control shouldn't reach this point. If it does, its a bug! */
  error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix the "
        "problem. Control should not reach the end of this function",
        __func__, PACKAGE_BUGREPORT);
  return NAN;
}




















/**********************************************************************/
/****************      Convert decimal to string      *****************/
/**********************************************************************/

/* Max-length of output string. */
#define UNITS_RADECSTR_MAXLENGTH 50

/* Parse the right ascension input as a decimal to a string in form of
   hh:mm:ss.ss . */
char *
gal_units_degree_to_ra(double decimal, int usecolon)
{
  size_t nchars;
  int hours=0, minutes=0;
  float seconds=0.0; /* For sub-second accuracy. */

  /* Allocate a long string which is large enough for string of format
     hh:mm:ss.ss and sign. */
  char *ra=gal_pointer_allocate(GAL_TYPE_UINT8, UNITS_RADECSTR_MAXLENGTH,
                                0, __func__, "ra");

  /* Check if decimal value is within bounds otherwise return error. */
  if (decimal<0 || decimal>360)
    {
      error (0, 0, "%s: value of decimal should be between be 0 and 360, "
             "but is %g\n", __func__, decimal);
      return NULL;
    }

  /* Divide decimal value by 15 and extract integer part of decimal value
     to obtain hours. */
  decimal /= 15.0;
  hours = (int)decimal;

  /* Subtract hours from decimal and multiply remaining value by 60 to
     obtain minutes. */
  minutes = (int)((decimal - hours) * 60);

  /* Subtract hours and minutes from decimal and multiply remaining value
     by 3600 to obtain seconds. */
  seconds = (decimal - hours - minutes / 60.0) * 3600;

  /* Format the extracted hours, minutes and seconds as a string with
     leading zeros if required, in hh:mm:ss format. */
  nchars = snprintf(ra, UNITS_RADECSTR_MAXLENGTH-1,
                    usecolon ? "%02d:%02d:%g" : "%02dh%02dm%g",
                    hours, minutes, seconds);
  if(nchars>UNITS_RADECSTR_MAXLENGTH)
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to address "
          "the problem. The output string has an unreasonable length of "
          "%zu characters", __func__, PACKAGE_BUGREPORT, nchars);

  /* Return the final string. */
  return ra;
}





/* Parse the declination input as a decimal to a string in form of dd:mm:ss . */
char *
gal_units_degree_to_dec(double decimal, int usecolon)
{
  size_t nchars;
  float arc_seconds=0.0;
  int sign, degrees=0, arc_minutes=0;

  /* Allocate string of fixed length which is large enough for string of
   * format hh:mm:ss.ss and sign. */
  char *dec=gal_pointer_allocate(GAL_TYPE_UINT8, UNITS_RADECSTR_MAXLENGTH,
                                 0, __func__, "ra");

  /* Check if decimal value is within bounds otherwise return error. */
  if(decimal<-90 || decimal>90)
    {
      error (0, 0, "%s: value of decimal should be between be -90 and 90, "
             "but is %g\n", __func__, decimal);
      return NULL;
    }

  /* If declination is negative, we set 'sign' equal to -1. We multiply the
     decimal by to make sure it is positive. We then extract degrees,
     arc-minutes and arc-seconds from the decimal. Finally, we add a minus
     sign in beginning of string if input was negative. */
  sign = decimal<0.0 ? -1 : 1;
  decimal *= sign;

  /* Extract integer part of decimal value to obtain degrees. */
  degrees=(int)decimal;

  /* Subtract degrees from decimal and multiply remaining value by 60 to
     obtain arc-minutes. */
  arc_minutes=(int)( (decimal - degrees) * 60 );

  /* Subtract degrees and arc-minutes from decimal and multiply remaining
     value by 3600 to obtain arc-seconds. */
  arc_seconds = (decimal - degrees - arc_minutes / 60.0) * 3600;

  /* Format the extracted degrees, arc-minutes and arc-seconds as a string
     with leading zeros if required, in hh:mm:ss format with correct
     sign. */
  nchars = snprintf(dec, UNITS_RADECSTR_MAXLENGTH-1,
                    usecolon ? "%s%02d:%02d:%g" : "%s%02dd%02dm%g",
                    sign<0?"-":"+", degrees, arc_minutes, arc_seconds);
  if(nchars>UNITS_RADECSTR_MAXLENGTH)
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to address "
          "the problem. The output string has an unreasonable length of "
          "%zu characters", __func__, PACKAGE_BUGREPORT, nchars);

  /* Return the final string. */
  return dec;
}




















/**********************************************************************/
/****************          Flux conversions           *****************/
/**********************************************************************/

/* Convert counts to magnitude using the given zeropoint. */
double
gal_units_counts_to_mag(double counts, double zeropoint)
{
  return ( counts > 0.0f
           ? ( -2.5f * log10(counts) + zeropoint )
           : NAN );
}





/* Convert magnitude to counts using the given zeropoint. */
double
gal_units_mag_to_counts(double mag, double zeropoint)
{
  return pow(10, (mag - zeropoint)/(-2.5f));
}



/* Convert apparent magnitude to luminosity. The absolute magnitude of the
   sun for different filters can be taken from Table 3 of
   https://arxiv.org/abs/1804.07788
*/
double
gal_units_mag_to_luminosity(double mag, double mag_absolute_sun,
                            double distance_modulus)
{
  return pow(10.0, (mag_absolute_sun - (mag - distance_modulus)) / 2.5);
}




double
gal_units_luminosity_to_mag(double luminosity, double mag_absolute_sun,
                            double distance_modulus)
{
  return mag_absolute_sun - 2.5*log10(luminosity) + distance_modulus;
}



double
gal_units_mag_to_sb(double mag, double area_arcsec2)
{
  return mag+2.5*log10(area_arcsec2);
}





double
gal_units_sb_to_mag(double sb, double area_arcsec2)
{
  return sb-2.5*log10(area_arcsec2);
}





double
gal_units_counts_to_sb(double counts, double zeropoint,
                       double area_arcsec2)
{
  return gal_units_mag_to_sb(
                   gal_units_counts_to_mag(counts, zeropoint),
                   area_arcsec2);
}





double
gal_units_sb_to_counts(double sb, double zeropoint,
                       double area_arcsec2)
{
  return gal_units_mag_to_counts(
                   gal_units_sb_to_mag(sb, area_arcsec2),
                   zeropoint);
}




/* Convert Pixel values to Janskys with an AB-magnitude based
   zero-point. See the "Brightness, Flux, Magnitude and Surface
   brightness". */
double
gal_units_counts_to_jy(double counts, double zeropoint_ab)
{
  return counts * 3631 * pow(10, -1 * zeropoint_ab / 2.5);
}





/* Convert Janskys to pixel values with an AB-magnitude based
   zero-point. See the "Brightness, Flux, Magnitude and Surface
   brightness". */
double
gal_units_jy_to_counts(double jy, double zeropoint_ab)
{
  return jy / 3631 / pow(10, -1 * zeropoint_ab / 2.5);
}





/* Convert counts to a custom zero point. The job of this function is
   equivalent to the double-call bellow. We just don't want to repeat some
   extra multiplication/divisions.

     gal_units_jy_to_counts(gal_units_counts_to_jy(counts, zeropoint_in),
                            custom_out)
*/
double
gal_units_zeropoint_change(double counts, double zeropoint_in,
                           double zeropoint_out)
{
  return ( counts
           * pow(10, -1 * zeropoint_in  / 2.5)
           / pow(10, -1 * zeropoint_out / 2.5) );
}





double
gal_units_counts_to_nanomaggy(double counts, double zeropoint_ab)
{
  return gal_units_zeropoint_change(counts, zeropoint_ab, 22.5);
}





double
gal_units_nanomaggy_to_counts(double counts, double zeropoint_ab)
{
  return gal_units_zeropoint_change(counts, 22.5, zeropoint_ab);
}





double
gal_units_jy_to_mag(double jy)
{
  double zp=0;
  return gal_units_counts_to_mag(gal_units_jy_to_counts(jy, zp),zp);
}





/* Converting Janskys ($f(\nu)$ or flux in units of frequency) to
   $f(\lambda)$ (or wavelength flux density). See the description of this
   operator in the book for its derivation.*/
double
gal_units_jy_to_wavelength_flux_density(double jy, double angstrom)
{
  return jy * 2.99792458e-05 / (angstrom*angstrom);
}





double
gal_units_wavelength_flux_density_to_jy(double wfd, double angstrom)
{
  return wfd * (angstrom*angstrom) / 2.99792458e-05;
}




double
gal_units_mag_to_jy(double mag)
{
  double zp=0;
  return gal_units_counts_to_jy(gal_units_mag_to_counts(mag, zp),zp);
}





double
gal_units_sblim_diff(double r_frac, double t_frac)
{
  return 2.5*log10(r_frac*sqrt(t_frac));
}




















/**********************************************************************/
/****************         Distance conversions        *****************/
/**********************************************************************/
/* Convert Astronomical Units (AU) to Parsecs (PC). From the definition of
   Parsecs, 648000/pi AU = 1 PC. The mathematical constant 'PI' is imported
   from the GSL as M_PI. So: */
double
gal_units_au_to_pc(double au)
{
  return au / 648000.0f * M_PI;
}





/* Convert Parsecs (PC) to Astronomical units (AU), see comment of
   'gal_units_au_to_pc'. */
double
gal_units_pc_to_au(double au)
{
  return au * 648000.0f / M_PI;
}





/* Convert Light-years to Parsecs, according to
   https://en.wikipedia.org/wiki/Light-year#Definitions:

   1 light-year = 9460730472580800 metres (exactly)
                ~ 9.461 petametres
                ~ 9.461 trillion kilometres (5.879 trillion miles)
                ~ 63241.077 astronomical units
                ~ 0.306601 parsecs  */
double
gal_units_ly_to_pc(double ly)
{
  return ly * 0.306601f;
}





/* Convert Parsecs to Light-years (see comment of gal_units_ly_to_pc). */
double
gal_units_pc_to_ly(double pc)
{
  return pc / 0.306601f;
}





/* Convert Astronomical Units to Light-years (see comment of
   gal_units_ly_to_pc). */
double
gal_units_au_to_ly(double au)
{
  return au / 63241.077f;
}





/* Convert Light-years to Astronomical Units (see comment of
   gal_units_ly_to_pc). */
double
gal_units_ly_to_au(double ly)
{
  return ly * 63241.077f;
}






/**********************************************************************/
/****************           Time conversions          *****************/
/**********************************************************************/

/* Fill the 'tm' structure (defined in 'time.h') with the values derived
   from a FITS format date-string and return the (optional) sub-second
   information as a double.

   The basic FITS string is defined under the 'DATE' keyword in the FITS
   standard. For the more complete format which includes timezones, see the
   W3 standard: https://www.w3.org/TR/NOTE-datetime */
static char *
gal_units_date_to_struct_tm(char *datestr, struct tm *tp)
{
  char *C, *cc, *c=NULL, *cf, *subsec=NULL, *nosubsec=datestr;
  int hasT=0, hassq=0, usesdash=0, usesslash=0, hasZ=0, hasoffset=0;

  /* Initialize the 'tm' structure to all-zero elements. In particular, The
     FITS standard times are written in UTC, so, the time zone ('tm_zone'
     element, which specifies number of seconds to shift for the time zone)
     has to be zero. The day-light saving flag ('isdst' element) also has
     to be set to zero. */
  tp->tm_sec=tp->tm_min=tp->tm_hour=tp->tm_mday=tp->tm_mon=tp->tm_year=0;
  tp->tm_wday=tp->tm_yday=tp->tm_isdst=tp->tm_gmtoff=0;
  tp->tm_zone=NULL;

  /* According to the FITS standard the 'T' in the middle of the date and
     time of day is optional (the time is not mandatory). */
  cf=(c=datestr)+strlen(datestr);
  do
    switch(*c)
      {
      case  'T': hasT=1;      break; /* With 'T' HH:MM:SS are defined.    */
      case  '/': usesslash=1; break; /* Day definition(old): DD/MM/YY.    */
      case '\'': hassq=1;     break; /* Wholly Wrapped in a single-quote. */
      case  '+': hasoffset=1; break; /* An offset, see the case for '-'.  */
      case  'Z': hasZ=1;      break; /* When ends in 'Z', means UTC. See  */
                                   /* https://www.w3.org/TR/NOTE-datetime */

      /* In case we have sub-seconds in the string, we need to remove it
         because 'strptime' doesn't recognize sub-second resolution. */
      case '.':

        /* Allocate space (by copying the remaining full string and adding
           a '\0' where necessary. */
        gal_checkset_allocate_copy(c, &subsec);
        gal_checkset_allocate_copy(datestr, &nosubsec);

        /* Parse the sub-seconds part and end the 'subsec' string when the
           sub-seconds digits finishs (there may be characters after it. */
        for(C=subsec+1;*C!='\0';C++) if(!isdigit(*C)) {*C='\0'; break;}

        /* Copy everything after the sub-seconds into the 'nosubsec'
           string. */
        cc=nosubsec+(c-datestr);
        for(C=c+(C-subsec); *C!='\0'; ++C) *cc++=*C;
        *cc='\0';
        break;

      /* A dash (hyphen) can be interpreted in two ways

           - Before any possible 'T': if we see a dash, it shows the date
             in the YYYY-MM-DD format.

           - After a possible 'T': according to the ISO 8601 standard (more
             general than FITS: https://en.wikipedia.org/wiki/ISO_8601), an
             offset (timezone), it will be given at the end of the string,
             starting with a '+HH:MM' or '-HH:MM'. So if we see a hyphen
             after the 'T', it is actually an offset/negative value.
              */
      case '-':
        if(hasT) /* Already passed the 'T': this is a negative offset. */
          hasoffset=1;
        else /* Not passed a 'T' yet: this is the date. */
          usesdash=1;
        break;
      }
  while(++c<cf);

  /* Sanity checks: */
  if(usesdash && usesslash)
    error(EXIT_FAILURE, 0, "'%s' is not readable because the date "
          "includes a hyphen (dash: '-') and a slash ('/'). When "
          "a hyphen is used, the date is expected to be in ISO8601 "
          "format (or YYYY-MM-DD), but when a dash is used it is "
          "interpreted as DD/MM/YY", datestr);
  if(hasZ && hasoffset)
    error(EXIT_FAILURE, 0, "'%s' is not readable because it contains "
          "both the 'Z' character (implying no HH:MM offset) and "
          "an offset, it should be either 'YYYY-MM-DDThh:mm:ssZ' (no "
          "offset) or 'YYYY-MM-DDThh:mm:ss+OH:OM' (with an offset of "
          "'OH' hours and 'OM' minutes", datestr);

  /* Convert this date into seconds since 1970/01/01, 00:00:00. */
  c = ( hasoffset
        ? ( usesdash
            ? ( hasT
                ? strptime(nosubsec, hassq?"'%FT%T%z'":"%FT%T%z", tp)
                : strptime(nosubsec, hassq?"'%F%z'":"%F%z", tp) )
            : ( hasT
                ? strptime(nosubsec,
                           hassq?"'%d/%m/%yT%T%z'":"%d/%m/%yT%T%z", tp)
                : strptime(nosubsec, hassq?"'%d/%m/%y%z'":"%d/%m/%y%z",
                           tp) )
            )
        : ( usesdash
            ? ( hasT
                ? ( hasZ
                    ? strptime(nosubsec, hassq?"'%FT%TZ'":"%FT%TZ", tp)
                    : strptime(nosubsec, hassq?"'%FT%T'":"%FT%T", tp) )
                : strptime(nosubsec, hassq?"'%F'":"%F", tp) )
            : ( hasT
                ? ( hasZ
                    ? strptime(nosubsec,
                               hassq?"'%d/%m/%yT%TZ'":"%d/%m/%yT%TZ", tp)
                    : strptime(nosubsec,
                               hassq?"'%d/%m/%yT%T'":"%d/%m/%yT%T", tp))
                : strptime(nosubsec, hassq?"'%d/%m/%y'"   :"%d/%m/%y", tp)
                )
            )
        );

  /* The value might have sub-seconds. In that case, 'c' will point to a
     '.' and we'll have to parse it as double. */
  if( c==NULL || (*c!='.' && *c!='\0') )
    error(EXIT_FAILURE, 0, "'%s' isn't in the correct FITS or ISO 8601 "
          "date format. Please run 'info gnuastro date' for more on "
          "the expected input date format", datestr);

  /* If the subseconds were removed (and a new string was allocated), free
     that extra new string. */
  if(nosubsec!=datestr) free(nosubsec);

  /* Return the subsecond value. */
  return subsec;
}





/* Convert the FITS standard or ISO 8601 date format (as a string, already
   read from the keywords) into number of seconds since 1970/01/01,
   00:00:00. Very useful to avoid calendar issues like number of days in a
   different months or leap years and etc. The remainder of the format
   string (sub-seconds) will be put into the two pointer arguments:
   'subsec' is in double-precision floating point format but it will only
   get a value when 'subsecstr!=NULL'. */
int64_t
gal_units_date_to_unix_seconds(char *datestr, char **subsecstr,
                               double *subsec)
{
  time_t t;
  char *tmp;
  struct tm tp;
  void *subsecptr=subsec;
  int64_t seconds, gmtoff;

  /* If the string is blank, return a blank value. */
  if( strcmp(datestr, GAL_BLANK_STRING)==0 )
    {
      if(subsec) *subsec=NAN;
      if(subsecstr) *subsecstr=NULL;
      return GAL_BLANK_INT64;
    }

  /* Fill in the 'tp' elements with values read from the string. */
  tmp=gal_units_date_to_struct_tm(datestr, &tp);

  /* If the user wanted a possible sub-second string/value, then
     'subsecstr!=NULL'. */
  if(subsecstr)
    {
      /* Set the output pointer. Note that it may be NULL if there was no
         sub-second string, but that is fine (and desired because the user
         can use this to check if there was a sub-string or not). */
      *subsecstr=tmp;

      /* If there was a sub-second string, then also read it as a
         double-precision floating point. */
      if(tmp)
        {
          if(subsec)
            if( gal_type_from_string(&subsecptr, tmp, GAL_TYPE_FLOAT64) )
              error(EXIT_FAILURE, 0, "%s: the sub-second portion of '%s' "
                    "(or '%s') couldn't be read as a number", __func__,
                    datestr, tmp);
        }
      else { if(subsec) *subsec=NAN; }
    }

  /* Convert the contents of the 'tm' structure to 'time_t' (a positive
     integer). Note that by design (as described in its GNU C Library
     manual), 'timegm' (like 'mktime') will ignore 'tm_gmtoff' (offset from
     GMT/UTC) and set it to zero. However, this offset is part of the ISO
     8601 convention (that we accept here), so we need to read it before
     calling 'timegm' and apply (subtract) it afterwards. */
  gmtoff=tp.tm_gmtoff;
  t=timegm(&tp); /* should be after extracting 'gmtoff'.*/
  seconds = (t == (time_t)(-1)) ? GAL_BLANK_INT64 : t-gmtoff;
  return seconds;
}





/* Convert Unix-time (possibly with sub-seconds) to a UTC-based date
   string. */
char *
gal_units_unix_seconds_to_date(int64_t unixsec, double subsec,
                               int subsecdigits)
{
  size_t nchr;
  time_t t=unixsec;
  struct tm tp, *tpo;
  char *tmp, *tout, strfmt[5], *out=NULL;

  /* Sanity check. */
  if(subsec<0.0f)
    error(EXIT_FAILURE, 0, "%s: the value of the 'subsec' argument "
          "should be larger or equal to 0.0. However, its value is "
          "%e", __func__, subsec);

  /* Break-down the Unix-seconds into YYYY-MM-DD and hh:mm:ss
     components. We are using 'gmtime' to have an output in UTC, and its
     '_r' variant to enable parallel processing. */
  tpo=gmtime_r(&t, &tp);
  if( tpo != &tp )
    error(EXIT_FAILURE, 0, "%s: the unix-time '%ld' could not "
          "be broken down into a calendar time", __func__, unixsec);

  /* Allocate the output string and fill it up. Note that the desired
     output format of 'YYYY-MM-DDThh:mm:ss' has 20 characters (including
     the string termination character. But for B.C. dates, a '-' will be
     added on the year and for dates too far into the future, the year can
     be longer, so to avoid strange outputs, we'll allow 29 characaters in
     the printable string. */
  out=gal_pointer_allocate(GAL_TYPE_UINT8, 30, 0, __func__, "out");
  nchr=strftime(out, 30, "%FT%T", &tp);
  if( nchr > 29 ) /* Number of chars ('nchr') does not include '/0'. */
    error(EXIT_FAILURE, 0, "%s: the broken-down time could not "
          "printed in the desired format for the input Unix-second "
          "'%ld'", __func__, unixsec);

  /* If a non-zero sub-second string is given, add it to the output. */
  if(subsecdigits && subsec!=0.0f)
    {
      sprintf(strfmt, "%%.%df", subsecdigits);
      asprintf(&tmp, strfmt, subsec);
      tout=gal_checkset_malloc_cat(out, tmp+1);
      free(out);
      out=tout;
    }

  /* For a check.
  printf("%s: %s\n", __func__, out); exit(0);
  //*/

  /* Return the output string. */
  return out;
}





/* Unix-seconds (US) to Julian days (JD; epoch: 12:00 January 1, 4713 BC)

       US = (JD - 2440587.5) × 86400
       JD = US/86400 + 2440587.5

   If there are sub-seconds given, that is also divided by 86400 and added
   to the JD of the equation above. For more, see
   https://en.wikipedia.org/wiki/Julian_day */
double
gal_units_unix_seconds_to_jd(int64_t unixsec, double subsec)
{
  double jd=((double)unixsec)/86400+2440587.5;
  return isnan(subsec) ? jd : jd+subsec/86400;
}

int64_t
gal_units_jd_to_unix_seconds(double jd, double *subsec)
{
  double d=(jd-2440587.5)*86400;
  int64_t i=d;
  *subsec=d-i;
  return i;
}





/* Unix-seconds to Reduced Julian date. Epoch: 12:00 November 16, 1858. See
   https://en.wikipedia.org/wiki/Julian_day . The conversion */
double
gal_units_jd_to_rjd(double jd)
{
  return jd-2400000;
}

int64_t
gal_units_rjd_to_unix_seconds(double rjd, double *subsec)
{
  double d=rjd*86400;
  int64_t i=d;
  *subsec=d-i;
  return i - 3506760000;
}





/* Unix-seconds to Modified Julian date. Epoch: 0:00 November 17, 1858. See
   https://en.wikipedia.org/wiki/Julian_day */
double
gal_units_jd_to_mjd(double jd)
{
  return jd-2400000.5;
}

int64_t
gal_units_mjd_to_unix_seconds(double mjd, double *subsec)
{
  double d=mjd*86400;
  int64_t i=d;
  *subsec=d-i;
  return i - 3506716800;
}





/* Unix-seconds to Truncated Julian date. Epoch: 0:00 May 24, 1968. See
   https://en.wikipedia.org/wiki/Julian_day */
double
gal_units_jd_to_tjd(double jd)
{
  return jd-2440000.5;
}

int64_t
gal_units_tjd_to_unix_seconds(double tjd, double *subsec)
{
  double d=tjd*86400;
  int64_t i=d;
  *subsec=d-i;
  return i - 50716800;
}





/* Unix-seconds to Dublin Julian date. Epoch: 12:00 December 31, 1899. See
   https://en.wikipedia.org/wiki/Julian_day */
double
gal_units_jd_to_djd(double jd)
{
  return jd-2415020;
}

int64_t
gal_units_djd_to_unix_seconds(double djd, double *subsec)
{
  double d=djd*86400;
  int64_t i=d;
  *subsec=d-i;
  return i - 2209032000;
}





/* Unix-seconds to CNES Julian date. Epoch: 0:00 January 1, 1950. See
   https://en.wikipedia.org/wiki/Julian_day */
double
gal_units_jd_to_cjd(double jd)
{
  return jd-2433282.5;
}

int64_t
gal_units_cjd_to_unix_seconds(double cjd, double *subsec)
{
  double d=cjd*86400;
  int64_t i=d;
  *subsec=d-i;
  return i - 631152000;
}





/* Unix-seconds to CCSDS Julian date. Epoch: 0:00 January 1, 1958. See
   https://en.wikipedia.org/wiki/Julian_day */
double
gal_units_jd_to_ccjd(double jd)
{
  return jd-2436204.5;
}

int64_t
gal_units_ccjd_to_unix_seconds(double ccjd, double *subsec)
{
  double d=ccjd*86400;
  int64_t i=d;
  *subsec=d-i;
  return i - 378691200;
}





/* Unix-seconds to Modified Julian date 2000. Epoch: 0:00 January 1,
   2000. See https://en.wikipedia.org/wiki/Julian_day */
double
gal_units_jd_to_mjd2000(double jd)
{
  return jd-2451544.5;
}

int64_t
gal_units_mjd2000_to_unix_seconds(double mjd2000, double *subsec)
{
  double d=mjd2000*86400;
  int64_t i=d;
  *subsec=d-i;
  return i + 946684800;
}





/* Unix-seconds to Lilian date. Similar to JD, but with day 1 being October
   15, 1582 (first day of the Gregorian calendar: days have been properly
   counted since then): see https://en.wikipedia.org/wiki/Julian_day */
double
gal_units_jd_to_lilian(double jd)
{
  return jd-2299159.5;
}

int64_t
gal_units_lilian_to_unix_seconds(double lilian, double *subsec)
{
  double d=lilian*86400;
  int64_t i=d;
  *subsec=d-i;
  return i - 12219379200;
}





/* Unix-seconds to Rata Die. Day 1 = January 1, 0001. See
   https://en.wikipedia.org/wiki/Julian_day */
double
gal_units_jd_to_rata(double jd)
{
  return jd-1721424.5;
}

int64_t
gal_units_rata_to_unix_seconds(double rata, double *subsec)
{
  double d=rata*86400;
  int64_t i=d;
  *subsec=d-i;
  return i - 62135683200;
}





/* Unix-seconds to Mars Sol Date. 12:00 December 29, 1873. See
   https://en.wikipedia.org/wiki/Timekeeping_on_Mars#Mars_Sol_Date and
   https://www.giss.nasa.gov/tools/mars24/help/algorithm.html (equation for
   MST: three constants merged into one). */
double
gal_units_jd_to_mars_sol(double jd)
{
  return (jd-2.49514695148767e+06)/1.0274912517;
}

int64_t
gal_units_mars_sol_to_unix_seconds(double sol, double *subsec)
{
  double d=sol*1.0274912517*86400;
  int64_t i=d;
  *subsec=d-i;
  return i + 4713936608;
}
