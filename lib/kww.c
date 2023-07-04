/* kww.c:
 *   Calculation of the Kohlrausch-Williams-Watts spectrum, i.e.
 *   Laplace-Fourier transform of the stretched exponential function exp(-t^b).
 *   Frequently used to describe relaxation in disordered systems.
 *
 * Copyright:
 *   (C) 2009, 2012 Joachim Wuttke
 *
 * Licence:
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published
 *   by the Free Software Foundation; either version 3 of the License, or
 *   (at your option) any later version. Alternative licenses can be
 *   obtained through written agreement from the author.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but without any warranty; without even the implied warranty of
 *   merchantability or fitness for a particular purpose.
 *   See the GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * Author:
 *   Joachim Wuttke
 *   Forschungszentrum Jülich, Germany
 *   j.wuttke@fz-juelich.de
 *
 * Website:
 *   https://jugit.fz-juelich.de/mlz/kww
 *
 * Reference:
 *   Wuttke, Algorithms 5, 604-628 (2012), doi:10.3390/a5040604
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include "kww.h"
#include "kww_lowlevel.h"

#ifdef __MINGW32__
#define printf __mingw_printf
#endif

#define PI           3.14159265358979323846L  /* pi */
#define SQR(x) ((x)*(x))

/*****************************************************************************/
/*  Approximate limits of asymptotic regimes                                 */
/*****************************************************************************/

double kwwc_lim_low( const double b )
{
    if ( b>1.024 )
        return -0.8774954*b + 3.5873*b*b -2.083*b*b*b +0.3796*b*b*b*b;
    else
        return exp( -0.02194/b/b -4.130/b +2.966189 +0.030104*b +1.062*b*b );
}

double kwws_lim_low( const double b )
{
    if ( b>1.024 )
        return -1.68725*b + 4.8108*b*b -2.561*b*b*b +0.442*b*b*b*b;
    else
        return exp( -0.03208/b/b -4.314/b +3.516200 -0.50287*b +1.240*b*b );
}

double kwwp_lim_low( const double b )
{
    if ( b>1.085 )
        return -10.49909 +19.23618*b -9.234064*b*b +1.553016*b*b*b;
    else
        return exp( -0.02259971/b/b -4.099837/b +3.100445
                    -0.1838126*b +1.118149*b*b );
}


double kwwc_lim_hig( const double b )
{
    if ( b<0.82 )
        return exp( 0.006923209/b/b -1.321692/b -1.44582
                    +2.516339*b +0.2973773*b*b );
    else
        return exp( -0.746496154631 +6.057558*(b-.82) -3.41052*SQR(b-.82)
                    +0.7932314*pow(b-.82,3) );
}

double kwws_lim_hig( const double b )
{
    if ( b<0.82 )
        return exp( 0.07847516/b/b -2.585876/b +4.999414
                    -8.460926*b +6.289183*b*b );
    else
        return exp( -0.962597724393 +5.818057*(b-.82) -3.026212*SQR(b-.82)
                    +0.5485754*pow(b-.82,3) );
}

double kwwp_lim_hig( const double b )
{
    if ( b<0.82 )
        return exp( 0.003809101/b/b -1.955504/b -1.938468
                    +5.893199*b -2.197289*b*b );
    else
        return exp( -0.962597724393 +7.074977*(b-.82) -5.231151*SQR(b-.82)
                    +1.717068*pow(b-.82,3) );
}


/*****************************************************************************/
/*  High-level wrapper functions                                             */
/*****************************************************************************/

/* \int_0^\infty dt cos(w*t) exp(-t^beta) */
double kwwc( const double w_in, const double beta )
{
    double w, res;
    /* check input data */
    if ( beta<0.1 ) {
        fprintf( stderr, "kww: beta smaller than 0.1\n" );
        exit( EDOM );
    }
    if ( beta>2.0 ) {
        fprintf( stderr, "kww: beta larger than 2.0\n" );
        exit( EDOM );
    }
    /* it's an even function; the value at w=0 is well known */
    if ( w_in==0 )
        return tgamma(1.0/beta)/beta;
    w = fabs( w_in );
    /* special case: Gaussian for b=2 */
    if ( beta==2 )
        return sqrt(PI)/2*exp(-SQR((double)w)/4);
    /* try series expansion */
    if        ( w<kwwc_lim_low( beta ) ) {
        Xdouble s = kwwc_low( w, beta );
        if ( s>0 )
            return s;
        res = s;
    } else if ( w>kwwc_lim_hig( beta ) ) {
        Xdouble s = kwwc_hig( w, beta );
        if ( s>0 )
            return s;
        res = s;
    }
    /* fall back to numeric integration */
    res = kwwc_mid( w, beta );
    if ( res<0 ) {
        if( beta>1.9 )
            return 0; // must be tested by the user
        // otherwise it ought to be considered a bug
        fprintf( stderr, "kwwc: numeric integration failed for"
                 " omega=%25.18g, beta=%25.18g; error code %g\n",
                 w_in, beta, res );
        exit( ENOSYS );
    }
    return res;
}

/* \int_0^\infty dt sin(w*t) exp(-t^beta) */
double kwws( const double w_in, const double beta )
{
    double w, res;
    int sign_out;
    /* check input data */
    if ( beta<0.1 ) {
        fprintf( stderr, "kww: beta smaller than 0.1\n" );
        exit( EDOM );
    }
    if ( beta>2.0 ) {
        fprintf( stderr, "kww: beta larger than 2.0\n" );
        exit( EDOM );
    }
    /* it's an odd function */
    if ( w_in==0 )
        return 0;
    if ( w_in<0 ) {
        w = - w_in;
        sign_out = -1;
    } else {
        w = w_in;
        sign_out = 1;
    }
    /* try series expansion */
    if        ( w<kwws_lim_low( beta ) ) {
        Xdouble s = kwws_low( w, beta );
        if ( s>0 )
            return sign_out*s;
        res = s;
    } else if ( w>kwws_lim_hig( beta ) ) {
        Xdouble s = kwws_hig( w, beta );
        if ( s>0 )
            return sign_out*s;
        res = s;
    }
    /* fall back to numeric integration */
    res = kwws_mid( w, beta );
    if ( res<0 ) {
        fprintf( stderr, "kwws: numeric integration failed for"
                 " omega=%25.18g, beta=%25.18g; error code %g\n",
                 w_in, beta, res );
        exit( ENOSYS );
    }
    return sign_out*res;
}

/* \int_0^w dw' \int_0^\infty dt cos(w'*t) exp(-t^beta) */
double kwwp( const double w_in, const double beta )
{
    double w, res;
    int sign_out;
    /* check input data */
    if ( beta<0.1 ) {
        fprintf( stderr, "kww: beta smaller than 0.1\n" );
        exit( EDOM );
    }
    if ( beta>2.0 ) {
        fprintf( stderr, "kww: beta larger than 2.0\n" );
        exit( EDOM );
    }
    /* it's an odd function */
    if ( w_in==0 )
        return 0;
    if ( w_in<0 ) {
        w = - w_in;
        sign_out = -1;
    } else {
        w = w_in;
        sign_out = 1;
    }
    /* try series expansions */
    if        ( w<kwwp_lim_low( beta ) ) {
        Xdouble s = kwwp_low( w, beta );
        if ( s>0 )
            return sign_out*s;
        res = s;
    } else if ( w>kwwp_lim_hig( beta ) ) {
        Xdouble s = kwwp_hig( w, beta );
        if ( s>0 )
            return sign_out*s;
        res = s;
    }
    /* fall back to numeric integration */
    res = kwwp_mid( w, beta );
    if ( res<0 ) {
        fprintf( stderr, "kwwp: numeric integration failed for"
                 " omega=%25.18g, beta=%25.18g; error code %g\n",
                 w_in, beta, res );
        exit( ENOSYS );
    }
    return sign_out*res;
}
