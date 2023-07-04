/* kww_counterms.c
 *
 * Copyright (C) 2009 Joachim Wuttke
 *
 * Licence: GNU General Public License, version 3 or later
 *
 * Author:
 *   Joachim Wuttke, Forschungszentrum Jülich, Germany <j.wuttke@fz-juelich.de>
 *
 * Purpose:
 *   Count number of terms in the trapezoid sum in libkww.
 *   Used to determine the optimum double-exponential transformation.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kww.h"
#include "kww_lowlevel.h"

KWW_IMPORT extern int kww_num_of_terms;

double kwwc_lim_low( const double beta );
double kwwc_lim_hig( const double beta );
double kwws_lim_low( const double beta );
double kwws_lim_hig( const double beta );

int main( int argc, char **argv )
{
    int nb, nw, j, k, i;
    long int s0, s1;
    double b, r, wl, wh, w;

    if( argc!=3 && argc!=4 ){
        fprintf( stderr,  "usage:\n" );
        fprintf( stderr,  "   %s <nb> <nw> [<b>]\n", argv[0] );
        fprintf( stderr,  "with arguments:\n" );
        fprintf( stderr,  "   <nb>: number of different beta's\n" );
        fprintf( stderr,  "   <nw>: number of different omega's\n" );
        fprintf( stderr,  "   <b>:  value of beta if nb=1\n" );
        fprintf( stderr,  "output:\n" );
        fprintf( stderr,  "   number of terms w\n" );
        exit(-1);
    }

    nb = atoi(argv[1]);
    nw = atoi(argv[2]);

    s0 = 0;
    for( j=0; j<nb; ++j ){
        if( nb== 1 )
            b = atof(argv[3]);
        else
            b = 0.1 + 1.8 * j/(nb-1.0);
        s1 = 0;
        for( k=0; k<2; ++k ){
            if( k==0 ){
                wl = kwwc_lim_low( b );
                wh = kwwc_lim_hig( b );
            } else {
                wl = kwws_lim_low( b );
                wh = kwws_lim_hig( b );
            }
            wl /= 1.02;
            wh *= 1.02;
            for( i=0; i<nw; ++i ){
                w = wl * pow(wh/wl, i/(nw-1.0));
                if( k==0 )
                    r = kwwc_mid( w, b );
                else
                    r = kwws_mid( w, b );
                if( r<0 ){
                    fprintf( stderr,
                             "integration %1i failed %25.18g %25.18g -> %g\n",
                             k, b, w, r );
                    exit(1);
                }
                s1 += kww_num_of_terms;
            }
        }
        printf( "%15.9g %15.9g\n", b, (double)s1/nw );
        s0 += s1;
    }
    printf( "total: %15.9g\n", (double)s0/nb/nw );
    return 0;
}
