/* kww_echolims.c
 * 
 * Copyright (C) 2009 Joachim Wuttke
 * 
 * Licence: GNU General Public License, version 3 or later
 *
 * Author:
 *   Joachim Wuttke, Forschungszentrum Jülich, Germany <j.wuttke@fz-juelich.de>
 *
 * Purpose:
 *   Print hard-coded limits of series expansion domains from kww.c
 *   Needed to produce the black lines in figure kww-fig-lims.
 *   Also print KWW values.
 *   Needed to determine S_min, to fix maximum allowed absolute error.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kww.h"

double kwwc_lim_low( const double beta );
double kwwc_lim_hig( const double beta );
double kwws_lim_low( const double beta );
double kwws_lim_hig( const double beta );
double kwwp_lim_low( const double beta );
double kwwp_lim_hig( const double beta );

extern int kww_algorithm, kww_num_of_terms;

int main( int argc, char **argv )
{
    char dir, lim;
    int nb, i;
    double b,w,s,e=1.2;

    if( argc!=4 && argc!=5 ){
        fprintf( stderr,  "usage:\n" );
        fprintf( stderr,  "   kww_echolims c|s|p l|h <nb> [<b>]\n" );
        fprintf( stderr,  "with arguments:\n" );
        fprintf( stderr,  "   <nb>:   number of beta values\n" );
        fprintf( stderr,  "   <b>:    fixed beta if nb=1\n" );
        fprintf( stderr,  "output lines contain:\n" );
        fprintf( stderr,  "   beta, w, S~(w)\n" );
        exit(-1);
    }

    dir = argv[1][0];
    lim = argv[2][0];

    nb = atoi(argv[3]);

    for( i=0; i<nb; ++i ){
        if( nb==1 )
            b = atof( argv[4] );
        else
            b = 0.1 * pow(1.999/0.1, i/(nb-1.0));
        if        ( dir=='c' ) {
            w = lim == 'l' ? kwwc_lim_low(b) : kwwc_lim_hig(b);
            s = kwwc(w/e,b);
        } else if ( dir=='s' ) {
            w = lim == 'l' ? kwws_lim_low(b) : kwws_lim_hig(b);
            s = kwws(w/e,b);
        } else if ( dir=='p' ) {
            w = lim == 'l' ? kwwp_lim_low(b) : kwwp_lim_hig(b);
            s = kwwp(w/e,b);
        }
        printf( "%12.5g %12.5g %12.5g\n", b, w, w/3.14*s );
    }
    return 0;
}
