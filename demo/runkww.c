/* runkww.c
 *
 * Copyright (C) 2009 Joachim Wuttke
 *
 * Licence: GNU General Public License, version 3 or later
 *
 * Author:
 *   Joachim Wuttke, Forschungszentrum Jülich, Germany <j.wuttke@fz-juelich.de>
 *
 * Purpose:
 *   Command-line interface to libkww:
 *      return one function value, along with some info on the used algorithm
 */

#include <stdio.h>
#include <stdlib.h>
#include "kww.h"
#include "kww_lowlevel.h"

KWW_IMPORT extern int kww_algorithm;
KWW_IMPORT extern int kww_num_of_terms;
KWW_IMPORT extern int kww_debug;

int main( int argc, char **argv )
{
    char dir, alg;
    double w, b, ret;

    if( argc!=6 ){
        fprintf( stderr,  "usage:\n" );
        fprintf( stderr,  "   runkww <deb> c|s|p a|l|m|h <b> <w>\n" );
        fprintf( stderr,  "with arguments:\n" );
        fprintf( stderr,  "   <deb>: integer debug code\n" );
        fprintf( stderr,  "   flag1: c: cos transform\n" );
        fprintf( stderr,  "          s: sin transform\n" );
        fprintf( stderr,  "          p: primitive of cos transform\n" );
        fprintf( stderr,  "   flag2: a: automatic determination of algorithm\n" );
        fprintf( stderr,  "          l: low-omega series expansion\n" );
        fprintf( stderr,  "          m: mid-omega numeric integration\n" );
        fprintf( stderr,  "          h: hig-omega series expansion\n" );
        fprintf( stderr,  "   b: stretching exponent (between 0.1 and 2)\n" );
        fprintf( stderr,  "   w: omega\n" );
        fprintf( stderr,  "output:\n" );
        fprintf( stderr,  "   value1: the function value of kwws or kwwc\n" );
        fprintf( stderr,  "   value2: algorithm used (1=l, 2=m, 3=h)\n" );
        fprintf( stderr,  "   value3: the number of summed terms\n" );
        exit(-1);
    }

    kww_debug = atoi( argv[1] );

    dir = argv[2][0];
    if( dir!='c' && dir!='s' && dir!='p' ){
        fprintf( stderr, " choose transform 'c' or 's' or 'p'\n" );
        exit(-1);
    }
    alg = argv[3][0];
    if( alg!='a' && alg!='l' && alg!='m' && alg!='h' ){
        fprintf( stderr, " choose algorithm 'a' or 'l' or 'm' or 'h'\n" );
        exit(-1);
    }
    b     = atof( argv[4] );
    w     = atof( argv[5] );

    ret = -111;
    if     ( alg=='a' ) {
        if      ( dir=='c' )
            ret = kwwc( w, b );
        else if ( dir=='s' )
            ret = kwws( w, b );
        else if ( dir=='p' )
            ret = kwwp( w, b );
    } else if( alg=='l' ) {
        if      ( dir=='c' )
            ret = kwwc_low( w, b );
        else if ( dir=='s' )
            ret = kwws_low( w, b );
        else if ( dir=='p' )
            ret = kwwp_low( w, b );
    } else if( alg=='m' ) {
        if      ( dir=='c' )
            ret = kwwc_mid( w, b );
        else if ( dir=='s' )
            ret = kwws_mid( w, b );
        else if ( dir=='p' )
            ret = kwwp_mid( w, b );
    } else if( alg=='h' ) {
        if      ( dir=='c' )
            ret = kwwc_hig( w, b );
        else if ( dir=='s' )
            ret = kwws_hig( w, b );
        else if ( dir=='p' )
            ret = kwwp_hig( w, b );
    } else {
        fprintf( stderr, "invalid alg flag\n" );
        exit( -1 );
    }
    printf( "%25.19g %1i %6i\n", ret, kww_algorithm, kww_num_of_terms );
    return 0;
}
