/* lowlevel.c:
 *   Computation either by series expansion (low/hig) or by integration (mid).
 *
 * Copyright:
 *   (C) 2009, 2012, 2023 Joachim Wuttke
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
#include <float.h>
#include <errno.h>
#include "kww.h"
#include "kww_lowlevel.h"

#ifdef __MINGW32__
#define printf __mingw_printf
#endif

#define PI           3.14159265358979323846L  /* pi */
#define PI_2         1.57079632679489661923L  /* pi/2 */
#define SQR(x) ((x)*(x))

// for external analysis:
KWW_EXPORT int kww_algorithm;
KWW_EXPORT int kww_num_of_terms;
KWW_EXPORT int kww_debug=0;

/*****************************************************************************/
/*  Numeric precision and maximum number of terms                            */
/*****************************************************************************/

const double kww_delta=2.2e-16, kww_eps=5.5e-20;
const int max_terms=200;

/*****************************************************************************/
/*  Low-level implementation: series expansion for low frequencies           */
/*****************************************************************************/

Xdouble kww_low( const double w, const double beta,
                const int kappa, const int mu )
{
    int kk;               // this is 2*k+kappa
    int isig=1;           // alternating sign
    Xdouble S=0;      // summed series
    Xdouble T=0;      // sum of absolute values
    Xdouble u;        // precomputed common factors
    Xdouble u_next=0; // - next value [initialized to avoid warning]
    Xdouble gl;       // local variable

    // set diagnostic variable
    kww_algorithm = 1;

    // check input
    if ( beta<0.1 || beta>2.0 ) {
        fprintf( stderr, "invalid call to kww_low: beta out of range\n" );
        exit( EDOM );
    }
    if ( w<=0 ) {
        fprintf( stderr, "invalid call to kww_low: w out of range\n" );
        exit( EDOM );
    }

    // sum the expansion
    kk = kappa;
    for ( int i=0; i<max_terms; ++i ) {
        kww_num_of_terms = i;
        // t_n must be computed in advance
        u = u_next;
        // use log gamma instead of gamma to avoid overflow
        gl = lgammaX((Xdouble)(kk+1)/(Xdouble)beta)-
            lgammaX((Xdouble)kk+1)+(kk+mu)*logX((Xdouble)w);
        if ( gl>DBL_MAX_EXP/2 )
            return -3; // gamma function overflow
        u_next = expX( gl );
        if( mu )
            u_next /= (kk+1);
        kk += 2;
        if( !i )
            continue;
        // now we use t_{n-1} to compute S_n
        S += isig*u;
        T += u;
        // termination criteria
        if ( kww_eps*T+u_next <= kww_delta*S )
            return S / beta; // reached required precision
        else if ( kww_eps*T >= kww_delta*S )
            return -6; // too much cancellation
        else if ( beta<1 && u_next>u )
            return -5; // asymptotic expansion diverges too early
        else if ( S<DBL_MIN )
            return -7; // underflow
        isig = -isig;
    }
    return -9; // too many terms
}

Xdouble kwwc_low( const double w, const double beta )
{
    return kww_low( w, beta, 0, 0 );
}

Xdouble kwws_low( const double w, const double beta )
{
    return kww_low( w, beta, 1, 0 );
}

Xdouble kwwp_low( const double w, const double beta )
{
    return kww_low( w, beta, 0, 1 );
}

/*****************************************************************************/
/*  Low-level implementation: series expansion for high frequencies          */
/*****************************************************************************/

Xdouble kww_hig( const double w, const double beta,
                const int kappa, const int mu )
{
    int k;           // in computation of A_k w^k
    int isig=1;      // alternating sign
    int alternating; // has factor (-)^k
    Xdouble b;        // either beta or 2-beta
    Xdouble sinphi;   // to compute r from u
    Xdouble truncfac; // for termination criterion
    Xdouble rfac;     // for computation of remainder
    Xdouble S=0;      // summed series
    Xdouble Sabs;     // absolute value thereof
    Xdouble T=0;      // sum of absolute values
    Xdouble u;        // precomputed common factors
    Xdouble u_next=0; // - next value [initialized to avoid warning]
    Xdouble s;        // full term (with trigonometric factor)
    Xdouble x, gl;    // local variables

    // set diagnostic variable
    kww_algorithm = 3;

    // check input
    if ( beta<0.1 || beta>2.0 ) {
        fprintf( stderr, "invalid call to kww_hig: beta out of range\n" );
        exit( EDOM );
    }
    if ( w<=0 ) {
        fprintf( stderr, "invalid call to kww_hig: w out of range\n" );
        exit( EDOM );
    }

    // set some beta-dependent constants
    if ( beta<1 ) {
        b = beta;
        alternating = 1;
        sinphi = 1;
        truncfac = 1;
    } else {
        b = 2.0-beta;
        alternating = 0;
        sinphi = sinX( PI_2/(Xdouble)beta );
        truncfac = powX( sinphi, -(Xdouble)beta );
    }
    rfac = 1/sinphi;

    if( kww_debug & 2 ) {
        printf( "sinphi %20.14Le truncfac %20.14Le\n", sinphi, truncfac );
    }
    if( kww_debug & 1 )
        printf( "%3s %20s %20s %12s %12s %12s %12s %12s %12s %12s\n",
                "k+2", "S=sum(s)", "T=sum(|s|)",
                "s", "s/u", "u(k)", "u(k+1)", "rfac",
                "eps*T+u(k+1)*r", "delta*|S|" );

    // sum the expansion
    k=1-kappa;
    if( k )
        rfac *= truncfac;
    for ( int i=0; i<max_terms; ++i ) {
        kww_num_of_terms = i;
        // t_n must be computed in advance
        u = u_next;
        x = k*(Xdouble)beta+1;
        // use log gamma instead of gamma to avoid overflow
        gl = lgammaX(x)-lgammaX((Xdouble)k+1)+(mu-x)*logX((Xdouble)w);
        if ( gl>DBL_MAX_EXP/2 )
            return -3; // gamma function overflow
        u_next = expX( gl );
        if( mu )
            u_next /= (k*beta);
        ++k;
        if( !i )
            continue;
        // now we use t_{n-1} to compute S_n (k is even 2 ahead)
        s = u * isig * ( kappa ? cosX(PI_2*(k-2)*b) : sinX(PI_2*(k-2)*b) );
        S += s;
        Sabs = fabsX(S);
        T += fabsX(s);
        rfac *= truncfac; // sin(phi)^(-1-k*beta)
        if( kww_debug & 1 )
            printf( "%3i %20.13Le %20.13Le %12.5Le %12.5Le %12.5Le %12.5Le"
                    " %12.5Le %12.5Le %12.5Le\n",
                    k, S, T, s, s/u, u, u_next, rfac,
                    kww_eps*T+u_next*rfac, kww_delta*Sabs );
        // termination criteria
        if ( kww_eps*T+u_next*rfac <= kww_delta*Sabs )
            return S; // reached required precision
        else if ( beta>1 && u_next*truncfac>u )
            return -5; // asymptotic expansion diverges too early
        else if ( Sabs<DBL_MIN )
            return -7; // underflow
        if ( alternating )
            isig = -isig;
    }
    return -9; // not converged
}

Xdouble kwwc_hig( const double w, const double beta )
{
    return kww_hig( w, beta, 0, 0 );
}

Xdouble kwws_hig( const double w, const double beta )
{
    return kww_hig( w, beta, 1, 0 );
}

Xdouble kwwp_hig( const double w, const double beta )
{
    double res = kww_hig( w, beta, 0, 1 );
    if ( res>=PI_2 ) {
        fprintf( stderr, "kwwp: invalid result %g <= 0\n", res );
        exit( ENOSYS );
    }
    return res<0 ? res : PI_2-res;
}


/*****************************************************************************/
/*  Low-level implementation: integration for intermediate frequencies       */
/*****************************************************************************/

#define max_iter_int 12
#define num_range 6
Xdouble kww_mid( const double w, const double beta,
                const int kind, const int mu )
// kind: 0 cos, 1 sin transform (precomputing arrays[2] depend on this)
{
    int iter;
    int kaux;
    int isig;
    int N;
    int j;               // range
    int diffmode;        // subtract Gaussian ?
    Xdouble S=0;     // trapezoid sum
    Xdouble S_last;  // - in last iteration
    Xdouble s;       // term contributing to S
    Xdouble T;       // sum of abs(s)
    // precomputed coefficients
    static int firstCall=1;
    static int iterDone[2][num_range]; // Nm,Np,ak,bk are precomputed up to this
    static int NN[num_range][max_iter_int];
    static Xdouble *ak[2][num_range][max_iter_int];
    static Xdouble *bk[2][num_range][max_iter_int];
    // auxiliary for computing ak and bk
    Xdouble u;
    Xdouble e;
    Xdouble tk;
    Xdouble chi;
    Xdouble dchi;
    Xdouble h;
    Xdouble k;
    Xdouble f;
    Xdouble ahk;
    Xdouble chk;
    Xdouble dhk;
    double p;
    double q;
    const double Smin=2e-20; // to assess worst truncation error

    // dynamic initialization upon first call
    if ( firstCall ) {
        for ( j=0; j<num_range; ++ j ) {
            iterDone[0][j] = -1;
            iterDone[1][j] = -1;
        }
        firstCall = 0;
    }

    // check input
    if ( !( kind==0 || kind==1 ) ) {
        fprintf( stderr, "invalid call to kww_mid: invalid kind\n" );
        exit( EDOM );
    } else if ( beta<0.1 || beta>2.0 ) {
        fprintf( stderr, "invalid call to kww_mid: beta out of range\n" );
        exit( EDOM );
    } else if ( w<=0 ) {
        fprintf( stderr, "invalid call to kww_mid: w out of range\n" );
        exit( EDOM );
    }

    // cosine transform needs special care for beta->2
    if ( kind==0 ) {
        if ( beta==2 )
            return sqrt(PI)/2*exp(-SQR(w)/4);
        diffmode = beta>1.75;
    } else {
        diffmode = 0;
    }

    // determine range, set p,q
    if        ( beta<0.15 ) {
        j=0; p=1.8; q=0.2;
    } else if ( beta<0.25 ) {
        j=1; p=1.6; q=0.4;
    } else if ( beta<1 ) {
        j=2; p=1.4; q=0.6;
    } else if ( beta<1.75 ) {
        j=3; p=1.0; q=0.2;
    } else if ( beta<1.95 ) {
        j=4; p=.75; q=0.2;
    } else {
        j=5; p=.15; q=0.4;
    }

    // iterative integration
    kww_algorithm = 2;
    kww_num_of_terms = 0;
    if( kww_debug & 4 )
        // do not iterate, inspect just one sum
        N = 100;
    else
        N = 40;

    for ( iter=0; iter<max_iter_int; ++iter ) {
        // static initialisation of NN, ak, bk for given 'iter'
        if ( iter>iterDone[kind][j] ) {
            if ( N>1e6 )
                return -3; // integral limits overflow
            NN[j][iter] = N;
            if ( !( ak[kind][j][iter]=malloc((sizeof(Xdouble))*(2*N+1)) ) ||
                 !( bk[kind][j][iter]=malloc((sizeof(Xdouble))*(2*N+1)) )) {
                fprintf( stderr, "kww: Workspace allocation failed\n" );
                exit( ENOMEM );
            }
            iterDone[kind][j] = iter;
            h = logX( logX( 42*N/kww_delta/Smin ) / p ) / N; // 42=(pi+1)*10
            isig=1-2*(NN[j][iter]&1);
            if( kww_debug & 8 ) {
                printf( "init iter %i kind %i j %i siz %i\n",
                        iter, kind, j, 2*N+1 );
            }
            for ( kaux=-NN[j][iter]; kaux<=NN[j][iter]; ++kaux ) {
                k = kaux;
                if( !kind )
                    k -= 0.5;
                u = k*h;
                chi  = 2*p*sinhX(u) + 2*q*u;
                dchi = 2*p*coshX(u) + 2*q;
                if ( u==0 ) {
                    if ( k!=0 )
                        return -4; // integration variable underflow
                    // special treatment to bridge singularity at u=0
                    ahk = PI/h/dchi;
                    dhk = 0.5;
                    chk = sin( ahk );
                } else {
                    if ( -chi>DBL_MAX_EXP/2 )
                        return -5; // integral transformation overflow
                    e = expX( -chi );
                    ahk = PI/h * u/(1-e);
                    dhk = 1/(1-e) - u*e*dchi/SQR(1-e);
                    chk = e>1 ?
                        ( kind ? sinX( PI*k/(1-e) ) : cosX( PI*k/(1-e) ) ) :
                        isig * sinX( PI*k*e/(1-e) );
                }
                ak[kind][j][iter][kaux+NN[j][iter]] = ahk;
                bk[kind][j][iter][kaux+NN[j][iter]] = dhk * chk;
                isig = -isig;
            }
        }
        // integrate according to trapezoidal rule
        S_last = S;
        S = 0;
        T = 0;
        for ( kaux=-NN[j][iter]; kaux<=NN[j][iter]; ++kaux ) {
            tk = ak[kind][j][iter][kaux+NN[j][iter]] / w;
            f = expX(-powX(tk,(Xdouble)beta));
            if ( diffmode )
                f -= expX(-SQR(tk));
            if ( mu )
                f /= tk;
            s = bk[kind][j][iter][kaux+NN[j][iter]] * f;
            S += s;
            T += fabsX(s);
            if( kww_debug & 2 )
                printf( "%2i %6i %12.4Lg %12.4Lg"
                        " %12.4Lg %12.4Lg %12.4Lg %12.4Lg\n",
                        iter, kaux, ak[kind][j][iter][kaux+NN[j][iter]],
                        bk[kind][j][iter][kaux+NN[j][iter]], f, s, S, T );
        }
        if( kww_debug & 1 )
            printf( "%23.17Le  %23.17Le\n", S, T );
        kww_num_of_terms += 2*NN[j][iter]+1;
        if ( diffmode )
            S += w/sqrt(PI)/2*exp(-SQR(w)/4);
        // termination criteria
        if      ( kww_debug & 4 )
            return -1; // we want to inspect just one sum
        else if ( S < 0 && !diffmode )
            return -6; // cancelling terms lead to negative S
        else if ( kww_eps*T > kww_delta*fabsX(S) )
            return -2; // cancellation
        else if ( iter && fabsX(S-S_last) + kww_eps*T < kww_delta*fabsX(S) )
            return S * PI / w; // success (for factor pi/w see my eq. 48)
        N *= 2; // retry with more points
    }
    return -9; // not converged
}

Xdouble kwwc_mid( const double w, const double beta )
{
    return kww_mid( w, beta, 0, 0 );
}

Xdouble kwws_mid( const double w, const double beta )
{
    return kww_mid( w, beta, 1, 0 );
}

Xdouble kwwp_mid( const double w, const double beta )
{
    return kww_mid( w, beta, 1, 1 );
}
