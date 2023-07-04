/* kww.h:
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

#ifndef __KWW_H__
#define __KWW_H__
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif
__BEGIN_DECLS

#if _WIN32
#define KWW_EXPORT __declspec(dllexport)
#define KWW_IMPORT __declspec(dllimport)
#else
#define KWW_EXPORT
#define KWW_IMPORT
#endif

/*****************************************************************************/
/*  High-level calls                                                         */
/*****************************************************************************/

/* \int_0^\infty dt cos(w*t) exp(-t^beta) */
KWW_EXPORT double kwwc( const double w, const double beta );

/* \int_0^\infty dt sin(w*t) exp(-t^beta) */
KWW_EXPORT double kwws( const double w, const double beta );

/* \int_0^w dw' kwwc(w') */
KWW_EXPORT double kwwp( const double w, const double beta );


/*****************************************************************************/
/*  Low-level calls                                                          */
/*****************************************************************************/

/* range limits */
KWW_EXPORT double kwwc_lim_low( const double b );
KWW_EXPORT double kwwc_lim_hig( const double b );
KWW_EXPORT double kwws_lim_low( const double b );
KWW_EXPORT double kwws_lim_hig( const double b );
KWW_EXPORT double kwwp_lim_low( const double b );
KWW_EXPORT double kwwp_lim_hig( const double b );

__END_DECLS
#endif /* __KWW_H__ */
