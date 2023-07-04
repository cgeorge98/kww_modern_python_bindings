/* lowlevel.h:
 *   Computation either by series expansion (low/hig) or by integration (mid).
 *   Exported for use in demokww only.
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

#ifndef __KWW_LOWLEVEL_H__
#define __KWW_LOWLEVEL_H__

#include "extended_double.h"

/* low-w expansion */
KWW_EXPORT Xdouble kwwc_low( const double w, const double beta );
KWW_EXPORT Xdouble kwws_low( const double w, const double beta );
KWW_EXPORT Xdouble kwwp_low( const double w, const double beta );

/* high-w expansion */
KWW_EXPORT Xdouble kwwc_hig( const double w, const double beta );
KWW_EXPORT Xdouble kwws_hig( const double w, const double beta );
KWW_EXPORT Xdouble kwwp_hig( const double w, const double beta );

/* mid-w integration */
KWW_EXPORT Xdouble kwwc_mid( const double w, const double beta );
KWW_EXPORT Xdouble kwws_mid( const double w, const double beta );
KWW_EXPORT Xdouble kwwp_mid( const double w, const double beta );

__END_DECLS
#endif /* __KWW_LOWLEVEL_H__ */
