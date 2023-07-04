/* kww.i:
 *   SWIG interface file for the kww library .
 * 
 * Copyright:
 *   (C) 2014 Antti Soininen
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
 * Contact:
 *   Joachim Wuttke
 *   Forschungszentrum Jülich, Germany
 *   j.wuttke@fz-juelich.de
 *
 * Website:
 *   http://apps.jcns.fz-juelich.de/doku/sc/kww
 *
 */

%module kww
%{
extern double kwwc( const double w, const double beta );
extern double kwws( const double w, const double beta );
extern double kwwp( const double w, const double beta );
%}

extern double kwwc( const double w, const double beta );
extern double kwws( const double w, const double beta );
extern double kwwp( const double w, const double beta );

