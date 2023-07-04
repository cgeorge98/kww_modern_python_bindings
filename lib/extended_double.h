/* extended_double.h
 *   Macros that expand to either "long double" or to "__float128".
 *   To facilitate cross-platform computations with at least extended precision.
 *
 * Copyright:
 *   (C) 2022 Joachim Wuttke
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
 */

#ifndef EXTENDED_DOUBLE_H
#define EXTENDED_DOUBLE_H

#ifdef USE_FLOAT128

    #include <quadmath.h>
    #define Xdouble __float128

    #define cosX cosq
    #define coshX coshq
    #define expX expq
    #define fabsX fabsq
    #define lgammaX lgammaq
    #define logX logq
    #define powX powq
    #define sinX sinq
    #define sinhX sinhq

#else // use long double

    #include <assert.h>
    static_assert(sizeof(long double)>=10, "long double shorter than 80 bits");

    #define Xdouble long double

    #define cosX cosl
    #define coshX coshl
    #define expX expl
    #define fabsX fabsl
    #define lgammaX lgammal
    #define logX logl
    #define powX powl
    #define sinX sinl
    #define sinhX sinhl

#endif

#endif // EXTENDED_DOUBLE_H
