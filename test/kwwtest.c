/* runkww.c
 *
 * Copyright (C) 2009-2019 Joachim Wuttke
 *
 * Licence: GNU General Public License, version 3 or later
 *
 * Author:
 *   Joachim Wuttke, Forschungszentrum Jülich, Germany <j.wuttke@fz-juelich.de>
 *
 * Purpose:
 *   Test whether numeric results agree with those established under Linux/gcc.
 */

#include "kww.h"

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <float.h>

/******************************************************************************/
/*  Auxiliary routines                                                        */
/******************************************************************************/

// compute relative error |b-a|/|a|, handling case of NaN and Inf,
static double relerr(double a, double b) {
    assert(a!=0);
    assert(b!=0);
    return fabs((b-a) / a);
}

void test_one(int* fail, double limit, double a, double b)
{
    static int num = 0;
    ++num;
    double re = relerr(a,b);
    if (re==0)
        return;
    int ok = re<=limit;
    printf(ok ? "   " : "ERR");
    printf(" test case %i: found=%g, expected=%g, relerr=%g\n", num, a, b, re);
    if (!ok)
        ++(*fail);
}


/******************************************************************************/
/*  Main: test sequence                                                       */
/******************************************************************************/

int main(void) {
    // Test specific function values.
    int fail = 0;

    test_one(&fail, 1e-14, kwwc(1e-6, 1.), 0.9999999999990000221);
    test_one(&fail, 1e-14, kwwc(1e-3, 1.), 0.9999990000009999491);
    test_one(&fail, 1e-14, kwwc(1   , 1.), .5);
    test_one(&fail, 1e-14, kwwc(1e3 , 1.), 9.999990000010000613e-07);
    test_one(&fail, 1e-14, kwwc(1e6 , 1.), 9.999999999989999315e-13);

    test_one(&fail, 1e-14, kwwc(1e-8, .623), 1.435159133351523009);
    test_one(&fail, 1e-14, kwwc(1e-6, .623), 1.435159133336882054);
    test_one(&fail, 1e-14, kwwc(1e-4, .623), 1.435158986927334901);
    test_one(&fail, 1e-14, kwwc(1e-2, .623), 1.433698427082435778);     // mode=1
    test_one(&fail, 1e-14, kwwc(1e-1, .623), 1.314878935071708499);     // mode=1
    test_one(&fail, 1e-14, kwwc(1   , .623), 0.3308323687099594124);    // mode=3
    test_one(&fail, 1e-14, kwwc(1e2 , .623), 0.0004053330090102800066);
    test_one(&fail, 1e-14, kwwc(1e4 , .623), 2.390040041093056597e-07);

    test_one(&fail, 1e-14, kwwc(2e-5, .314), 7.602900889248060956); // mode=1
    test_one(&fail, 1e-14, kwwc(2e-4, .314), 7.594626504743104078); // mode=2
    test_one(&fail, 1e-14, kwwc(2e-3, .314), 7.148958376075823296); // mode=2
    test_one(&fail, 1e-14, kwwc(2e-2, .314), 3.922292835319648674); // mode=3
    test_one(&fail, 1e-14, kwwc(2e-1, .314), 0.8172678260275950679);
    test_one(&fail, 1e-14, kwwc(2e0 , .314), 0.0837197392901634224);

    test_one(&fail, 1e-14, kwws(2e-5, .314), 0.01452905498350009518); // mode=1
    test_one(&fail, 1e-14, kwws(2e-3, .314), 1.11472540368966655); // mode=2
    test_one(&fail, 1e-14, kwws(2e-1, .314), 1.202467631193444353); // mode=3

    test_one(&fail, 1e-14, kwwp(3e-3, .459), 0.007116055704011668009); // mode=1
    test_one(&fail, 1e-14, kwwp(5e-3, .459), 0.01185130685163975767); // mode=2
    test_one(&fail, 5e-12, kwwp(2e-2, .459), 0.04668285680895551543); // mode=1

    printf("\n");
    if (fail) {
        printf("IN TOTAL, FAILURE IN %i TESTS\n", fail);
        return 1;
    } else {
        printf("OVERALL SUCCESS\n");
        return 0;
    }
}
