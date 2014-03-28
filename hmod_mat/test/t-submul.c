/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <mpir.h>
#include "flint.h"
#include "hmod_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    long i;
    flint_rand_t state;
    flint_randinit(state);

    printf("submul....");
    fflush(stdout);

    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        hmod_mat_t A, B, C, D, T, E;
        mp_limb_t mod = hmod_randmod(state);

        long m, k, n;

        m = n_randint(state, 100);
        k = n_randint(state, 100);
        n = n_randint(state, 100);

        /* Force Strassen test */
        if (i < 5)
        {
            m += 300;
            k += 300;
            n += 300;
        }

        hmod_mat_init(A, m, k, mod);
        hmod_mat_init(B, k, n, mod);
        hmod_mat_init(C, m, n, mod);
        hmod_mat_init(D, m, n, mod);
        hmod_mat_init(T, m, n, mod);
        hmod_mat_init(E, m, n, mod);

        hmod_mat_randtest(A, state);
        hmod_mat_randtest(B, state);
        hmod_mat_randtest(C, state);

        hmod_mat_submul(D, C, A, B);

        hmod_mat_mul(T, A, B);
        hmod_mat_sub(E, C, T);

        if (!hmod_mat_equal(D, E))
        {
            printf("FAIL: results not equal\n");
            hmod_mat_print_pretty(A);
            hmod_mat_print_pretty(B);
            hmod_mat_print_pretty(C);
            hmod_mat_print_pretty(D);
            hmod_mat_print_pretty(E);
            abort();
        }

        /* Check aliasing */
        hmod_mat_submul(C, C, A, B);

        if (!hmod_mat_equal(C, E))
        {
            printf("FAIL: results not equal (aliasing)\n");
            hmod_mat_print_pretty(A);
            hmod_mat_print_pretty(B);
            hmod_mat_print_pretty(C);
            hmod_mat_print_pretty(D);
            hmod_mat_print_pretty(E);
            abort();
        }

        hmod_mat_clear(A);
        hmod_mat_clear(B);
        hmod_mat_clear(C);
        hmod_mat_clear(D);
        hmod_mat_clear(E);
        hmod_mat_clear(T);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
