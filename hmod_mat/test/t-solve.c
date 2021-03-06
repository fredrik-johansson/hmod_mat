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
#include <mpir.h>
#include "flint.h"
#include "nmod_vec.h"
#include "hmod_mat.h"
#include "ulong_extras.h"


int
main(void)
{
    hmod_mat_t A, X, B, AX;
    long i, m, n, r;
    mp_limb_t mod;
    int solved;
    flint_rand_t state;
    flint_randinit(state);

    printf("solve....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 20);
        n = n_randint(state, 20);
        mod = hmod_randmod(state);

        hmod_mat_init(A, m, m, mod);
        hmod_mat_init(B, m, n, mod);
        hmod_mat_init(X, m, n, mod);
        hmod_mat_init(AX, m, n, mod);

        hmod_mat_randrank(A, state, m);
        hmod_mat_randtest(B, state);

        /* Dense */
        if (n_randint(state, 2))
            hmod_mat_randops(A, 1+n_randint(state, 1+m*m), state);

        solved = hmod_mat_solve(X, A, B);

        hmod_mat_mul(AX, A, X);

        if (!hmod_mat_equal(AX, B) || !solved)
        {
            printf("FAIL:\n");
            printf("AX != B!\n");
            printf("A:\n");
            hmod_mat_print_pretty(A);
            printf("B:\n");
            hmod_mat_print_pretty(B);
            printf("X:\n");
            hmod_mat_print_pretty(X);
            printf("AX:\n");
            hmod_mat_print_pretty(AX);
            printf("\n");
            abort();
        }

        hmod_mat_clear(A);
        hmod_mat_clear(B);
        hmod_mat_clear(X);
        hmod_mat_clear(AX);
    }

    /* Test singular systems */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        m = 1 + n_randint(state, 20);
        n = 1 + n_randint(state, 20);
        r = n_randint(state, m);
        mod = hmod_randmod(state);

        hmod_mat_init(A, m, m, mod);
        hmod_mat_init(B, m, n, mod);
        hmod_mat_init(X, m, n, mod);
        hmod_mat_init(AX, m, n, mod);

        hmod_mat_randrank(A, state, r);
        hmod_mat_randtest(B, state);

        /* Dense */
        if (n_randint(state, 2))
            hmod_mat_randops(A, 1+n_randint(state, 1+m*m), state);

        solved = hmod_mat_solve(X, A, B);

        if (solved)
        {
            printf("FAIL:\n");
            printf("singular system was 'solved'\n");
            hmod_mat_print_pretty(A);
            hmod_mat_print_pretty(X);
            hmod_mat_print_pretty(B);
            abort();
        }

        hmod_mat_clear(A);
        hmod_mat_clear(B);
        hmod_mat_clear(X);
        hmod_mat_clear(AX);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
