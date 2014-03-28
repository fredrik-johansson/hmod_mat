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

    Copyright (C) 2010,2011 Fredrik Johansson

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

    printf("solve_triu_classical....");
    fflush(stdout);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        hmod_mat_t A, X, B, Y;
        mp_limb_t m;
        long rows, cols;
        int unit;

        m = hmod_randmod(state);
        rows = n_randint(state, 100);
        cols = n_randint(state, 100);
        unit = n_randint(state, 2);

        hmod_mat_init(A, rows, rows, m);
        hmod_mat_init(B, rows, cols, m);
        hmod_mat_init(X, rows, cols, m);
        hmod_mat_init(Y, rows, cols, m);

        hmod_mat_randtriu(A, state, unit);
        hmod_mat_randtest(X, state);
        hmod_mat_mul(B, A, X);

        /* Check Y = A^(-1) * (A * X) = X */
        hmod_mat_solve_triu_classical(Y, A, B, unit);
        if (!hmod_mat_equal(Y, X))
        {
            printf("FAIL!\n");
            printf("A:\n");
            hmod_mat_print_pretty(A);
            printf("X:\n");
            hmod_mat_print_pretty(X);
            printf("B:\n");
            hmod_mat_print_pretty(B);
            printf("Y:\n");
            hmod_mat_print_pretty(Y);
            abort();
        }

        /* Check aliasing */
        hmod_mat_solve_triu_classical(B, A, B, unit);
        if (!hmod_mat_equal(B, X))
        {
            printf("FAIL!\n");
            printf("aliasing test failed");
            printf("A:\n");
            hmod_mat_print_pretty(A);
            printf("B:\n");
            hmod_mat_print_pretty(B);
            abort();
        }

        hmod_mat_clear(A);
        hmod_mat_clear(B);
        hmod_mat_clear(X);
        hmod_mat_clear(Y);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
