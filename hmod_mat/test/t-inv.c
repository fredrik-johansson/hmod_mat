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
    hmod_mat_t A, B, C, I;
    long i, j, m, r;
    mp_limb_t mod;
    int result;
    flint_rand_t state;
    flint_randinit(state);

    printf("inv....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 20);
        mod = hmod_randmod(state);

        hmod_mat_init(A, m, m, mod);
        hmod_mat_init(B, m, m, mod);
        hmod_mat_init(C, m, m, mod);
        hmod_mat_init(I, m, m, mod);

        for (j = 0; j < m; j++)
            I->rows[j][j] = 1UL;

        /* Verify that A * A^-1 = I for random matrices */

        hmod_mat_randrank(A, state, m);
        /* Dense or sparse? */
        if (n_randint(state, 2))
            hmod_mat_randops(A, 1+n_randint(state, 1+m*m), state);

        result = hmod_mat_inv(B, A);
        hmod_mat_mul(C, A, B);

        if (!hmod_mat_equal(C, I) || !result)
        {
            printf("FAIL:\n");
            printf("A * A^-1 != I!\n");
            printf("A:\n");
            hmod_mat_print_pretty(A);
            printf("A^-1:\n");
            hmod_mat_print_pretty(B);
            printf("A * A^-1:\n");
            hmod_mat_print_pretty(C);
            printf("\n");
            abort();
        }

        /* Test aliasing */
        hmod_mat_set(C, A);
        hmod_mat_inv(A, A);
        hmod_mat_mul(B, A, C);

        if (!hmod_mat_equal(B, I))
        {
            printf("FAIL:\n");
            printf("aliasing failed!\n");
            hmod_mat_print_pretty(C);
            abort();
        }

        hmod_mat_clear(A);
        hmod_mat_clear(B);
        hmod_mat_clear(C);
        hmod_mat_clear(I);
    }

    /* Test singular systems */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        m = 1 + n_randint(state, 20);
        mod = hmod_randmod(state);
        r = n_randint(state, m);

        hmod_mat_init(A, m, m, mod);
        hmod_mat_init(B, m, m, mod);

        hmod_mat_randrank(A, state, r);

        /* Dense */
        if (n_randint(state, 2))
            hmod_mat_randops(A, 1+n_randint(state, 1+m*m), state);

        result = hmod_mat_inv(B, A);

        if (result)
        {
            printf("FAIL:\n");
            printf("singular matrix reported as invertible\n");
            abort();
        }

        /* Aliasing */
        result = hmod_mat_inv(A, A);
        if (result)
        {
            printf("FAIL:\n");
            printf("singular matrix reported as invertiblen");
            abort();
        }

        hmod_mat_clear(A);
        hmod_mat_clear(B);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
