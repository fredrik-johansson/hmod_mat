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

void
hmod_mat_mul_check(hmod_mat_t C, const hmod_mat_t A, const hmod_mat_t B)
{
    long i, j, k;

    mp_limb_t s0, s1, s2;
    mp_limb_t t0, t1;

    for (i = 0; i < A->r; i++)
    {
        for (j = 0; j < B->c; j++)
        {
            s0 = s1 = s2 = 0UL;

            for (k = 0; k < A->c; k++)
            {
                umul_ppmm(t1, t0, (mp_limb_t) A->rows[i][k], (mp_limb_t) B->rows[k][j]);
                add_sssaaaaaa(s2, s1, s0, s2, s1, s0, 0, t1, t0);
            }

            NMOD_RED(s2, s2, C->mod);
            NMOD_RED3(s0, s2, s1, s0, C->mod);
            C->rows[i][j] = s0;
        }
    }
}

int
main(void)
{
    long i;
    flint_rand_t state;
    flint_randinit(state);

    printf("mul....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        hmod_mat_t A, B, C, D;
        mp_limb_t mod;

        long m, k, n;

        m = n_randint(state, 50);
        k = n_randint(state, 50);
        n = n_randint(state, 50);

        /* We want to generate matrices with many entries close to half
           or full limbs with high probability, to stress overflow handling */
        switch (n_randint(state, 3))
        {
            case 0:
                mod = hmod_randmod(state);
                break;
            case 1:
                mod = UINT_MAX/2 + 1 - n_randbits(state, 4);
                break;
            case 2:
            default:
                mod = UINT_MAX - n_randbits(state, 4);
                break;
        }

        hmod_mat_init(A, m, n, mod);
        hmod_mat_init(B, n, k, mod);
        hmod_mat_init(C, m, k, mod);
        hmod_mat_init(D, m, k, mod);

        if (n_randint(state, 2))
            hmod_mat_randtest(A, state);
        else
            hmod_mat_randfull(A, state);

        if (n_randint(state, 2))
            hmod_mat_randtest(B, state);
        else
            hmod_mat_randfull(B, state);

        hmod_mat_randtest(C, state);  /* make sure noise in the output is ok */

        hmod_mat_mul(C, A, B);
        hmod_mat_mul_check(D, A, B);

        if (!hmod_mat_equal(C, D))
        {
            printf("FAIL: results not equal\n");
            hmod_mat_print_pretty(A);
            hmod_mat_print_pretty(B);
            hmod_mat_print_pretty(C);
            hmod_mat_print_pretty(D);
            abort();
        }

        hmod_mat_clear(A);
        hmod_mat_clear(B);
        hmod_mat_clear(C);
        hmod_mat_clear(D);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
