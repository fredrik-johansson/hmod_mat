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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <mpir.h>
#include "flint.h"
#include "hmod_mat.h"
#include "nmod_vec.h"
#include "ulong_extras.h"

int
main(void)
{
    long m, n, mod, rep;
    flint_rand_t state;
    flint_randinit(state);

    printf("scalar_mul....");
    fflush(stdout);

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        hmod_mat_t A, B, C, D;
        mp_limb_t c;

        m = n_randint(state, 20);
        n = n_randint(state, 20);
        mod = hmod_randmod(state);

        c = n_randint(state, mod);

        hmod_mat_init(A, m, n, mod);
        hmod_mat_init(B, m, n, mod);
        hmod_mat_init(C, m, n, mod);
        hmod_mat_init(D, m, n, mod);

        hmod_mat_randtest(A, state);
        hmod_mat_randtest(B, state);

        hmod_mat_scalar_mul(C, A, c);
        hmod_mat_scalar_mul(D, A, nmod_sub(c, 1UL, A->mod));

        /* c*A - (c-1)*A == A */
        hmod_mat_sub(D, C, D);

        if (!hmod_mat_equal(A, D))
        {
            printf("FAIL\n");
            abort();
        }

        /* Aliasing */
        hmod_mat_scalar_mul(C, A, c);
        hmod_mat_scalar_mul(A, A, c);

        if (!hmod_mat_equal(A, C))
        {
            printf("FAIL\n");
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
