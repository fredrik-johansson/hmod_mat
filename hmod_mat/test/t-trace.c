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
#include <mpir.h>
#include "flint.h"
#include "hmod_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    long i;
    flint_rand_t state;

    printf("trace....");
    fflush(stdout);

    flint_randinit(state);

    /* Test trace(AB) = trace(BA) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        hmod_mat_t A, B, AB, BA;
        mp_limb_t mod, trab, trba;
        long m, n;

        mod = hmod_randmod(state);
        m = n_randint(state, 10);
        n = n_randint(state, 10);

        hmod_mat_init(A, m, n, mod);
        hmod_mat_init(B, n, m, mod);
        hmod_mat_init(AB, m, m, mod);
        hmod_mat_init(BA, n, n, mod);

        hmod_mat_randtest(A, state);
        hmod_mat_randtest(B, state);

        hmod_mat_mul(AB, A, B);
        hmod_mat_mul(BA, B, A);

        trab = hmod_mat_trace(AB);
        trba = hmod_mat_trace(BA);

        if (trab != trba)
        {
            printf("FAIL:\n");
            hmod_mat_print_pretty(A), printf("\n");
            hmod_mat_print_pretty(B), printf("\n");
            hmod_mat_print_pretty(AB), printf("\n");
            hmod_mat_print_pretty(BA), printf("\n");
            printf("tr(AB): %lu\n", trab);
            printf("tr(BA): %lu\n", trba);
            abort();
        }

        hmod_mat_clear(A);
        hmod_mat_clear(B);
        hmod_mat_clear(AB);
        hmod_mat_clear(BA);
    }

    flint_randclear(state);
    printf("PASS\n");
    return 0;
}
