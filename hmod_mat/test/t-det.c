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
#include "fmpz.h"
#include "fmpz_mat.h"

int
main(void)
{
    long m, rep;
    flint_rand_t state;
    flint_randinit(state);

    printf("det....");
    fflush(stdout);

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        hmod_mat_t A;
        fmpz_mat_t B;
        mp_limb_t Adet, mod;
        fmpz_t Bdet;
        ulong t;

        m = n_randint(state, 30);
        mod = hmod_randmod(state);

        hmod_mat_init(A, m, m, mod);
        fmpz_mat_init(B, m, m);

        switch (rep % 3)
        {
            case 0:
                hmod_mat_randrank(A, state, m);
                hmod_mat_randops(A, n_randint(state, 2*m + 1), state);
                break;
            case 1:
                t = n_randint(state, m);
                t = FLINT_MIN(t, m);
                hmod_mat_randrank(A, state, m);
                hmod_mat_randops(A, n_randint(state, 2*m + 1), state);
                break;
            default:
                hmod_mat_randtest(A, state);
        }

/*        fmpz_mat_set_hmod_mat_unsigned(B, A); */
        {
            long j, k;

            for (j = 0; j < A->r; j++)
                for (k = 0; k < A->c; k++)
                    fmpz_set_ui(B->rows[j] + k, A->rows[j][k]);
        }

        Adet = hmod_mat_det(A);

        fmpz_init(Bdet);
        fmpz_mat_det_bareiss(Bdet, B);
        fmpz_mod_ui(Bdet, Bdet, mod);

        if (Adet != fmpz_get_ui(Bdet))
        {
            printf("FAIL\n");
            abort();
        }

        hmod_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_clear(Bdet);
    }

    flint_randclear(state);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
