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

#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "hmod_mat.h"

long
hmod_mat_nullspace(hmod_mat_t X, const hmod_mat_t A)
{
    long i, j, k, m, n, rank, nullity;
    long * p;
    long * pivots;
    long * nonpivots;
    hmod_mat_t tmp;

    m = A->r;
    n = A->c;

    p = flint_malloc(sizeof(long) * FLINT_MAX(m, n));

    hmod_mat_init_set(tmp, A);
    rank = hmod_mat_rref(tmp);
    nullity = n - rank;

    hmod_mat_zero(X);

    if (rank == 0)
    {
        for (i = 0; i < nullity; i++)
            hmod_mat_entry(X, i, i) = 1UL;
    }
    else if (nullity)
    {
        pivots = p;            /* length = rank */
        nonpivots = p + rank;  /* length = nullity */

        for (i = j = k = 0; i < rank; i++)
        {
            while (hmod_mat_entry(tmp, i, j) == 0UL)
            {
                nonpivots[k] = j;
                k++;
                j++;
            }
            pivots[i] = j;
            j++;
        }
        while (k < nullity)
        {
            nonpivots[k] = j;
            k++;
            j++;
        }

        for (i = 0; i < nullity; i++)
        {
            for (j = 0; j < rank; j++)
            {
                hlimb_t c = hmod_mat_entry(tmp, j, nonpivots[i]);
                hmod_mat_entry(X, pivots[j], i) = nmod_neg(c, A->mod);
            }

            hmod_mat_entry(X, nonpivots[i], i) = 1UL;
        }
    }

    flint_free(p);
    hmod_mat_clear(tmp);

    return nullity;
}
