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
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "hmod_mat.h"


int
hmod_mat_solve(hmod_mat_t X, const hmod_mat_t A, const hmod_mat_t B)
{
    long i, rank, *perm;
    hmod_mat_t LU;
    int result;

    if (A->r == 0 || B->c == 0)
        return 1;

    hmod_mat_init_set(LU, A);
    perm = flint_malloc(sizeof(long) * A->r);
    for (i = 0; i < A->r; i++)
        perm[i] = i;

    rank = hmod_mat_lu(perm, LU, 1);

    if (rank == A->r)
    {
        hmod_mat_t PB;
        hmod_mat_window_init(PB, B, 0, 0, B->r, B->c);
        for (i = 0; i < A->r; i++)
            PB->rows[i] = B->rows[perm[i]];

        hmod_mat_solve_tril(X, LU, PB, 1);
        hmod_mat_solve_triu(X, LU, X, 0);

        hmod_mat_window_clear(PB);
        result = 1;
    }
    else
    {
        result = 0;
    }

    hmod_mat_clear(LU);
    flint_free(perm);

    return result;
}
