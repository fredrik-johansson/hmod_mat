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
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "hmod_mat.h"


static __inline__ int
hmod_mat_pivot(hmod_mat_t A, long * P, long start_row, long col)
{
    long j, t;
    hlimb_t * u;

    if (hmod_mat_entry(A, start_row, col) != 0)
        return 1;

    for (j = start_row + 1; j < A->r; j++)
    {
        if (hmod_mat_entry(A, j, col) != 0)
        {
            u = A->rows[j];
            A->rows[j] = A->rows[start_row];
            A->rows[start_row] = u;

            t = P[j];
            P[j] = P[start_row];
            P[start_row] = t;

            return -1;
        }
    }
    return 0;
}


long 
hmod_mat_lu_classical(long * P, hmod_mat_t A, int rank_check)
{
    hlimb_t d, e, **a;
    nmod_t mod;
    long i, m, n, rank, length, row, col;

    m = A->r;
    n = A->c;
    a = A->rows;
    mod = A->mod;

    rank = row = col = 0;

    for (i = 0; i < m; i++)
        P[i] = i;

    while (row < m && col < n)
    {
        if (hmod_mat_pivot(A, P, row, col) == 0)
        {
            if (rank_check)
                return 0;
            col++;
            continue;
        }

        rank++;

        d = a[row][col];
        d = n_invmod(d, mod.n);
        length = n - col - 1;

        for (i = row + 1; i < m; i++)
        {
            e = n_mulmod2_preinv(a[i][col], d, mod.n, mod.ninv);
            if (length != 0)
                _hmod_vec_scalar_addmul_hmod(a[i] + col + 1,
                    a[row] + col + 1, length, nmod_neg(e, mod), mod);

            a[i][col] = 0;
            a[i][rank - 1] = e;
        }
        row++;
        col++;
    }

    return rank;
}
