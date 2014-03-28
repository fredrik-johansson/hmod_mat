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

    Copyright (C) 2010,2012 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "hmod_mat.h"
#include "nmod_vec.h"

/*
with op = 0, computes D = A*B
with op = 1, computes D = C + A*B
with op = -1, computes D = C - A*B
*/

static __inline__ void
_hmod_mat_addmul_basic(hlimb_t ** D, hlimb_t ** const C, hlimb_t ** const A,
    hlimb_t ** const B, long m, long k, long n, int op, nmod_t mod, int nlimbs)
{
    long i, j;
    hlimb_t c;

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            c = _hmod_vec_dot_ptr(A[i], B, j, k, mod, nlimbs);

            if (op == 1)
                c = nmod_add(C[i][j], c, mod);
            else if (op == -1)
                c = nmod_sub(C[i][j], c, mod);

            D[i][j] = c;
        }
    }
}

static __inline__ void
_hmod_mat_addmul_transpose(hlimb_t ** D, hlimb_t ** const C, hlimb_t ** const A,
    hlimb_t ** const B, long m, long k, long n, int op, nmod_t mod, int nlimbs)
{
    hlimb_t * tmp;
    hlimb_t c;
    long i, j;

    tmp = flint_malloc(sizeof(hlimb_t) * k * n);

    for (i = 0; i < k; i++)
        for (j = 0; j < n; j++)
            tmp[j*k + i] = B[i][j];

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            c = _hmod_vec_dot(A[i], tmp + j*k, k, mod, nlimbs);

            if (op == 1)
                c = nmod_add(C[i][j], c, mod);
            else if (op == -1)
                c = nmod_sub(C[i][j], c, mod);

            D[i][j] = c;
        }
    }

    flint_free(tmp);
}

void
_hmod_mat_mul_classical(hmod_mat_t D, const hmod_mat_t C,
                                const hmod_mat_t A, const hmod_mat_t B, int op)
{
    long m, k, n;
    int nlimbs;
    nmod_t mod;

    mod = A->mod;
    m = A->r;
    k = A->c;
    n = B->c;

    if (k == 0)
    {
        if (op == 0)
            hmod_mat_zero(D);
        else
            hmod_mat_set(D, C);
        return;
    }

    nlimbs = _hmod_vec_dot_bound_limbs(k, mod);

    if (m < HMOD_MAT_MUL_TRANSPOSE_CUTOFF
        || n < HMOD_MAT_MUL_TRANSPOSE_CUTOFF
        || k < HMOD_MAT_MUL_TRANSPOSE_CUTOFF)
    {
        _hmod_mat_addmul_basic(D->rows, (op == 0) ? NULL : C->rows,
            A->rows, B->rows, m, k, n, op, D->mod, nlimbs);
    }
    else
    {
        _hmod_mat_addmul_transpose(D->rows, (op == 0) ? NULL : C->rows,
            A->rows, B->rows, m, k, n, op, D->mod, nlimbs);
    }
}


void
hmod_mat_mul_classical(hmod_mat_t C, const hmod_mat_t A, const hmod_mat_t B)
{
    _hmod_mat_mul_classical(C, NULL, A, B, 0);
}
