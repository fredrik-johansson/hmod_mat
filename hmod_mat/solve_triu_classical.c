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

#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "hmod_mat.h"
#include "nmod_vec.h"

void
hmod_mat_solve_triu_classical(hmod_mat_t X, const hmod_mat_t U,
                                                const hmod_mat_t B, int unit)
{
    int nlimbs;
    long i, j, n, m;
    nmod_t mod;
    hlimb_t * inv, * tmp;

    n = U->r;
    m = B->c;
    mod = U->mod;

    if (!unit)
    {
        inv = _nmod_vec_init(n);
        for (i = 0; i < n; i++)
            inv[i] = n_invmod(hmod_mat_entry(U, i, i), mod.n);
    }
    else
        inv = NULL;

    nlimbs = _hmod_vec_dot_bound_limbs(n, mod);
    tmp = _nmod_vec_init(n);

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
            tmp[j] = hmod_mat_entry(X, j, i);

        for (j = n - 1; j >= 0; j--)
        {
            hlimb_t s;
            s = _hmod_vec_dot(U->rows[j] + j + 1,
                              tmp + j + 1, n - j - 1, mod, nlimbs);
            s = nmod_sub(hmod_mat_entry(B, j, i), s, mod);
            if (!unit)
                s = n_mulmod2_preinv(s, inv[j], mod.n, mod.ninv);
            tmp[j] = s;
        }

        for (j = 0; j < n; j++)
            hmod_mat_entry(X, j, i) = tmp[j];
    }

    _nmod_vec_clear(tmp);
    if (!unit)
        _nmod_vec_clear(inv);
}
