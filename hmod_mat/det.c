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
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "hmod_mat.h"
#include "perm.h"


hlimb_t
_hmod_mat_det(hmod_mat_t A)
{
    hlimb_t det;
    long * P;

    long m = A->r;
    long rank;
    long i;

    P = flint_malloc(sizeof(long) * m);
    rank = hmod_mat_lu(P, A, 1);

    det = 0UL;

    if (rank == m)
    {
        det = 1UL;
        for (i = 0; i < m; i++)
            det = n_mulmod2_preinv(det, hmod_mat_entry(A, i, i),
                A->mod.n, A->mod.ninv);
    }

    if (_perm_parity(P, m) == 1)
        det = nmod_neg(det, A->mod);

    flint_free(P);
    return det;
}

hlimb_t
hmod_mat_det(const hmod_mat_t A)
{
    hmod_mat_t tmp;
    hlimb_t det;
    long dim = A->r;

    if (dim != A->c)
    {
        printf("hmod_mat_det: nonsquare matrix");
        abort();
    }

    if (dim == 0) return 1UL;
    if (dim == 1) return hmod_mat_entry(A, 0, 0);

    hmod_mat_init_set(tmp, A);
    det = _hmod_mat_det(tmp);
    hmod_mat_clear(tmp);

    return det;
}
