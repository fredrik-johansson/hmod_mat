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


int hmod_mat_inv(hmod_mat_t B, const hmod_mat_t A)
{
    hmod_mat_t I;
    long i, dim;
    int result;

    dim = A->r;

    switch (dim)
    {
        case 0:
            result = 1;
            break;

        case 1:
            if (hmod_mat_entry(A, 0, 0) == 0UL)
            {
                result = 0;
            }
            else
            {
                hmod_mat_entry(B, 0, 0) = 
                    n_invmod(hmod_mat_entry(A, 0, 0), B->mod.n);
                result = 1;
            }
            break;

        default:
            hmod_mat_init(I, dim, dim, B->mod.n);
            for (i = 0; i < dim; i++)
                hmod_mat_entry(I, i, i) = 1UL;
            result = hmod_mat_solve(B, A, I);
            hmod_mat_clear(I);
    }

    return result;
}
