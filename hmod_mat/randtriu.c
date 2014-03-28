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

void
hmod_mat_randtriu(hmod_mat_t mat, flint_rand_t state, int unit)
{
    long i, j;

    for (i = 0; i < mat->r; i++)
    {
        for (j = 0; j < mat->c; j++)
        {
            if (j > i)
            {
                hmod_mat_entry(mat, i, j) = n_randlimb(state) % (mat->mod.n);
            }
            else if (i == j)
            {
                hmod_mat_entry(mat, i, j) = n_randlimb(state) % (mat->mod.n);
                if (unit || hmod_mat_entry(mat, i, j) == 0UL)
                    hmod_mat_entry(mat, i, j) = 1UL;
            }
            else
            {
                hmod_mat_entry(mat, i, j) = 0UL;
            }
        }
    }
}
