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
#include "nmod_vec.h"

hlimb_t
hmod_mat_trace(const hmod_mat_t mat)
{
    hlimb_t t;
    long i, n = hmod_mat_nrows(mat);

    if (n == 0)
        return 0;

    t = hmod_mat_entry(mat, 0, 0);

    for (i = 1; i < n; i++)
        t = nmod_add(t, hmod_mat_entry(mat, i, i), mat->mod);

    return t;
}
