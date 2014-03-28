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

    Copyright (C) 2010 William Hart
    Copyright (C) 2010,2011 Fredrik Johansson

******************************************************************************/

#ifndef HMOD_MAT_H
#define HMOD_MAT_H

#undef ulong /* interferes with system includes */
#include <stdlib.h>
#define ulong unsigned long

#include <mpir.h>
#include "flint.h"
#include "longlong.h"
#include "ulong_extras.h"
#include "nmod_vec.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef unsigned int hlimb_t;

static __inline__ mp_limb_t hmod_randmod(flint_rand_t state)
{
    mp_limb_t p;
    do {
        p = n_randtest_prime(state, 0);
    } while (p >= (1UL << (FLINT_BITS / 2)));
    return p;
}

static __inline__ void _hmod_vec_neg(hlimb_t * a, const hlimb_t * b, long n, nmod_t mod)
{
    long i;
    for (i = 0; i < n; i++)
        a[i] = n_negmod(b[i], mod.n);
}

static __inline__ void _hmod_vec_add(hlimb_t * a, const hlimb_t * b,  const hlimb_t * c, long n, nmod_t mod)
{
    long i;
    for (i = 0; i < n; i++)
        a[i] = n_addmod(b[i], c[i], mod.n);
}

static __inline__ void _hmod_vec_sub(hlimb_t * a, const hlimb_t * b,  const hlimb_t * c, long n, nmod_t mod)
{
    long i;
    for (i = 0; i < n; i++)
        a[i] = n_submod(b[i], c[i], mod.n);
}

static __inline__ void _hmod_vec_scalar_addmul_hmod(hlimb_t * res, const hlimb_t * vec, 
                            long len, hlimb_t c, nmod_t mod)
{
    long i;
    for (i = 0; i < len; i++)
        NMOD_ADDMUL(res[i], vec[i], c, mod);
}

/* right now we only care about this case */
static __inline__ int _hmod_vec_dot_bound_limbs(long len, nmod_t mod)
{
    return 2;
}

static __inline__ hlimb_t
_hmod_vec_dot(const hlimb_t * vec1, const hlimb_t * vec2, long len, nmod_t mod, int nlimbs)
{
    hlimb_t res;
    long i;
    mp_limb_t s0, s1, t0;
    s0 = s1 = 0UL;

    for (i = 0; i < len; i++)
    {
        t0 = (mp_limb_t) vec1[i] * (mp_limb_t) vec2[i];
        add_ssaaaa(s1, s0, s1, s0, 0, t0);
    }

    NMOD2_RED2(s0, s1, s0, mod);
    res = s0;
    return res;
}

static __inline__ hlimb_t
_hmod_vec_dot_ptr(const hlimb_t * vec1, hlimb_t ** const vec2, long offset, long len, nmod_t mod, int nlimbs)
{
    hlimb_t res;
    long i;
    mp_limb_t s0, s1, t0;
    s0 = s1 = 0UL;

    for (i = 0; i < len; i++)
    {
        t0 = (mp_limb_t) vec1[i] * (mp_limb_t) vec2[i][offset];
        add_ssaaaa(s1, s0, s1, s0, 0, t0);
    }

    NMOD2_RED2(s0, s1, s0, mod);
    res = s0;
    return res;
}


typedef struct
{
    hlimb_t * entries;
    long r;
    long c;
    hlimb_t ** rows;
    nmod_t mod;
}
hmod_mat_struct;

/* hmod_mat_t allows reference-like semantics for hmod_mat_struct */
typedef hmod_mat_struct hmod_mat_t[1];

#define hmod_mat_entry(mat,i,j) ((mat)->rows[(i)][(j)])
#define hmod_mat_nrows(mat) ((mat)->r)
#define hmod_mat_ncols(mat) ((mat)->c)

static __inline__
void
_hmod_mat_set_mod(hmod_mat_t mat, mp_limb_t n)
{
    mat->mod.n = n;
    mat->mod.ninv = n_preinvert_limb(n);
    count_leading_zeros(mat->mod.norm, n);
}

/* Memory management */
void hmod_mat_init(hmod_mat_t mat, long rows, long cols, hlimb_t n);
void hmod_mat_init_set(hmod_mat_t mat, const hmod_mat_t src);
void hmod_mat_clear(hmod_mat_t mat);

void hmod_mat_window_init(hmod_mat_t window, const hmod_mat_t mat, long r1, long c1, long r2, long c2);
void hmod_mat_window_clear(hmod_mat_t window);

/* Random matrix generation */
void hmod_mat_randtest(hmod_mat_t mat, flint_rand_t state);
void hmod_mat_randfull(hmod_mat_t mat, flint_rand_t state);
int hmod_mat_randpermdiag(hmod_mat_t mat, flint_rand_t state, 
                 mp_srcptr diag, long n);
void hmod_mat_randrank(hmod_mat_t, flint_rand_t state, long rank);
void hmod_mat_randops(hmod_mat_t mat, long count, flint_rand_t state);
void hmod_mat_randtril(hmod_mat_t mat, flint_rand_t state, int unit);
void hmod_mat_randtriu(hmod_mat_t mat, flint_rand_t state, int unit);


void hmod_mat_print_pretty(const hmod_mat_t mat);

int hmod_mat_equal(const hmod_mat_t mat1, const hmod_mat_t mat2);

void hmod_mat_zero(hmod_mat_t mat);

int hmod_mat_is_zero(const hmod_mat_t mat);

static __inline__ int
hmod_mat_is_empty(const hmod_mat_t mat)
{
    return (mat->r == 0) || (mat->c == 0);
}

static __inline__ int
hmod_mat_is_square(const hmod_mat_t mat)
{
    return (mat->r == mat->c);
}


void hmod_mat_set(hmod_mat_t B, const hmod_mat_t A);
void hmod_mat_transpose(hmod_mat_t B, const hmod_mat_t A);

/* Addition and subtraction */

void hmod_mat_add(hmod_mat_t C, const hmod_mat_t A, const hmod_mat_t B);
void hmod_mat_sub(hmod_mat_t C, const hmod_mat_t A, const hmod_mat_t B);
void hmod_mat_neg(hmod_mat_t B, const hmod_mat_t A);

/* Matrix-scalar arithmetic */

void hmod_mat_scalar_mul(hmod_mat_t B, const hmod_mat_t A, hlimb_t c);

/* Matrix multiplication */

void hmod_mat_mul(hmod_mat_t C, const hmod_mat_t A, const hmod_mat_t B);
void hmod_mat_mul_classical(hmod_mat_t C, const hmod_mat_t A, const hmod_mat_t B);
void hmod_mat_mul_strassen(hmod_mat_t C, const hmod_mat_t A, const hmod_mat_t B);

void
_hmod_mat_mul_classical(hmod_mat_t D, const hmod_mat_t C,
                                const hmod_mat_t A, const hmod_mat_t B, int op);

void hmod_mat_addmul(hmod_mat_t D, const hmod_mat_t C,
                                const hmod_mat_t A, const hmod_mat_t B);

void hmod_mat_submul(hmod_mat_t D, const hmod_mat_t C,
                                const hmod_mat_t A, const hmod_mat_t B);

/* Trace */

hlimb_t hmod_mat_trace(const hmod_mat_t mat);

/* Determinant */

hlimb_t _hmod_mat_det(hmod_mat_t A);
hlimb_t hmod_mat_det(const hmod_mat_t A);

/* Rank */

long hmod_mat_rank(const hmod_mat_t A);

/* Inverse */

int hmod_mat_inv(hmod_mat_t B, const hmod_mat_t A);

/* Triangular solving */

void hmod_mat_solve_tril(hmod_mat_t X, const hmod_mat_t L, const hmod_mat_t B, int unit);
void hmod_mat_solve_tril_recursive(hmod_mat_t X, const hmod_mat_t L, const hmod_mat_t B, int unit);
void hmod_mat_solve_tril_classical(hmod_mat_t X, const hmod_mat_t L, const hmod_mat_t B, int unit);

void hmod_mat_solve_triu(hmod_mat_t X, const hmod_mat_t U, const hmod_mat_t B, int unit);
void hmod_mat_solve_triu_recursive(hmod_mat_t X, const hmod_mat_t U, const hmod_mat_t B, int unit);
void hmod_mat_solve_triu_classical(hmod_mat_t X, const hmod_mat_t U, const hmod_mat_t B, int unit);

/* LU decomposition */

long hmod_mat_lu(long * P, hmod_mat_t A, int rank_check);
long hmod_mat_lu_classical(long * P, hmod_mat_t A, int rank_check);
long hmod_mat_lu_recursive(long * P, hmod_mat_t A, int rank_check);

/* Nonsingular solving */

int hmod_mat_solve(hmod_mat_t X, const hmod_mat_t A, const hmod_mat_t B);
int hmod_mat_solve_vec(hlimb_t * x, const hmod_mat_t A, const hlimb_t * b);

/* Reduced row echelon form */

long hmod_mat_rref(hmod_mat_t A);

/* Nullspace */

long hmod_mat_nullspace(hmod_mat_t X, const hmod_mat_t A);


/* Tuning parameters *********************************************************/

/* Size at which pre-transposing becomes faster in classical multiplication */
#define HMOD_MAT_MUL_TRANSPOSE_CUTOFF 20

/* Strassen multiplication */
#define HMOD_MAT_MUL_STRASSEN_CUTOFF 256

/* Cutoff between classical and recursive triangular solving */
#define HMOD_MAT_SOLVE_TRI_ROWS_CUTOFF 64
#define HMOD_MAT_SOLVE_TRI_COLS_CUTOFF 64

/* Cutoff between classical and recursive LU decomposition */
#define HMOD_MAT_LU_RECURSIVE_CUTOFF 4

#ifdef __cplusplus
}
#endif

#endif

