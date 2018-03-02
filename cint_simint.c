/*
 * Copyright (c) 2013-2018 Georgia Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * The GNU Lesser General Public License is included in this distribution
 * in the file COPYING.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <simint/simint.h>
#include <CInt.h>
#include <cint_basisset.h>
#include <cint_config.h>
#include <cint_def.h>

typedef unsigned long long TimerType;
#define CLOCK(ticks) do { \
    volatile unsigned int a__, d__; \
    __asm__ __volatile__("rdtsc" : "=a" (a__), "=d" (d__) : );   \
    (ticks) = ((TimerType) a__)|(((TimerType)d__)<<32); \
  } while(0)

struct SIMINT
{
    int nthreads;
    int max_am;
    int workmem_per_thread;
    int outmem_per_thread;
    double *workbuf;
    double *outbuf;

    int nshells;
    struct simint_shell *shells;
    struct simint_multi_shellpair *shellpairs;

    int    screen_method;
    double screen_tol;

    // only master thread will write to these
    TimerType time_outer;
    TimerType time_inner;
};

// CInt_createSIMINT is called by all nodes.
// All nodes have a copy of the BasisSet_t structure here and will form and 
//   store the Simint shells for all shells of the molecule.

CIntStatus_t CInt_createSIMINT(BasisSet_t basis, SIMINT_t *simint, int nthreads)
{
    CINT_ASSERT(nthreads > 0);

    SIMINT_t s = (SIMINT_t) calloc(1, sizeof(struct SIMINT));
    CINT_ASSERT(s != NULL);

    simint_init();

    // allocate workbuf for all threads on this node
    s->nthreads = nthreads;
    s->max_am = basis->max_momentum;
    s->workmem_per_thread = simint_ostei_workmem(0, s->max_am); // consider aligning
    s->workbuf = (double *) malloc(s->workmem_per_thread*nthreads*sizeof(double));
    CINT_ASSERT(s->workbuf != NULL);

    // allocate outbuf for all threads on this node
    int max_ncart = ( (s->max_am+1)*(s->max_am+2) )/2;
    int maxsize = max_ncart * max_ncart * max_ncart * max_ncart; // consider aligning
    s->outmem_per_thread = maxsize * _SIMINT_NSHELL_SIMD;
    s->outbuf = (double *) malloc(s->outmem_per_thread * nthreads * sizeof(double));
    CINT_ASSERT(s->outbuf != NULL);

    // form and store simint shells for all shells of this molecule
    s->nshells = basis->nshells;
    s->shells = (struct simint_shell *) malloc(sizeof(struct simint_shell)*basis->nshells);
    CINT_ASSERT(s->shells != NULL);

    struct simint_shell *shell_p = s->shells;
    for (int i=0; i<basis->nshells; i++)
    {
        // initialize variables in structure
        simint_initialize_shell(shell_p); 

        // allocate space for alpha and coef for the shell
        simint_allocate_shell(basis->nexp[i], shell_p);

        shell_p->am    = basis->momentum[i];
        shell_p->nprim = basis->nexp[i];
        shell_p->x     = basis->xyz0[i*4+0];
        shell_p->y     = basis->xyz0[i*4+1];
        shell_p->z     = basis->xyz0[i*4+2];

        for (int j=0; j<basis->nexp[i]; j++)
        {
            shell_p->alpha[j] = basis->exp[i][j];
            shell_p->coef[j]  = basis->cc[i][j];
        }

        shell_p++;
    }

    // here we assume there are no unit shells (shells with zero orbital exponent)
    simint_normalize_shells(basis->nshells, s->shells);

    s->screen_method = SIMINT_SCREEN_FASTSCHWARZ;
    s->screen_tol = 1.e-14;
    printf("Screen method: %d\n", s->screen_method);
    printf("Screen tol:    %e\n", s->screen_tol);

    // precompute all shell pairs
    // could just do this as needed
    s->shellpairs = (struct simint_multi_shellpair *)
        malloc(sizeof(struct simint_multi_shellpair)*basis->nshells*basis->nshells);
    CINT_ASSERT(s->shellpairs != NULL);
    // UNDONE: exploit symmetry
    for (int i=0; i<basis->nshells; i++)
    {
        for (int j=0; j<basis->nshells; j++)
        {
            struct simint_multi_shellpair *pair;
            pair = &s->shellpairs[i*basis->nshells+j];
            simint_initialize_multi_shellpair(pair);
            simint_create_multi_shellpair(1, s->shells+i, 1, s->shells+j, pair, s->screen_method);
        }
    }

    s->time_outer = 0;
    s->time_inner = 0;

    *simint = s;
    return CINT_STATUS_SUCCESS;
}

CIntStatus_t CInt_destroySIMINT(SIMINT_t simint)
{
    printf("Time outer: %f\n", simint->time_outer/2.3e9);
    printf("Time inner: %f\n", simint->time_inner/2.3e9);

    struct simint_multi_shellpair *shellpair_p = simint->shellpairs;
    for (int i=0; i<simint->nshells*simint->nshells; i++)
        simint_free_multi_shellpair(shellpair_p++);

    struct simint_shell *shell_p = simint->shells;
    for (int i=0; i<simint->nshells; i++)
        simint_free_shell(shell_p++);

    free(simint->shellpairs);
    free(simint->shells);
    free(simint->workbuf);
    free(simint->outbuf);
    free(simint);

    simint_finalize();
    return CINT_STATUS_SUCCESS;
}

// for Simint, caller provides memory where integrals will be stored;
// for ERD, library returns pointer to where integrals are stored;

/* ---------- huangh223 modification part start ---------- */

int CInt_SIMINT_getShellpairAMIndex(SIMINT_t simint, int P, int Q)
{
    struct simint_shell *shells = simint->shells;
    return shells[P].am * ((_SIMINT_OSTEI_MAXAM) + 1) + shells[Q].am;
}

// Compute ( M N | P_list[i] Q_list[i] ), i = 0, ..., npairs - 1
// AM(P_list[:]) are the same, and AM(Q_list[:]) are the same
// (P_list[:], Q_list[:]) should be packed as a simint_multi_shelpair,
// but I haven't finish that yet.
CIntStatus_t 
CInt_computeShellQuartetBatch_SIMINT(
    BasisSet_t basis, SIMINT_t simint, int tid,
    int M, int N, int *P_list, int *Q_list,
    int npairs, dbl_ptr *batch_integrals, int *batch_nints
)
{
    TimerType start0, stop0;
    TimerType start1, stop1;

    int size, ret;
    struct simint_multi_shellpair *bra_pair_p;
    struct simint_multi_shellpair *ket_pair_p;

    if (tid==0) CLOCK(start0);

    bra_pair_p = &simint->shellpairs[M * basis->nshells + N];
    int simint_max_outbuf_size = simint->outmem_per_thread / _SIMINT_NSHELL_SIMD;
    
    for (int ipair = 0; ipair < npairs; ipair++)
    {
        int P = P_list[ipair];
        int Q = Q_list[ipair];
        ket_pair_p = &simint->shellpairs[P * basis->nshells + Q];

        if (tid==0) CLOCK(start1);
        ret = simint_compute_eri(
            bra_pair_p, ket_pair_p, simint->screen_tol,
            &simint->workbuf[tid*simint->workmem_per_thread],
            &simint->outbuf [tid*simint->outmem_per_thread + ipair * simint_max_outbuf_size]
        );
        if (tid==0) CLOCK(stop1);
        
        if (ret < 0)
        {
            size = 0; // return zero size to caller; output buffer is not initialized
        } else {
            CINT_ASSERT(ret == 1); // single shell quartet
            struct simint_shell *shells = simint->shells;
            size = (shells[M].am+1)*(shells[M].am+2)/2 *
                   (shells[N].am+1)*(shells[N].am+2)/2 *
                   (shells[P].am+1)*(shells[P].am+2)/2 *
                   (shells[Q].am+1)*(shells[Q].am+2)/2;
        }

        batch_integrals[ipair] = &simint->outbuf[tid*simint->outmem_per_thread + ipair * simint_max_outbuf_size];
        batch_nints[ipair]     = size;

        if (tid==0) CLOCK(stop0);
        if (tid==0)
        {
            simint->time_outer += (stop0-start0);
            simint->time_inner += (stop1-start1);
        }
    }

    return CINT_STATUS_SUCCESS;
}

/* ----------- huangh223 modification part end ----------- */

CIntStatus_t 
CInt_computeShellQuartet_SIMINT(BasisSet_t basis, SIMINT_t simint, int tid,
                                int A, int B, int C, int D,
                                double **integrals, int *nints)
{
    TimerType start0, stop0;
    TimerType start1, stop1;

    int size, ret;
    struct simint_multi_shellpair *bra_pair_p;
    struct simint_multi_shellpair *ket_pair_p;

    if (tid==0) CLOCK(start0);

    bra_pair_p = &simint->shellpairs[A*basis->nshells+B];
    ket_pair_p = &simint->shellpairs[C*basis->nshells+D];

    if (tid==0) CLOCK(start1);
    ret = simint_compute_eri(bra_pair_p, ket_pair_p, simint->screen_tol,
      &simint->workbuf[tid*simint->workmem_per_thread],
      &simint->outbuf [tid*simint->outmem_per_thread]);
    if (tid==0) CLOCK(stop1);
    if (ret < 0)
        size = 0; // return zero size to caller; output buffer is not initialized
    else
    {
        CINT_ASSERT(ret == 1); // single shell quartet
        struct simint_shell *shells = simint->shells;
        size = (shells[A].am+1)*(shells[A].am+2)/2 *
               (shells[B].am+1)*(shells[B].am+2)/2 *
               (shells[C].am+1)*(shells[C].am+2)/2 *
               (shells[D].am+1)*(shells[D].am+2)/2;
    }

    *integrals = &simint->outbuf[tid*simint->outmem_per_thread];
    *nints = size;

    if (tid==0) CLOCK(stop0);
    if (tid==0)
    {
        simint->time_outer += (stop0-start0);
        simint->time_inner += (stop1-start1);
    }

    return CINT_STATUS_SUCCESS;
}

// interface has tid argument, different from OED version.
CIntStatus_t
CInt_computePairOvl_SIMINT(BasisSet_t basis, SIMINT_t simint, int tid,
                           int A, int B,
                           double **integrals, int *nints)
{
    int size, ret;
    struct simint_shell *shells = simint->shells;

    ret = simint_compute_overlap(&shells[A], &shells[B],
       &simint->outbuf[tid*simint->outmem_per_thread]);
    CINT_ASSERT(ret == 1);
    size = (shells[A].am+1)*(shells[A].am+2)/2 *
           (shells[B].am+1)*(shells[B].am+2)/2;

    *integrals = &simint->outbuf[tid*simint->outmem_per_thread];
    *nints = size;

    return CINT_STATUS_SUCCESS;
}

// interface has tid argument, different from OED version.
CIntStatus_t
CInt_computePairCoreH_SIMINT(BasisSet_t basis, SIMINT_t simint, int tid,
                           int A, int B,
                           double **integrals, int *nints)
{
    int size, ret;
    struct simint_shell *shells = simint->shells;

    // number of scalar quantities computed
    size = (shells[A].am+1)*(shells[A].am+2)/2 *
           (shells[B].am+1)*(shells[B].am+2)/2;

    // allocate temporary buffer
    double *temp = (double *) malloc(size*sizeof(double));
    CINT_ASSERT(temp != NULL);

    ret = simint_compute_ke(&shells[A], &shells[B], temp);
    CINT_ASSERT(ret == 1);

    ret = simint_compute_potential(basis->natoms, basis->charge,
       basis->xn, basis->yn, basis->zn,
       &shells[A], &shells[B],
       &simint->outbuf[tid*simint->outmem_per_thread]);
    CINT_ASSERT(ret == 1);

    *integrals = &simint->outbuf[tid*simint->outmem_per_thread];
    *nints = size;

    // sum outputs
    double *p = *integrals;
    for (int i=0; i<size; i++)
        *p++ += temp[i];

    free(temp);

    return CINT_STATUS_SUCCESS;
}
