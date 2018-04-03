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

/*
typedef unsigned long long TimerType;
#define CLOCK(ticks) do { \
    volatile unsigned int a__, d__; \
    __asm__ __volatile__("rdtsc" : "=a" (a__), "=d" (d__) : );   \
    (ticks) = ((TimerType) a__)|(((TimerType)d__)<<32); \
  } while(0)
*/

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

    // For timer, only master thread will write to these
    double ostei_actual, ostei_setup, fock_update_F;

    // For statistic
    double *num_multi_shellpairs, *sum_nprim;
    double *num_screened_prim, *num_unscreened_prim, *num_screened_vec, *num_unscreened_vec;
};

// CInt_createSIMINT is called by all nodes.
// All nodes have a copy of the BasisSet_t structure here and will form and 
// store the Simint shells for all shells of the molecule.

// for Simint, caller provides memory where integrals will be stored;
// for ERD, library returns pointer to where integrals are stored;

CIntStatus_t CInt_createSIMINT(BasisSet_t basis, SIMINT_t *simint, int nthreads)
{
    CINT_ASSERT(nthreads > 0);

    SIMINT_t s = (SIMINT_t) calloc(1, sizeof(struct SIMINT));
    CINT_ASSERT(s != NULL);

    simint_init();

    // Allocate workbuf for all threads on this node
    s->nthreads = nthreads;
    s->max_am = basis->max_momentum;
    s->workmem_per_thread = simint_ostei_workmem(0, s->max_am);
    s->workmem_per_thread = (s->workmem_per_thread + 7) / 8 * 8;  // Align to 8 double (64 bytes)
    s->workbuf = (double *) _mm_malloc(s->workmem_per_thread*nthreads*sizeof(double), 64);
    CINT_ASSERT(s->workbuf != NULL);

    // Allocate outbuf for all threads on this node
    int max_ncart = ( (s->max_am+1)*(s->max_am+2) )/2;
    int maxsize = max_ncart * max_ncart * max_ncart * max_ncart;
    maxsize = (maxsize + 7) / 8 * 8;   // Align to 8 double (64 bytes)

    // Output buffer should holds SIMINT_NSHELL_SIMD ERI results
    // +8 for Simint primitive screening statistic info 
    s->outmem_per_thread = maxsize * _SIMINT_NSHELL_SIMD + 8;  
    s->outbuf = (double *) _mm_malloc(s->outmem_per_thread * nthreads * sizeof(double), 64);
    CINT_ASSERT(s->outbuf != NULL);

    // Form and store simint shells for all shells of this molecule
    s->nshells = basis->nshells;
    s->shells = (struct simint_shell *) malloc(sizeof(struct simint_shell)*basis->nshells);
    CINT_ASSERT(s->shells != NULL);

    struct simint_shell *shell_p = s->shells;
    for (int i=0; i<basis->nshells; i++)
    {
        // Initialize variables in structure
        simint_initialize_shell(shell_p); 

        // Allocate space for alpha and coef for the shell
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

    // Here we assume there are no unit shells (shells with zero orbital exponent)
    simint_normalize_shells(basis->nshells, s->shells);

    s->screen_method = SIMINT_SCREEN_FASTSCHWARZ;
    s->screen_tol = 1.e-14;
    printf("Screen method: %d\n", s->screen_method);
    printf("Screen tol:    %e\n", s->screen_tol);

    // Precompute all shell pairs
    // Will be used by CInt_SIMINT_fillMultishellpairByShellList(), DO NOT SKIP it!!!
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
    
    // Reset timer
    s->ostei_setup   = 0.0;
    s->ostei_actual  = 0.0;
    s->fock_update_F = 0.0;

    // Allocate space for statistic info
    int stat_info_size = sizeof(double) * nthreads;
    s->num_multi_shellpairs = (double*) malloc(stat_info_size);
    s->sum_nprim            = (double*) malloc(stat_info_size);
    s->num_screened_prim    = (double*) malloc(stat_info_size);
    s->num_unscreened_prim  = (double*) malloc(stat_info_size);
    s->num_screened_vec     = (double*) malloc(stat_info_size);
    s->num_unscreened_vec   = (double*) malloc(stat_info_size);
    CINT_ASSERT(s->num_multi_shellpairs != NULL && s->sum_nprim           != NULL);
    CINT_ASSERT(s->num_screened_prim    != NULL && s->num_unscreened_prim != NULL);
    CINT_ASSERT(s->num_screened_vec     != NULL && s->num_unscreened_vec  != NULL);
    memset(s->num_multi_shellpairs, 0, stat_info_size);
    memset(s->sum_nprim,            0, stat_info_size);
    memset(s->num_screened_prim,    0, stat_info_size);
    memset(s->num_unscreened_prim,  0, stat_info_size);
    memset(s->num_screened_vec,     0, stat_info_size);
    memset(s->num_unscreened_vec,   0, stat_info_size);
    
    *simint = s;
    return CINT_STATUS_SUCCESS;
}

CIntStatus_t CInt_destroySIMINT(SIMINT_t simint, int show_stat)
{
    // Generate final statistic info
    double sum_msp = 0, sum_nprim = 0;
    double total_prim = 0, unscreened_prim = 0;
    double total_vec  = 0, unscreened_vec  = 0;
    for (int i = 0; i < simint->nthreads; i++)
    {
        sum_msp    += (double) simint->num_multi_shellpairs[i];
        sum_nprim  += (double) simint->sum_nprim[i];
        total_prim      += simint->num_screened_prim[i] + simint->num_unscreened_prim[i];
        unscreened_prim += simint->num_unscreened_prim[i];
        total_vec       += simint->num_screened_vec[i] + simint->num_unscreened_vec[i];
        unscreened_vec  += simint->num_unscreened_vec[i];
    }
    double avg_nprim = sum_nprim / sum_msp;
    double prim_unscreen_ratio = unscreened_prim / total_prim;
    double vec_unscreen_ratio  = unscreened_vec  / total_vec;
    
    // Print timer and statistic info
    if (show_stat)
    {
        printf(
            "Timer: OSTEI setup, OSTEI actual, fock_task update_F = %lf, %lf, %lf sec\n", 
            simint->ostei_setup, simint->ostei_actual, simint->fock_update_F
        );
        printf(
            "Simint statistic: avg. ket-side nprim, prim unscreened ratio, SIMD unscreened ratio = %.1lf, %.1lf %%, %.1lf %%\n",
            avg_nprim, prim_unscreen_ratio * 100.0, vec_unscreen_ratio * 100.0
        );
    }

    // Free shell pair info
    struct simint_multi_shellpair *shellpair_p = simint->shellpairs;
    for (int i=0; i<simint->nshells*simint->nshells; i++)
        simint_free_multi_shellpair(shellpair_p++);

    struct simint_shell *shell_p = simint->shells;
    for (int i=0; i<simint->nshells; i++)
        simint_free_shell(shell_p++);

    // Free memory
    free(simint->shellpairs);
    free(simint->shells);
    _mm_free(simint->workbuf);
    _mm_free(simint->outbuf);
    free(simint->num_multi_shellpairs);
    free(simint->sum_nprim);
    free(simint->num_screened_prim);
    free(simint->num_unscreened_prim);
    free(simint->num_screened_vec);
    free(simint->num_unscreened_vec);
    free(simint);

    simint_finalize();
    return CINT_STATUS_SUCCESS;
}

void CInt_SIMINT_addupdateFtimer(SIMINT_t simint, double sec)
{
    simint->fock_update_F += sec;
}

int CInt_SIMINT_getShellpairAMIndex(SIMINT_t simint, int P, int Q)
{
    struct simint_shell *shells = simint->shells;
    return shells[P].am * ((_SIMINT_OSTEI_MAXAM) + 1) + shells[Q].am;
}

void CInt_SIMINT_createThreadMultishellpair(void **thread_multi_shellpair)
{
    struct simint_multi_shellpair *multi_shellpair;
    multi_shellpair = (struct simint_multi_shellpair *) malloc(sizeof(struct simint_multi_shellpair));
    CINT_ASSERT(multi_shellpair != NULL);
    
    // Need not to worry about memory allocation, it will be handled later
    simint_initialize_multi_shellpair(multi_shellpair);
    
    *thread_multi_shellpair = multi_shellpair;
}

void CInt_SIMINT_freeThreadMultishellpair(void **thread_multi_shellpair)
{
    struct simint_multi_shellpair *multi_shellpair = *thread_multi_shellpair;
    CINT_ASSERT(multi_shellpair != NULL);
    
    simint_free_multi_shellpair(multi_shellpair);
    
    free(multi_shellpair);
}

static void CInt_SIMINT_fillMultishellpairByShellList(
    SIMINT_t simint, int npairs, int *P_list, int *Q_list, 
    struct simint_multi_shellpair *multi_shellpair
)
{
    // Put the original multi_shellpairs corresponding to the shell
    // pairs (P_list[i], Q_list[i]) into the list
    struct simint_multi_shellpair *Pin[_SIMINT_NSHELL_SIMD];
    for (int ipair = 0; ipair < npairs; ipair++)
    {
        int P = P_list[ipair];
        int Q = Q_list[ipair];
        Pin[ipair] = &simint->shellpairs[P * simint->nshells + Q];
    }
    
    // Reset output multi_shellpair and copy from existing multi_shellpairs.
    // simint_cat_multi_shellpair() will check and allocate memory for output
    multi_shellpair->nprim = 0;
    simint_cat_shellpairs(
        npairs, (const struct simint_multi_shellpair **) Pin, 
        multi_shellpair, simint->screen_method
    );
}

// Compute ( M N | P_list[i] Q_list[i] ), i = 0, ..., npairs - 1
// AM(P_list[:]) are the same, and AM(Q_list[:]) are the same
// (P_list[:], Q_list[:]) will be packed as a simint_multi_shelpair.
CIntStatus_t 
CInt_computeShellQuartetBatch_SIMINT(
    SIMINT_t simint, int tid,
    int M, int N, int *P_list, int *Q_list,
    int npairs, double **thread_batch_integrals, int *thread_batch_nints,
    void **thread_multi_shellpair
)
{
    double setup_start, setup_end, ostei_start, ostei_end;
    
    int ret, size;

    if (tid == 0) setup_start = CInt_get_walltime_sec();

    struct simint_multi_shellpair *bra_pair_p = &simint->shellpairs[M * simint->nshells + N];
    
    struct simint_multi_shellpair *multi_shellpair = (struct simint_multi_shellpair *) *thread_multi_shellpair;
    assert(multi_shellpair != NULL);
    
    CInt_SIMINT_fillMultishellpairByShellList(simint, npairs, P_list, Q_list, multi_shellpair);
    
    simint->num_multi_shellpairs[tid] += 1.0;
    simint->sum_nprim[tid] += (double) multi_shellpair->nprim;
    
    if (tid == 0) 
    {
        setup_end   = CInt_get_walltime_sec();
        ostei_start = CInt_get_walltime_sec();
    }
    
    ret = simint_compute_eri(
        bra_pair_p, multi_shellpair, simint->screen_tol,
        &simint->workbuf[tid*simint->workmem_per_thread],
        &simint->outbuf [tid*simint->outmem_per_thread]
    );
    
    if (tid == 0) ostei_end = CInt_get_walltime_sec();
    
    if (ret <= 0)
    {
        size = 0; // Return zero size to caller; output buffer is not initialized
    } else {
        CINT_ASSERT(ret == npairs);
        struct simint_shell *shells = simint->shells;
        int P = P_list[0], Q = Q_list[0];
        size  = (shells[M].am+1)*(shells[M].am+2)/2 *
                (shells[N].am+1)*(shells[N].am+2)/2 *
                (shells[P].am+1)*(shells[P].am+2)/2 *
                (shells[Q].am+1)*(shells[Q].am+2)/2;
    }
    
    // Shells in P_list[] have same AM, shells in Q_list[] have same AM,
    // The result sizes for each quartets are the same
    *thread_batch_integrals = &simint->outbuf[tid*simint->outmem_per_thread];
    *thread_batch_nints     = size;
    
    double *prim_screen_stat_info = *thread_batch_integrals + size * npairs;
    simint->num_unscreened_prim[tid] += prim_screen_stat_info[0];
    simint->num_screened_prim[tid]   += prim_screen_stat_info[1];
    simint->num_unscreened_vec[tid]  += prim_screen_stat_info[2];
    simint->num_screened_vec[tid]    += prim_screen_stat_info[3];
    
    if (tid == 0)
    {
        simint->ostei_setup  += setup_end - setup_start;
        simint->ostei_actual += ostei_end - ostei_start;
    }

    return CINT_STATUS_SUCCESS;
}

// Compute a shell quartet (AB|CD) at a time
CIntStatus_t 
CInt_computeShellQuartet_SIMINT(SIMINT_t simint, int tid,
                                int A, int B, int C, int D,
                                double **integrals, int *nints)
{
    double setup_start, setup_end, ostei_start, ostei_end;
    
    int size, ret;
    struct simint_multi_shellpair *bra_pair_p;
    struct simint_multi_shellpair *ket_pair_p;

    if (tid == 0) setup_start = CInt_get_walltime_sec();

    bra_pair_p = &simint->shellpairs[A*simint->nshells + B];
    ket_pair_p = &simint->shellpairs[C*simint->nshells + D];
    
    simint->num_multi_shellpairs[tid] += 1.0;
    simint->sum_nprim[tid] += (double) ket_pair_p->nprim;

    if (tid == 0) 
    {
        setup_end   = CInt_get_walltime_sec();
        ostei_start = CInt_get_walltime_sec();
    }
    
    ret = simint_compute_eri(
        bra_pair_p, ket_pair_p, simint->screen_tol,
        &simint->workbuf[tid * simint->workmem_per_thread],
        &simint->outbuf [tid * simint->outmem_per_thread]
    );
    
    if (tid == 0) ostei_end = CInt_get_walltime_sec();
    
    if (ret < 0) 
    {
        size = 0; // Return zero size to caller; output buffer is not initialized
    } else {
        CINT_ASSERT(ret == 1); // Single shell quartet
        struct simint_shell *shells = simint->shells;
        size = (shells[A].am+1)*(shells[A].am+2)/2 *
               (shells[B].am+1)*(shells[B].am+2)/2 *
               (shells[C].am+1)*(shells[C].am+2)/2 *
               (shells[D].am+1)*(shells[D].am+2)/2;
    }

    *integrals = &simint->outbuf[tid*simint->outmem_per_thread];
    *nints = size;
    
    double *prim_screen_stat_info = *integrals + size;
    simint->num_unscreened_prim[tid] += prim_screen_stat_info[0];
    simint->num_screened_prim[tid]   += prim_screen_stat_info[1];
    simint->num_unscreened_vec[tid]  += prim_screen_stat_info[2];
    simint->num_screened_vec[tid]    += prim_screen_stat_info[3];

    if (tid == 0)
    {
        simint->ostei_setup  += setup_end - setup_start;
        simint->ostei_actual += ostei_end - ostei_start;
    }

    return CINT_STATUS_SUCCESS;
}

// Interface has tid argument, different from OED version.
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

// Interface has tid argument, different from OED version.
CIntStatus_t
CInt_computePairCoreH_SIMINT(BasisSet_t basis, SIMINT_t simint, int tid,
                           int A, int B,
                           double **integrals, int *nints)
{
    int size, ret;
    struct simint_shell *shells = simint->shells;

    // Number of scalar quantities computed
    size = (shells[A].am+1)*(shells[A].am+2)/2 *
           (shells[B].am+1)*(shells[B].am+2)/2;

    // Allocate temporary buffer
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

    // Sum outputs
    double *p = *integrals;
    for (int i=0; i<size; i++)
        *p++ += temp[i];

    free(temp);

    return CINT_STATUS_SUCCESS;
}
