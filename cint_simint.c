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
//#include <CInt.h>
//#include <cint_type.h>    // BasisSet definition
#include <cint_basisset.h>
#include <cint_config.h>  // CINT_PRINTF,...
#include <cint_def.h>

struct SIMINT
{
    int nthreads;
    int max_am;
    int workmem_per_thread;
    double *workspace;

    struct simint_shell *shells;
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

    // allocate workspace for all threads on this node
    s->nthreads = nthreads;
    s->max_am = basis->max_momentum;
    s->workmem_per_thread = simint_ostei_workmem(0, s->max_am);
    s->workspace = (double *) malloc(s->workmem_per_thread*nthreads*sizeof(double));

    // form and store simint shells for all shells of this molecule
    s->shells = (struct simint_shell *) malloc(sizeof(struct simint_shell)*basis->nshells);
    CINT_ASSERT(s->shells != NULL);

    struct simint_shell *shell_p = s->shells;
    for (int i=0; i<basis->nshells; i++)
    {
        // initialize variables in structure
        simint_initialize_shell(shell_p); 

        // do not allocate space for alpha and coef for the shell;
        // we will set the pointers to existing data;
        // simint_allocate_shell(nprim, shell_p);

        shell_p->am    = basis->momentum[i];
        shell_p->nprim = basis->nexp[i];
        shell_p->x     = basis->xyz0[i*4+0];
        shell_p->y     = basis->xyz0[i*4+1];
        shell_p->z     = basis->xyz0[i*4+2];

        shell_p->alpha = basis->bs_exp[i];
        shell_p->coef  = basis->bs_cc[i];

        // UNDONE: have the shells already been normalized?
        // we should use simint normalization....
        // so we need to make sure shells are not normalized beforehand

        shell_p++;
    }

    // here we assume there are no unit shells (shells with zero orbital exponent)
    simint_normalize_shells(basis->nshells, s->shells);

    *simint = s;
    return CINT_STATUS_SUCCESS;
}

CIntStatus_t CInt_destroySIMINT(SIMINT_t simint)
{
    free(simint->workspace);
    free(simint->shells);
    free(simint);

    simint_finalize();
    return CINT_STATUS_SUCCESS;
}

CIntStatus_t 
CInt_computeShellQuartet_SIMINT(BasisSet_t basis, SIMINT_t simint, int tid,
                                int A, int B, int C, int D,
                                double **integrals, int *nints)
{
    struct simint_shell *shells = simint->shells;

    // multishell pairs
    struct simint_multi_shellpair bra_pair;
    struct simint_multi_shellpair ket_pair;

    simint_initialize_multi_shellpair(&bra_pair);
    simint_initialize_multi_shellpair(&ket_pair);

    // final argument is screen_method
    simint_create_multi_shellpair(1, &shells[A], 1, &shells[B], &bra_pair, 0);
    simint_create_multi_shellpair(1, &shells[C], 1, &shells[D], &ket_pair, 0);

    // set workspace pointer
    // set targets...

    *nints =
        simint_compute_eri(&bra_pair, &ket_pair, 0.0, 
          &simint->workspace[tid*simint->workmem_per_thread], *integrals);

    simint_free_multi_shellpair(&bra_pair);
    simint_free_multi_shellpair(&ket_pair);

    *nints = 1; // UNDONE: meaning of return value is in shell quartets?
    return CINT_STATUS_SUCCESS;
}
