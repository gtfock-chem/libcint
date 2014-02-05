#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "erd_integral.h"
#include "basisset.h"
#include "config.h"
#include "cint_def.h"


#pragma offload_attribute(push, target(mic))

static void erd_max_scratch (BasisSet_t basis, ERD_t erd)
{
    int max_momentum;
    int max_primid;
    int maxnpgto;
    
    _maxMomentum (basis, &max_momentum);
    _maxPrimid (basis, &max_primid);

    maxnpgto = basis->nexp[max_primid];
        
    if (max_momentum < 2)
    {
        erd__memory_1111_csgto (maxnpgto, maxnpgto, maxnpgto, maxnpgto,
                                max_momentum, max_momentum,
                                max_momentum, max_momentum,
                                1.0, 1.0, 1.0, 2.0, 2.0, 2.0,
                                3.0, 3.0, 3.0, 4.0, 4.0, 4.0,
                                &(erd->int_memory_opt),
                                &(erd->fp_memory_opt));
    }
    else
    {
        erd__memory_csgto (maxnpgto, maxnpgto, maxnpgto, maxnpgto,
                           max_momentum, max_momentum,
                           max_momentum, max_momentum,
                           1.0, 1.0, 1.0, 2.0, 2.0, 2.0,
                           3.0, 3.0, 3.0, 4.0, 4.0, 4.0,
                           ERD_SPHERIC,
                           &(erd->int_memory_opt),
                           &(erd->fp_memory_opt));
    }
    
    printf ("max = %d %d\n", erd->int_memory_opt, erd->fp_memory_opt);
}


CIntStatus_t CInt_createERD (BasisSet_t basis, ERD_t *erd)
{
    ERD_t e;
    int tid;
    
    tid = omp_get_thread_num ();
    if (tid == 0)
        printf ("@@@ Optimized ERD code!!\n");
    
    e = (ERD_t)calloc (1, sizeof(struct ERD));
    if (NULL == e)
    {
#ifndef __INTEL_OFFLOAD
        CINT_PRINTF (1, "memory allocation failed\n");
#endif
        return CINT_STATUS_ALLOC_FAILED;
    }

    erd_max_scratch (basis, e);
    e->zcore = (double *)ALIGNED_MALLOC (e->fp_memory_opt * sizeof(double));
    e->icore = (int *)ALIGNED_MALLOC (e->int_memory_opt * sizeof(int));   
    if (NULL == e->zcore || NULL == e->icore)
    {
#ifndef __INTEL_OFFLOAD
        CINT_PRINTF (1, "memory allocation failed\n");
#endif
        return CINT_STATUS_ALLOC_FAILED;
    }
    e->zmax = e->fp_memory_opt;
    
    *erd = e;

    return CINT_STATUS_SUCCESS;
}


CIntStatus_t CInt_destroyERD (ERD_t erd)
{
    ALIGNED_FREE (erd->zcore);
    ALIGNED_FREE (erd->icore);
    free (erd);

    return CINT_STATUS_SUCCESS;
}


__attribute__((target(mic))) CIntStatus_t CInt_computeShellQuartet ( BasisSet_t basis, ERD_t erd,
                                        int A, int B, int C, int D,
                                        double **integrals, int *nints)
{
    int nfirst;
    int shell1;
    int shell2;
    int shell3;
    int shell4;
    int maxshell;

#if ( _DEBUG_LEVEL_ == 3 )
    if (A < 0 || A >= basis->nshells ||
        B < 0 || B >= basis->nshells ||
        C < 0 || C >= basis->nshells ||
        D < 0 || D >= basis->nshells)
    {
#ifndef __INTEL_OFFLOAD
        CINT_PRINTF (1, "invalid shell indices\n");
#endif
        *nints = 0;
        return CINT_STATUS_INVALID_VALUE;
    }
#endif

    shell1 = basis->momentum[A];
    shell2 = basis->momentum[B];
    shell3 = basis->momentum[C];
    shell4 = basis->momentum[D];
    maxshell = MAX(shell1, shell2);
    maxshell = MAX(maxshell, shell3);
    maxshell = MAX(maxshell, shell4);     
    if (maxshell < 2)
    {
	    erd__1111_csgto (erd->zmax, basis->nexp[A], basis->nexp[B],
                         basis->nexp[C], basis->nexp[D],
                         shell1, shell2, shell3, shell4,
                         basis->x[A], basis->y[A], basis->z[A],
                         basis->x[B], basis->y[B], basis->z[B],
                         basis->x[C], basis->y[C], basis->z[C],
                         basis->x[D], basis->y[D], basis->z[D],
                         basis->exp[A], basis->exp[B],
                         basis->exp[C], basis->exp[D],
                         basis->cc[A], basis->cc[B],
                         basis->cc[C], basis->cc[D],
                         ERD_SCREEN, erd->icore,
                         nints, &nfirst, erd->zcore);
    }
    else
    {
	    erd__csgto (erd->zmax, basis->nexp[A], basis->nexp[B],
                    basis->nexp[C], basis->nexp[D],
                    shell1, shell2, shell3, shell4,
                    basis->x[A], basis->y[A], basis->z[A],
                    basis->x[B], basis->y[B], basis->z[B],
                    basis->x[C], basis->y[C], basis->z[C],
                    basis->x[D], basis->y[D], basis->z[D],
                    basis->exp[A], basis->exp[B],
                    basis->exp[C], basis->exp[D],
                    basis->cc[A], basis->cc[B],
                    basis->cc[C], basis->cc[D],
                    ERD_SPHERIC, ERD_SCREEN, erd->icore,
                    nints, &nfirst, erd->zcore);
    }

    *integrals = &(erd->zcore[nfirst - 1]);

    return CINT_STATUS_SUCCESS;
}


int CInt_getMaxMemory (ERD_t erd)
{
    return (erd->fp_memory_opt + erd->int_memory_opt);
}

#pragma offload_attribute(pop)
