#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "erd_integral.h"
#include "basisset.h"
#include "config.h"
#include "cint_def.h"

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif


static void erd_max_scratch (BasisSet_t basis, ERD_t erd)
{
    int max_momentum;
    int max_primid;
    int maxnpgto;

    max_momentum = basis->max_momentum;
    max_primid = basis->max_nexp_id;
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
}


static CIntStatus_t create_vrrtable (BasisSet_t basis, ERD_t erd)
{
    int xf;
    int xfp;
    int xfmax;
    int yfend;
    int xyfp;
    int yf;
    int xyf;
    int sfend;
    int sf;
    int zf;
    int i;
    int idx;
    int nxyzf;
    int nxyzq;
    int nxyzft;
    int count;
    int max_shella;
    int tablesize;
    int **vrrtable;
    int shellp;
    int shella;
    
    max_shella = basis->max_momentum + 1;
    erd->max_shella = max_shella;
    tablesize = 2 * max_shella * max_shella;        
    vrrtable = (int **)malloc (sizeof(int *) * tablesize);
    if (NULL == vrrtable)
    {
#ifndef __INTEL_OFFLOAD
        CINT_PRINTF (1, "memory allocation failed\n");
#endif
        return CINT_STATUS_ALLOC_FAILED;
    }
    for (shella = 0; shella < max_shella; shella++)
    {
        for (shellp = shella; shellp < 2 * max_shella; shellp++)
        {
            i = shella * 2 * max_shella + shellp;
            nxyzq = (shellp + 1) * (shellp + 2) / 2;
            nxyzft = (shellp + 1) * (shellp + 2) * (shellp + 3) / 6 -
                     shella * (shella + 1) * (shella + 2) / 6;
            
            vrrtable[i] = (int *)ALIGNED_MALLOC (sizeof(int) * 4 * nxyzft);
            if (NULL == vrrtable[i])
            {
                return CINT_STATUS_ALLOC_FAILED;
            }
            // compute tables    
            xfp = nxyzft + 2;
            count = 0;
            for (xf = 0; xf <= shellp; ++xf)
            {
                xfp = xfp + xf - 2;
                xfmax = xf * shellp;
                yfend = shellp - xf;
                xyfp = xfp - xfmax;
                for (yf = 0; yf <= yfend; ++yf)
                {
                    xyf = xf + yf;
                    --xyfp;
                    sfend = MAX (shella, xyf);
                    nxyzf = nxyzq;
                    idx = xyfp;
                    for (sf = shellp; sf >= sfend; --sf)
                    {
                        zf = sf - xyf;                       
                        vrrtable[i][count + 0] = xf;
                        vrrtable[i][count + 1] = yf;
                        vrrtable[i][count + 2] = zf;
                        vrrtable[i][count + 3] = idx;
                        count += 4;
                        idx = idx - nxyzf + xf;
                        nxyzf = nxyzf - sf - 1;
                    }
                }
            }            
        }
    }

    erd->vrrtable = vrrtable;

    return CINT_STATUS_SUCCESS;
}


static CIntStatus_t destroy_vrrtable (ERD_t erd)
{
    int max_shella = erd->max_shella;
    int i;
    int j;
    int idx;
        
    for (i = 0; i < max_shella; i++)
    {
        for (j = i; j < 2 * max_shella; j++)
        {
            idx = i * 2 * max_shella + j;
            ALIGNED_FREE (erd->vrrtable[idx]);
        }
    }
    free (erd->vrrtable);
    return CINT_STATUS_SUCCESS;
}


CIntStatus_t CInt_createERD (BasisSet_t basis, ERD_t *erd, int nthreads)
{      
    ERD_t e;
    CIntStatus_t status;
    int i;

    if (nthreads <= 0)
    {
#ifndef __INTEL_OFFLOAD
        CINT_PRINTF (1, "invalid number of threads\n");
#endif
        return CINT_STATUS_INVALID_VALUE;
    }

    // malloc erd
    e = (ERD_t)calloc (1, sizeof(struct ERD));
    if (NULL == e)
    {
#ifndef __INTEL_OFFLOAD
        CINT_PRINTF (1, "memory allocation failed\n");
#endif
        return CINT_STATUS_ALLOC_FAILED;
    }
    erd_max_scratch (basis, e);

    // memory scratch memory
    e->nthreads = nthreads;
    e->zcore = (double **)malloc (nthreads * sizeof(double *));
    e->icore = (int **)malloc (nthreads * sizeof(int *));   
    if (NULL == e->zcore || NULL == e->icore)
    {
#ifndef __INTEL_OFFLOAD
        CINT_PRINTF (1, "memory allocation failed\n");
#endif
        return CINT_STATUS_ALLOC_FAILED;
    }    
    for (i = 0; i < nthreads; i++)
    {
        e->zcore[i] =
            (double *)ALIGNED_MALLOC (e->fp_memory_opt * sizeof(double));
        e->icore[i] =
            (int *)ALIGNED_MALLOC (e->int_memory_opt * sizeof(int));   
        if (NULL == e->zcore[i] || NULL == e->icore[i])
        {
    #ifndef __INTEL_OFFLOAD
            CINT_PRINTF (1, "memory allocation failed\n");
    #endif
            return CINT_STATUS_ALLOC_FAILED;
        }
    }
     
    // create vrr table
    status = create_vrrtable (basis, e);
    if (status != CINT_STATUS_SUCCESS)
    {
        return status;
    }
    CINT_INFO ("totally use %.3lf MB (%.3lf MB per thread)",
        (e->fp_memory_opt * sizeof(double)
        + e->int_memory_opt * sizeof(int)) * nthreads/1024.0/1024.0,
        (e->fp_memory_opt * sizeof(double)
        + e->int_memory_opt * sizeof(int))/1024.0/1024.0);
    *erd = e;
    return CINT_STATUS_SUCCESS;
}


CIntStatus_t CInt_destroyERD (ERD_t erd)
{
    int i;

    for (i = 0; i < erd->nthreads; i++)
    {
        ALIGNED_FREE (erd->zcore[i]);
        ALIGNED_FREE (erd->icore[i]);
    }
    free (erd->zcore);
    free (erd->icore);
    
    destroy_vrrtable (erd);

    free (erd);

    return CINT_STATUS_SUCCESS;
}


CIntStatus_t CInt_computeShellQuartet ( BasisSet_t basis, ERD_t erd, int tid,
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
    if (tid <= 0 ||
        tid >= erd->nthreads)
    {
#ifndef __INTEL_OFFLOAD
        CINT_PRINTF (1, "invalid thread id\n");
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
	    erd__1111_csgto (erd->fp_memory_opt, basis->nexp[A], basis->nexp[B],
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
                         basis->norm[A], basis->norm[B],
                         basis->norm[C], basis->norm[D],
                         ERD_SCREEN, erd->icore[tid],
                         nints, &nfirst, erd->zcore[tid]);
    }
    else
    {
	    erd__csgto (erd->fp_memory_opt, basis->nexp[A], basis->nexp[B],
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
                    basis->norm[A], basis->norm[B],
                    basis->norm[C], basis->norm[D],
                    erd->vrrtable, 2 * erd->max_shella,
                    ERD_SPHERIC, ERD_SCREEN,
                    erd->icore[tid], nints,
                    &nfirst, erd->zcore[tid]);
    }

    *integrals = &(erd->zcore[tid][nfirst - 1]);

    return CINT_STATUS_SUCCESS;
}


int CInt_getMaxMemory (ERD_t erd)
{
    return (erd->fp_memory_opt + erd->int_memory_opt);
}


#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif
