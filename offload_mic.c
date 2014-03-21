#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <assert.h>
#include <sys/time.h>
#include <string.h>
#include <ctype.h>
#ifdef __INTEL_OFFLOAD
#include <offload.h>
#endif

#include "config.h"
#include "basisset.h"


#ifdef __INTEL_OFFLOAD
__declspec (target (mic)) ERD_t erd_mic;
__declspec (target (mic)) BasisSet_t basis_mic;
#endif


#ifdef __INTEL_OFFLOAD

CIntStatus_t CInt_offload_createBasisSet (BasisSet_t * _basis)
{
    CIntStatus_t status;
    int i;
    int mic_numdevs;

    status = CInt_createBasisSet (_basis);
    if (status != CINT_STATUS_SUCCESS)
    {
        return status;
    }


    mic_numdevs = _Offload_number_of_devices ();
    _basis[0]->mic_numdevs = mic_numdevs;
    for (i = 0; i < mic_numdevs; i++)
    {
#pragma offload target(mic: i) \
                nocopy(basis_mic) out(status)
        {
            status = CInt_createBasisSet (&basis_mic);
        }
        if (status != CINT_STATUS_SUCCESS)
        {
            return CINT_STATUS_OFFLOAD_ERROR;
        }
    }

    return CINT_STATUS_SUCCESS;
}


CIntStatus_t CInt_offload_destroyBasisSet (BasisSet_t basis)
{
    CIntStatus_t status;
    int i;

    for (i = 0; i < basis->mic_numdevs; i++)
    {
#pragma offload target(mic: i)\
                nocopy(basis_mic) out(status)
        status = CInt_destroyBasisSet (basis_mic);
        if (status != CINT_STATUS_SUCCESS)
        {
            return status;
        }
    }

    status = CInt_destroyBasisSet (basis);
    if (status != CINT_STATUS_SUCCESS)
    {
        return status;
    }

    return CINT_STATUS_SUCCESS;
}


CIntStatus_t CInt_offload_loadBasisSet (BasisSet_t basis,
                                        char *bsfile, char *molfile)
{
    CIntStatus_t status;
    int i;
    char *buf;
    int bufsize;

    // read xyz file
    if ((status = import_molecule (molfile, basis)) != CINT_STATUS_SUCCESS)
    {
        return status;
    }
    // read basis set
    if ((status = import_basis (bsfile, basis)) != CINT_STATUS_SUCCESS)
    {
        return status;
    }

    //parse xyz
    if ((status = parse_molecule (basis)) != CINT_STATUS_SUCCESS)
    {
        return status;
    }

    CInt_packBasisSet (basis, (void **) &buf, &bufsize);
    for (i = 0; i < basis->mic_numdevs; i++)
    {
#pragma offload target(mic: i)\
                in(bufsize) in(buf: length(bufsize))\
                nocopy(basis_mic) out(status)
        {
            status = CInt_unpackBasisSet (basis_mic, buf);
        }
        if (status != CINT_STATUS_SUCCESS)
        {
            return status;
        }
    }
    free (buf);

    return CINT_STATUS_SUCCESS;
}


CIntStatus_t CInt_offload_createERD (BasisSet_t basis, ERD_t * erd,
                                     int nthreads, int nthreads_mic)
{
    CIntStatus_t status;
    int mic_id;
    printf ("create %d %d\n", nthreads, nthreads_mic);
    status = CInt_createERD (basis, erd, nthreads);
    if (status != CINT_STATUS_SUCCESS)
    {
        return status;
    }
    erd[0]->mic_numdevs = basis->mic_numdevs;

    for (mic_id = 0; mic_id < erd[0]->mic_numdevs; mic_id++)
    {
#pragma offload target(mic:mic_id)\
                nocopy(basis_mic, erd_mic)\
                in(nthreads_mic)\
                out(status)
        {
            status = CInt_createERD (basis_mic, &erd_mic, nthreads_mic);
        }
        if (status != CINT_STATUS_SUCCESS)
        {
            return status;
        }
    }

    return CINT_STATUS_SUCCESS;
}


CIntStatus_t CInt_offload_destroyERD (ERD_t erd)
{
    CIntStatus_t status;
    int mic_id;

    status = CInt_destroyERD (erd);
    if (status != CINT_STATUS_SUCCESS)
    {
        return status;
    }

    for (mic_id = 0; mic_id < erd->mic_numdevs; mic_id++)
    {
#pragma offload target(mic:mic_id)\
                nocopy(erd_mic)\
                out(status)
        {
            status = CInt_destroyERD (erd_mic);
        }
        if (status != CINT_STATUS_SUCCESS)
        {
            return status;
        }
    }

    return CINT_STATUS_SUCCESS;
}

#endif /* #ifdef __INTEL_OFFLOAD */
