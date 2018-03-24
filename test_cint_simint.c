#include <stdio.h>
#include <libgen.h>
#include <omp.h>

#include "CInt.h"
#include "simint/simint.h"

#define ABS(x) ((x)<0 ? -(x):(x))
void printvec(int nints, double *integrals)
{
    for (int i=0; i<nints; i++)
    {
        double temp = integrals[i];
        if (ABS(temp) < 1e-15)
            temp = 0.;
        printf("%3d    %e\n", i, temp);
    }
}

// permute 4-index array
void permute(int n1, int n2, int n3, int n4, const double *a, double *b)
{
    for (int i=0; i<n1; i++)
    for (int j=0; j<n2; j++)
    for (int k=0; k<n3; k++)
    for (int l=0; l<n4; l++)
        b[n1*n2*n3*l + n1*n2*k + n1*j + i] = a[n4*n3*n2*i + n4*n3*j + n4*k + l];
}

int main (int argc, char **argv)
{
    // create basis set    
    BasisSet_t basis;
    CInt_createBasisSet(&basis);

    int nshells;
    int natoms;
    int maxdim;

    CInt_loadBasisSet(basis, argv[1], argv[2]);
    nshells = CInt_getNumShells(basis);
    natoms = CInt_getNumAtoms(basis);
    maxdim = CInt_getMaxShellDim(basis);
    printf("Job information:\n");
    char *fname;
    fname = basename(argv[2]);
    printf("  molecule:  %s\n", fname);
    fname = basename(argv[1]);
    printf("  basisset:  %s\n", fname);
    printf("  charge     = %d\n", CInt_getTotalCharge(basis));
    printf("  #atoms     = %d\n", natoms);
    printf("  #shells    = %d\n", nshells);
    int nthreads = omp_get_max_threads();
    printf("  #nthreads_cpu = %d\n", nthreads);
    printf("  max dim = %d\n", maxdim);

    SIMINT_t simint;
    CInt_createSIMINT(basis, &simint, nthreads);

    double buffer[maxdim*maxdim*maxdim*maxdim];
    double *integrals;
    int nints;

    // test two electron integrals
    for (int i=0; i<nshells; i++)
    for (int j=0; j<nshells; j++)
    for (int k=0; k<nshells; k++)
    for (int l=0; l<nshells; l++)
    {
        CInt_computeShellQuartet_SIMINT(simint, /*tid*/0,
            i, j, k, l, &integrals, &nints);
        permute(CInt_getShellDim(basis, i),
                CInt_getShellDim(basis, j),
                CInt_getShellDim(basis, k),
                CInt_getShellDim(basis, l), integrals, buffer);
        printf("shell quartet %3d %3d %3d %3d\n", i, j, k, l);
        printvec(nints, buffer);
    }

    // test one electron integrals
    for (int i=0; i<nshells; i++)
    for (int j=0; j<nshells; j++)
    {
        CInt_computePairOvl_SIMINT(basis, simint, /*tid*/ 0, i, j, &integrals, &nints);

        CInt_computePairCoreH_SIMINT(basis, simint, /*tid*/ 0, i, j, &integrals, &nints);
        permute(CInt_getShellDim(basis, i),
                CInt_getShellDim(basis, j), 1, 1, integrals, buffer);
        printf("shell pair %3d %3d\n", i, j);
        printvec(nints, buffer);
    }

    CInt_destroySIMINT(simint, 1);
}
