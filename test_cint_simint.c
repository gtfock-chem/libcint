#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <unistd.h>
#include <libgen.h>
#include <omp.h>

#include "CInt.h"
#include "simint/simint.h"

int main (int argc, char **argv)
{
    // create basis set    
    BasisSet_t basis;
    CInt_createBasisSet(&basis);

    int nshells;
    int natoms;
    int nfunctions;

    CInt_loadBasisSet(basis, argv[1], argv[2]);
    nshells = CInt_getNumShells(basis);
    natoms = CInt_getNumAtoms(basis);
    nfunctions = CInt_getNumFuncs(basis);
    printf("Job information:\n");
    char *fname;
    fname = basename(argv[2]);
    printf("  molecule:  %s\n", fname);
    fname = basename(argv[1]);
    printf("  basisset:  %s\n", fname);
    printf("  charge     = %d\n", CInt_getTotalCharge(basis));
    printf("  #atoms     = %d\n", natoms);
    printf("  #shells    = %d\n", nshells);
    printf("  #functions = %d\n", nfunctions);
    int nthreads = omp_get_max_threads();
    printf("  #nthreads_cpu = %d\n", nthreads);

    ERD_t erd;
    SIMINT_t simint;
    CInt_createERD(basis, &erd, nthreads);
    CInt_createSIMINT(basis, &simint, nthreads);

    double *integrals;
    int nints;

    for (int i=0; i<nshells; i++)
    for (int j=0; j<nshells; j++)
    for (int k=0; k<1; k++)
    for (int l=0; l<1; l++)
    {
        //CInt_computeShellQuartet(basis, erd, /*tid*/0,
        //    i, j, k, l, &integrals, &nints);
        CInt_computeShellQuartet_SIMINT(basis, simint, /*tid*/0,
            i, j, k, l, &integrals, &nints);
        printf("%d %d %d %d num integrals computed = %d\n", i, j, k, l, nints);
    }
    // integrals cannot be compared because simint uses cartesian functions
    // and ERD uses spheric functions

    CInt_destroyERD(erd);
    CInt_destroySIMINT(simint);
}
