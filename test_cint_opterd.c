#include <stdio.h>
#include <libgen.h>
#include <omp.h>

#include "CInt.h"

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

int main (int argc, char **argv)
{
    // create basis set    
    BasisSet_t basis;
    CInt_createBasisSet(&basis);

    int nshells;
    int natoms;

    CInt_loadBasisSet(basis, argv[1], argv[2]);
    nshells = CInt_getNumShells(basis);
    natoms = CInt_getNumAtoms(basis);
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

    ERD_t erd;
    CInt_createERD(basis, &erd, nthreads);

    double *integrals;
    int nints;

    // test two electron integrals
    for (int i=0; i<nshells; i++)
    for (int j=0; j<nshells; j++)
    for (int k=0; k<nshells; k++)
    for (int l=0; l<nshells; l++)
    {
        CInt_computeShellQuartet(basis, erd, /*tid*/0,
            i, j, k, l, &integrals, &nints);
        printf("shell quartet %3d %3d %3d %3d\n", i, j, k, l);
        printvec(nints, integrals);
    }

    OED_t *oed = (OED_t *)malloc(sizeof(OED_t) * nthreads);
    for (int i=0; i<nthreads; i++)
        CInt_createOED(basis, &(oed[i]));

    // test one electron integrals
    // OptERD may return nints=0, i.e., values are screened
    for (int i=0; i<nshells; i++)
    for (int j=0; j<nshells; j++)
    {
        CInt_computePairOvl(basis, oed[0], i, j, &integrals, &nints);
        printf("shell pair %3d %3d\n", i, j);
        printvec(nints, integrals);

        CInt_computePairCoreH(basis, oed[0], i, j, &integrals, &nints);
        //printf("%d %d num 1e integrals computed = %d\n", i, j, nints);
    }

    for (int i=0; i<nthreads; i++)
        CInt_destroyOED(oed[i]);
    free(oed);

    CInt_destroyERD(erd);
}
