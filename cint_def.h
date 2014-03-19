#ifndef __CINT_DEF_H__
#define __CINT_DEF_H__

#include "cint_type.h"
#include <stdint.h>

typedef struct OED *OED_t;
typedef struct ERD *ERD_t;
typedef struct BasisSet *BasisSet_t;


typedef enum
{
    CINT_STATUS_SUCCESS = 0,
    CINT_STATUS_NOT_INITIALIZED = 1,
    CINT_STATUS_ALLOC_FAILED = 2,
    CINT_STATUS_INVALID_VALUE = 3,
    CINT_STATUS_EXECUTION_FAILED = 4,
    CINT_STATUS_INTERNAL_ERROR = 5,
    CINT_STATUS_FILEIO_FAILED = 6,
    CINT_STATUS_OFFLOAD_ERROR = 7
} CIntStatus_t;

#ifdef __INTEL_OFFLOAD
extern __declspec(target(mic)) BasisSet_t basis_mic;
#endif

// basisset
// TODO: develop basisset parser

CIntStatus_t CInt_createBasisSet( BasisSet_t *basis );

CIntStatus_t CInt_loadBasisSet( BasisSet_t basis,
                                char *bsfile,
                                char *xyzfile );

CIntStatus_t CInt_destroyBasisSet( BasisSet_t basis );

CIntStatus_t CInt_packBasisSet( BasisSet_t basis,
                                void **buf,
                                int *bufsize );

CIntStatus_t CInt_unpackBasisSet( BasisSet_t basis,
                                void *buf);

int CInt_getNumShells( BasisSet_t basis );

int CInt_getNumFuncs( BasisSet_t basis );

int CInt_getNumAtoms( BasisSet_t basis );

int CInt_getMaxShellDim( BasisSet_t basis );

int CInt_getNumOccOrb( BasisSet_t basis );

int CInt_getFuncEndInd( BasisSet_t basis,
                        int shellid );

int CInt_getAtomStartInd( BasisSet_t basis,
                          int atomid );

// one electron integrals

CIntStatus_t CInt_createOED( BasisSet_t basis,
                             OED_t *oed );

CIntStatus_t CInt_destroyOED( OED_t oed );
    
CIntStatus_t CInt_computePairKin( BasisSet_t basis,
                                  OED_t oed,
                                  int A,
                                  int B,
                                  double **integrals,
                                  int *nints );

CIntStatus_t CInt_computePairOvl( BasisSet_t basis,
                                  OED_t oed,
                                  int A,
                                  int B,
                                  double **integrals,
                                  int *nints );

CIntStatus_t CInt_computePairPot( BasisSet_t basis,
                                  OED_t oed,
                                  int A,
                                  int B,
                                  double **integrals,
                                  int *nints );

CIntStatus_t CInt_computePairCoreH( BasisSet_t basis,
                                    OED_t oed,
                                    int A,
                                    int B,
                                    double **integrals,
                                    int *nints );


// two electron integrals
#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

int CInt_getShellDim( BasisSet_t basis,
                      int shellid );

int CInt_getFuncStartInd( BasisSet_t basis,
                          int shellid );

CIntStatus_t CInt_createERD( BasisSet_t basis,
                             ERD_t *erd,
                             int nthreads );

CIntStatus_t CInt_destroyERD( ERD_t erd );


CIntStatus_t CInt_computeShellQuartet( BasisSet_t basis,
                                       ERD_t erd,
                                       int tid,
                                       int A,
                                       int B,
                                       int C,
                                       int D,
                                       double **integrals,
                                       int *nints );

CIntStatus_t CInt_computeShellQuartets(BasisSet_t basis,
                                       ERD_t erd,
                                       uint32_t threadId,
                                       uint32_t shellIndixA,
                                       const uint32_t*restrict shellIndicesB,
                                       uint32_t shellIndixC,
                                       const uint32_t*restrict shellIndicesD,
                                       uint32_t shellIndicesCount,
                                       double **integrals,
                                       int *integralsCount);

int CInt_getMaxMemory (ERD_t erd);

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif


#endif /* __CINT_DEF_H__ */
