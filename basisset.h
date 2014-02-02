#ifndef __BASISSET_H__
#define __BASISSET_H__


#include "cint_def.h"

__attribute__((target(mic))) void _maxMomentum (BasisSet_t basis, int *max_momentum);

__attribute__((target(mic))) void _maxPrimid (BasisSet_t basis, int *max_primid);

void _maxnumExp (BasisSet_t basis, int *max_nexp);

#endif /* __BASISSET_H__ */
