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

#pragma once

#include <stdint.h>
#include <stdbool.h>

struct OED
{
    int nalpha;
    int ncoeff;
    int ncgto1;
    int ncgto2;
    int npgto1;
    int npgto2;
    int shell1;
    int shell2;
    int natoms;
    int ncsum;
    int spheric;
    int screen;  
    double x1;
    double y1;
    double z1;
    double x2;
    double y2;
    double z2;
    double *xn;
    double *yn;
    double *zn;
    double *charge;
    double *cc;
    double *alpha;
    int cc_beg[2];
    int cc_end[2];
    int imax;
    int zmax;
    double *zcore;
    double *zcore2;
    int *icore;

    int fp_memory_opt;
    int int_memory_opt;
    int *coef_offset;
    int *exp_offset;
};


struct ERD
{
    /* The number of threads used for computation */
    uint32_t nthreads;
    size_t capacity;
    double **buffer;    
    /* Used for vrrtable */
    int max_shella;
    /* 2D array */
    int **vrrtable;
#ifdef __INTEL_OFFLOAD
    int mic_numdevs;
#endif

    double erd_wtime;
};

struct BasisSet
{
    // molecular information from xyz file

    int natoms;        // number of atoms in molecule
    int *eid;          // atomic numbers
    double *xn;        // x coords
    double *yn;        // y coords
    double *zn;        // z coords
    double *charge;    // double precision version of atomic numbers
    int nelectrons;    // sum of atomic numbers in molecule (not really num electrons)
    double **guess;    // initial guesses for each element in basis set (should not be in this section)
    int Q;             // net charge read from xyz file (not related to nelectrons)

    double ene_nuc;    // nuclear energy (computed)

    // basis set information from gbs file

    int bs_nelements;  // max number of elements supported in basis set
    int bs_natoms;     // number of elements in basis set
    int basistype;     // Cartesian or spherical
    int *bs_eid;       // atomic numbers of elements in basis set
    int *bs_eptr;      // map atomic number to entry in basis set (array of len bs_nelements)
    int *bs_atom_start;// start of element data in arrays of length nshells (array of length natoms+1)
    int bs_nshells;    // number of shells in the basis set (not the molecule)
    int bs_totnexp;    // total number of primitive functions in basis set
    int *bs_nexp;      // number of primitive functions for shell
    double **bs_exp;   // bs_exp[i] = orbital exponents for shell i
    double **bs_cc;    // bs_cc[i]  = contraction coefs for shell i
    double **bs_norm;  // bs_norm[i] = normalization constants for shell i
    int *bs_momentum;  // bs_momentum[i] = angular momentum for shell i
    
    // shell information for each shell in the given molecule

    uint32_t nshells;    // number of shells in given molecule
    uint32_t nfunctions; // number of basis functions for molecule
    uint32_t *f_start_id;// offset for first basis function for each shell
    uint32_t *f_end_id;
    uint32_t *s_start_id;// start of shell info for each atom
    uint32_t *nexp;      // number of primitives for each shell
    double **exp;        // exponents for each shell in molecule
    double *minexp;      // ?
    double **cc;         // contraction coefficients for each shell in molecule
    double **norm;       // ?
    uint32_t *momentum;  // angular momentum for each shell in molecule
    double *xyz0;        // centers for each shell in molecule, stored as linear array

    uint32_t maxdim;     // max number of functions among all shells in molecule
    uint32_t max_momentum;
    uint32_t max_nexp;
    uint32_t max_nexp_id;
    
    char str_buf[512];

#ifdef __INTEL_OFFLOAD
    int mic_numdevs;
#endif
};
