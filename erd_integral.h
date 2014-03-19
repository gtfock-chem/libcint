#pragma once
#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

#define ERD_SCREEN true
#define ERD_SPHERIC 1


#define MAX(a,b)    ((a) < (b) ? (b) : (a))


#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

extern void erd__1111_csgto(
    uint32_t npgto1, uint32_t npgto2, uint32_t npgto3, uint32_t npgto4,
    uint32_t shell1, uint32_t shell2, uint32_t shell3, uint32_t shell4,
    bool atomic,
    double x1, double y1, double z1,
    double x2, double y2, double z2,
    double x3, double y3, double z3,
    double x4, double y4, double z4,
    const double alpha1[restrict static npgto1], const double alpha2[restrict static npgto2], const double alpha3[restrict static npgto3], const double alpha4[restrict static npgto4],
    const double cc1[restrict static npgto1], const double cc2[restrict static npgto2], const double cc3[restrict static npgto3], const double cc4[restrict static npgto4],
    const double norm1[restrict static npgto1], const double norm2[restrict static npgto2], const double norm3[restrict static npgto3], const double norm4[restrict static npgto4],
    uint32_t integrals_count[restrict static 1], double integrals_ptr[restrict static 81]);

extern void erd__csgto(
    bool atomic,
    uint32_t npgto1, uint32_t npgto2, uint32_t npgto3, uint32_t npgto4,
    uint32_t shell1, uint32_t shell2, uint32_t shell3, uint32_t shell4,
    double x1, double y1, double z1,
    double x2, double y2, double z2,
    double x3, double y3, double z3,
    double x4, double y4, double z4,
    const double alpha1[restrict static npgto1], const double alpha2[restrict static npgto2], const double alpha3[restrict static npgto3], const double alpha4[restrict static npgto4],
    const double cc1[restrict static npgto1], const double cc2[restrict static npgto2], const double cc3[restrict static npgto3], const double cc4[restrict static npgto4],
    const double norm1[restrict static npgto1], const double norm2[restrict static npgto2], const double norm3[restrict static npgto3], const double norm4[restrict static npgto4],
    int **vrrtab, int ldvrrtab,
    bool spheric,
    uint32_t buffer_capacity, uint32_t output_length[restrict static 1], double output_buffer[restrict static 1]);

extern void erd__memory_csgto(uint32_t npgto1, uint32_t npgto2, uint32_t npgto3, uint32_t npgto4,
    uint32_t shell1, uint32_t shell2, uint32_t shell3, uint32_t shell4,
    double x1, double y1, double z1,
    double x2, double y2, double z2,
    double x3, double y3, double z3,
    double x4, double y4, double z4,
    bool spheric, size_t *iopt, size_t *zopt);

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif
