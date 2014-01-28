#ifndef __ERD_INTEGRAL_H__
#define __ERD_INTEGRAL_H__


#define ERD_SCREEN  1
#define ERD_SPHERIC 1


#define MAX(a,b)    ((a) < (b) ? (b) : (a))


extern int erd__1111_csgto (int zmax, int npgto1, int npgto2,
                            int npgto3, int npgto4,
                            int shell1, int shell2,
                            int shell3, int shell4,
                            double x1, double y1, double z1,
                            double x2, double y2, double z2,
                            double x3, double y3, double z3,
                            double x4, double y4, double z4,
                            double *alpha, double *cc, int screen,
                            int *icore, int *nbatch, int *nfirst, double *zcore);

extern int erd__csgto (int zmax, int npgto1, int npgto2,
                       int npgto3, int npgto4,
                       int shell1, int shell2,
                       int shell3, int shell4,
                       double x1, double y1, double z1,
                       double x2, double y2, double z2,
                       double x3, double y3, double z3,
                       double x4, double y4, double z4,
                       double *alpha, double *cc,
                       int spheric, int screen, int *icore,
                       int *nbatch, int *nfirst, double *zcore);

extern int erd__memory_1111_csgto (int npgto1, int npgto2,
                                   int npgto3, int npgto4,
                                   int shell1, int shell2,
                                   int shell3, int shell4,
                                   double x1, double y1, double z1,
                                   double x2, double y2, double z2,
                                   double x3, double y3, double z3,
                                   double x4, double y4, double z4,
                                   int *imin, int *iopt,
                                   int *zmin, int *zopt);

extern int erd__memory_csgto (int npgto1, int npgto2,
                              int npgto3, int npgto4,
                              int shell1, int shell2, int shell3, int shell4,
                              double x1, double y1, double z1,
                              double x2, double y2, double z2,
                              double x3, double y3, double z3,
                              double x4, double y4, double z4,
                              int spheric, int *imin, int *iopt,
                              int *zmin, int *zopt);


#endif /* __ERD_INTEGRAL_H__ */
