#ifndef __GUARD_util_h
#define __GUARD_util_h


void eye (double* m, int const n = 6);
void VeqMdotV (double *m, double *v, int const N = 6);
void MeqMdotN (double *m, double *n);
void MeqNdotM (double *m, double *n);
void LeqMdotN (double *l, double *m, double *n);
void rotmat (double* m, double const angle);

#endif
