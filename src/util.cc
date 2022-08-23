#include "util.h"

#include <cmath>

void eye (double* m, int const n)
{
  for (int iy = 0; iy != n; ++iy) {
    for (int ix = 0; ix != n; ++ix) {
      if (ix == iy) {
        m[iy*n + ix] = 1;
      } else {
        m[iy*n + ix] = 0;
      }
    }
  }
  return;
}



void VeqMdotV (double *m, double *v, int const N)
{
  // v = m .dot. v
  double w[N] = { 0 };
  for (int iy = 0; iy != N; ++iy) {
    for (int ix = 0; ix != N; ++ix) {
      w[iy] += m[iy*N + ix] * v[ix];
    }
  }
  for (int i = 0; i != N; ++i) {
    v[i] = w[i];
  }
  return;
}


void MeqMdotN (double *m, double *n)
{
  double v[36] = { 0 };
  for (int i = 0; i != 6; ++i) {
    for (int j = 0; j != 6; ++j) {
      for (int k = 0; k != 6; ++k) {
        v[i*6 + j] += m[i*6 + k] * n[k*6 + j];
      }
    }
  }
  for (int i = 0; i != 6; ++i) {
    for (int j = 0; j != 6; ++j) {
      m[i*6 + j] = v[i*6 + j];
    }
  }
  return;
}

void MeqNdotM (double *m, double *n)
{
  double v[36] = { 0 };
  for (int i = 0; i != 6; ++i) {
    for (int j = 0; j != 6; ++j) {
      for (int k = 0; k != 6; ++k) {
        v[i*6 + j] += n[i*6 + k] * m[k*6 + j];
      }
    }
  }
  for (int i = 0; i != 6; ++i) {
    for (int j = 0; j != 6; ++j) {
      m[i*6 + j] = v[i*6 + j];
    }
  }
  return;
}

void LeqMdotN (double *l, double *m, double *n)
{
  for (int i = 0; i != 6; ++i) {
    for (int j = 0; j != 6; ++j) {
      l[i*6 + j] = 0;
      for (int k = 0; k != 6; ++k) {
        l[i*6 + j] += m[i*6 + k] * n[k*6 + j];
      }
    }
  }
  return;
}

void rotmat (double* m, double const angle)
{
  double const c = cos(angle);
  double const s = sin(angle);

  double a[36] = { 0 };
  double n[36] = { 0 };

  for (int iy = 0; iy != 6; ++iy) {
    for (int ix = 0; ix != 6; ++ix) {
      if (ix == iy) {
        if (ix < 4) {
          a[iy*6 + ix] = c;
        } else {
          a[iy*6 + ix] = 1;
        }
      }
    }
  }
  a[0*6+2] = -s;
  a[1*6+3] = -s;
  a[2*6+0] = s;
  a[3*6+1] = s;

  for (int iy = 0; iy != 6; ++iy) {
    for (int ix = 0; ix != 6; ++ix) {
      for (int i = 0; i != 6; ++i) {
        n[iy*6+ix] += a[iy*6+i] * m[i*6+ix];
      }
    }
  }
  a[0*6+2] = s;
  a[1*6+3] = s;
  a[2*6+0] = -s;
  a[3*6+1] = -s;

  for (int iy = 0; iy != 6; ++iy) {
    for (int ix = 0; ix != 6; ++ix) {
      m[iy*6+ix] = 0;
      for (int i = 0; i != 6; ++i) {
        m[iy*6+ix] += n[iy*6+i] * a[i*6+ix];
      }
    }
  }

  return;
}

