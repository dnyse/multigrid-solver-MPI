/*
 * Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of nusif-solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 */
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "allocate.h"
#include "comm.h"
#include "solver.h"

#define FINEST_LEVEL 0
#define COARSEST_LEVEL (s->levels - 1)
#define S(i, j, k)                                                             \
  s[(k) * (imaxLvl + 2) * (jmaxLvl + 2) + (j) * (imaxLvl + 2) + (i)]
#define E(i, j, k)                                                             \
  e[(k) * (imaxLvl + 2) * (jmaxLvl + 2) + (j) * (imaxLvl + 2) + (i)]
#define R(i, j, k)                                                             \
  r[(k) * (imaxLvl + 2) * (jmaxLvl + 2) + (j) * (imaxLvl + 2) + (i)]
#define OLD(i, j, k)                                                           \
  old[(k) * (imaxLvl + 2) * (jmaxLvl + 2) + (j) * (imaxLvl + 2) + (i)]
#define P(i, j, k)                                                             \
  p[(k) * (imaxLvl + 2) * (jmaxLvl + 2) + (j) * (imaxLvl + 2) + (i)]
#define RHS(i, j, k)                                                           \
  rhs[(k) * (imaxLvl + 2) * (jmaxLvl + 2) + (j) * (imaxLvl + 2) + (i)]

static void restrictMG(Solver *s, int level, int imaxLvl, int jmaxLvl,
                       int kmaxLvl) {
  double *r = s->r[level + 1];
  double *old = s->r[level];
  int ic = imaxLvl / 2;
  int jc = jmaxLvl / 2;

  commExchangeLevel(s->comm, old, level);

  for (int k = 1; k < (kmaxLvl + 1) / 2; k++) {
    for (int j = 1; j < (jmaxLvl + 1) / 2; j++) {
      for (int i = 1; i < (imaxLvl + 1) / 2; ++i) {
        r[k * (ic + 2) * (jc + 2) + j * (ic + 2) + i] =
            (OLD(2 * i - 1, 2 * j - 1, 2 * k) +
             OLD(2 * i, 2 * j - 1, 2 * k) * 2 +
             OLD(2 * i + 1, 2 * j - 1, 2 * k) +
             OLD(2 * i - 1, 2 * j, 2 * k) * 2 + OLD(2 * i, 2 * j, 2 * k) * 8 +
             OLD(2 * i + 1, 2 * j, 2 * k) * 2 +
             OLD(2 * i - 1, 2 * j + 1, 2 * k) +
             OLD(2 * i, 2 * j + 1, 2 * k) * 2 +
             OLD(2 * i + 1, 2 * j + 1, 2 * k) +

             OLD(2 * i - 1, 2 * j - 1, 2 * k - 1) +
             OLD(2 * i, 2 * j - 1, 2 * k - 1) * 2 +
             OLD(2 * i + 1, 2 * j - 1, 2 * k - 1) +
             OLD(2 * i - 1, 2 * j, 2 * k - 1) * 2 +
             OLD(2 * i, 2 * j, 2 * k - 1) * 4 +
             OLD(2 * i + 1, 2 * j, 2 * k - 1) * 2 +
             OLD(2 * i - 1, 2 * j + 1, 2 * k - 1) +
             OLD(2 * i, 2 * j + 1, 2 * k - 1) * 2 +
             OLD(2 * i + 1, 2 * j + 1, 2 * k - 1) +

             OLD(2 * i - 1, 2 * j - 1, 2 * k + 1) +
             OLD(2 * i, 2 * j - 1, 2 * k + 1) * 2 +
             OLD(2 * i + 1, 2 * j - 1, 2 * k + 1) +
             OLD(2 * i - 1, 2 * j, 2 * k + 1) * 2 +
             OLD(2 * i, 2 * j, 2 * k + 1) * 4 +
             OLD(2 * i + 1, 2 * j, 2 * k + 1) * 2 +
             OLD(2 * i - 1, 2 * j + 1, 2 * k + 1) +
             OLD(2 * i, 2 * j + 1, 2 * k + 1) * 2 +
             OLD(2 * i + 1, 2 * j + 1, 2 * k + 1)) /
            64.0;
      }
    }
  }
}

// static void restrictMG(Solver *s, int level, int imaxLvl, int jmaxLvl,
//                        int kmaxLvl) {
//   // TODO: Exchange ghosts of residual before restricting
//   double *r = s->r[level + 1];
//   double *old = s->r[level];
//
//   commExchangeLevel(s->comm, old, level);
//
//   for (int k = 1; k < (kmaxLvl + 1) / 2; k++) {
//     for (int j = 1; j < (jmaxLvl + 1) / 2; j++) {
//       for (int i = 1; i < (imaxLvl + 1) / 2; ++i) {
//         R(i, j, k) =
//             (OLD(2 * i - 1, 2 * j - 1, 2 * k) +
//              OLD(2 * i, 2 * j - 1, 2 * k) * 2 +
//              OLD(2 * i + 1, 2 * j - 1, 2 * k) +
//              OLD(2 * i - 1, 2 * j, 2 * k) * 2 + OLD(2 * i, 2 * j, 2 * k) * 8
//              + OLD(2 * i + 1, 2 * j, 2 * k) * 2 + OLD(2 * i - 1, 2 * j + 1, 2
//              * k) + OLD(2 * i, 2 * j + 1, 2 * k) * 2 + OLD(2 * i + 1, 2 * j +
//              1, 2 * k) +
//
//              OLD(2 * i - 1, 2 * j - 1, 2 * k - 1) +
//              OLD(2 * i, 2 * j - 1, 2 * k - 1) * 2 +
//              OLD(2 * i + 1, 2 * j - 1, 2 * k - 1) +
//              OLD(2 * i - 1, 2 * j, 2 * k - 1) * 2 +
//              OLD(2 * i, 2 * j, 2 * k - 1) * 4 +
//              OLD(2 * i + 1, 2 * j, 2 * k - 1) * 2 +
//              OLD(2 * i - 1, 2 * j + 1, 2 * k - 1) +
//              OLD(2 * i, 2 * j + 1, 2 * k - 1) * 2 +
//              OLD(2 * i + 1, 2 * j + 1, 2 * k - 1) +
//
//              OLD(2 * i - 1, 2 * j - 1, 2 * k + 1) +
//              OLD(2 * i, 2 * j - 1, 2 * k + 1) * 2 +
//              OLD(2 * i + 1, 2 * j - 1, 2 * k + 1) +
//              OLD(2 * i - 1, 2 * j, 2 * k + 1) * 2 +
//              OLD(2 * i, 2 * j, 2 * k + 1) * 4 +
//              OLD(2 * i + 1, 2 * j, 2 * k + 1) * 2 +
//              OLD(2 * i - 1, 2 * j + 1, 2 * k + 1) +
//              OLD(2 * i, 2 * j + 1, 2 * k + 1) * 2 +
//              OLD(2 * i + 1, 2 * j + 1, 2 * k + 1)) /
//             64.0;
//       }
//     }
//   }
// }
//

static void prolongate(Solver *s, int level, int imaxLvl, int jmaxLvl,
                       int kmaxLvl) {
  double *old = s->e[level + 1];
  double *e = s->e[level];
  int ic = imaxLvl / 2;
  int jc = jmaxLvl / 2;

  for (int k = 2; k < kmaxLvl + 1; k++) {
    for (int j = 2; j < jmaxLvl + 1; j++) {
      for (int i = 2; i < imaxLvl + 1; i++) {
        E(i, j, k) = old[((k + 1) / 2) * (ic + 2) * (jc + 2) +
                         ((j + 1) / 2) * (ic + 2) + ((i + 1) / 2)];
      }
    }
  }
}
// static void prolongate(Solver *s, int level, int imaxLvl, int jmaxLvl,
//                        int kmaxLvl) {
//   // double *old = s->r[level + 1];
//   // double *e = s->r[level];
//   double *old = s->e[level + 1];
//   double *e = s->e[level];
//
//   for (int k = 2; k < kmaxLvl + 1; k++) {
//     for (int j = 2; j < jmaxLvl + 1; j++) {
//       for (int i = 2; i < imaxLvl + 1; i++) {
//         E(i, j, k) = OLD((i + 1) / 2, (j + 1) / 2, (k + 1) / 2);
//       }
//     }
//   }
// }

static void correct(Solver *s, double *p, int level, int imaxLvl, int jmaxLvl,
                    int kmaxLvl) {
  double *e = s->e[level];

  for (int k = 1; k < kmaxLvl + 1; ++k) {
    for (int j = 1; j < jmaxLvl + 1; ++j) {
      for (int i = 1; i < imaxLvl + 1; ++i) {
        P(i, j, k) += E(i, j, k);
      }
    }
  }
}

static void setBoundaryCondition(Solver *s, double *p, int imaxLvl, int jmaxLvl,
                                 int kmaxLvl) {
  if (commIsBoundary(s->comm, FRONT)) {
    for (int j = 1; j < jmaxLvl + 1; j++) {
      for (int i = 1; i < imaxLvl + 1; i++) {
        P(i, j, 0) = P(i, j, 1);
      }
    }
  }

  if (commIsBoundary(s->comm, BACK)) {
    for (int j = 1; j < jmaxLvl + 1; j++) {
      for (int i = 1; i < imaxLvl + 1; i++) {
        P(i, j, kmaxLvl + 1) = P(i, j, kmaxLvl);
      }
    }
  }

  if (commIsBoundary(s->comm, BOTTOM)) {
    for (int k = 1; k < kmaxLvl + 1; k++) {
      for (int i = 1; i < imaxLvl + 1; i++) {
        P(i, 0, k) = P(i, 1, k);
      }
    }
  }

  if (commIsBoundary(s->comm, TOP)) {
    for (int k = 1; k < kmaxLvl + 1; k++) {
      for (int i = 1; i < imaxLvl + 1; i++) {
        P(i, jmaxLvl + 1, k) = P(i, jmaxLvl, k);
      }
    }
  }

  if (commIsBoundary(s->comm, LEFT)) {
    for (int k = 1; k < kmaxLvl + 1; k++) {
      for (int j = 1; j < jmaxLvl + 1; j++) {
        P(0, j, k) = P(1, j, k);
      }
    }
  }

  if (commIsBoundary(s->comm, RIGHT)) {
    for (int k = 1; k < kmaxLvl + 1; k++) {
      for (int j = 1; j < jmaxLvl + 1; j++) {
        P(imaxLvl + 1, j, k) = P(imaxLvl, j, k);
      }
    }
  }
}

static void smooth(Solver *s, double *p, double *rhs, int level, int imaxLvl,
                   int jmaxLvl, int kmaxLvl) {
  double dx2 = s->grid->dx * s->grid->dx;
  double dy2 = s->grid->dy * s->grid->dy;
  double dz2 = s->grid->dz * s->grid->dz;
  double idx2 = 1.0 / dx2;
  double idy2 = 1.0 / dy2;
  double idz2 = 1.0 / dz2;
  double factor =
      s->omega * 0.5 * (dx2 * dy2 * dz2) / (dy2 * dz2 + dx2 * dz2 + dx2 * dy2);
  int pass, ksw, jsw, isw;

  ksw = 1;

  for (pass = 0; pass < 2; pass++) {
    jsw = ksw;

    // TODO: Exchange ghost cells after each red-black pass
    commExchangeLevel(s->comm, p, level);

    for (int k = 1; k < kmaxLvl + 1; k++) {
      isw = jsw;
      for (int j = 1; j < jmaxLvl + 1; j++) {
        for (int i = isw; i < imaxLvl + 1; i += 2) {

          P(i, j, k) -=
              factor *
              (RHS(i, j, k) -
               ((P(i + 1, j, k) - 2.0 * P(i, j, k) + P(i - 1, j, k)) * idx2 +
                (P(i, j + 1, k) - 2.0 * P(i, j, k) + P(i, j - 1, k)) * idy2 +
                (P(i, j, k + 1) - 2.0 * P(i, j, k) + P(i, j, k - 1)) * idz2));
        }
        isw = 3 - isw;
      }
      jsw = 3 - jsw;
    }
    ksw = 3 - ksw;
  }
}

static double calculateResidual(Solver *s, double *p, double *rhs, int level,
                                int imaxLvl, int jmaxLvl, int kmaxLvl) {
  double dx2 = s->grid->dx * s->grid->dx;
  double dy2 = s->grid->dy * s->grid->dy;
  double dz2 = s->grid->dz * s->grid->dz;
  double idx2 = 1.0 / dx2;
  double idy2 = 1.0 / dy2;
  double idz2 = 1.0 / dz2;
  double *r = s->r[level];
  int pass, ksw, jsw, isw;
  double res = 0.0; // WARN: Perhaps 0.0 would make more sense?

  ksw = 1;

  for (pass = 0; pass < 2; pass++) {
    jsw = ksw;

    for (int k = 1; k < kmaxLvl + 1; k++) {
      isw = jsw;
      for (int j = 1; j < jmaxLvl + 1; j++) {
        for (int i = isw; i < imaxLvl + 1; i += 2) {

          R(i, j, k) =
              (RHS(i, j, k) -
               ((P(i + 1, j, k) - 2.0 * P(i, j, k) + P(i - 1, j, k)) * idx2 +
                (P(i, j + 1, k) - 2.0 * P(i, j, k) + P(i, j - 1, k)) * idy2 +
                (P(i, j, k + 1) - 2.0 * P(i, j, k) + P(i, j, k - 1)) * idz2));

          res += (R(i, j, k) * R(i, j, k));
        }
        isw = 3 - isw;
      }
      jsw = 3 - jsw;
    }
    ksw = 3 - ksw;
  }

  // TODO: Global reduction to sum residuals across all processes
  commReduction(&res, SUM);
  res = res / (double)(s->grid->imax * s->grid->jmax * s->grid->kmax);

  return res;
}

static bool checkDimension(int imaxLvl, int jmaxLvl, int kmaxLvl) {
  if (imaxLvl < 2 || jmaxLvl < 2 || kmaxLvl < 2) {
    return true;
  }
  return false;
}

static double multiGrid(Solver *s, double *p, double *rhs, int level,
                        int imaxLevel, int jmaxLevel, int kmaxLevel) {

  /* imax, jmax, kmax function parameters tell you
   * How many active points exist at this level
   * What the loop bounds should be
   * How to scale indices when accessing parent/child grids
   *  At level 0: imax=64, jmax=64, kmax=64
   *  At level 1: imax=32, jmax=32, kmax=32
   *  At level 2: imax=16, jmax=16, kmax=16
   * ...
   *  In Computation scaling occurs as follows:
   *  2*i, 2*j, 2*k = scaling to map between coarse and fine grid positions
   */

  double res = 0.0;

  // coarsest level
  if (level == COARSEST_LEVEL ||
      checkDimension(imaxLevel, jmaxLevel, kmaxLevel)) {
    for (int i = 0; i < s->presmooth; i++) {
      smooth(s, p, rhs, level, imaxLevel, jmaxLevel, kmaxLevel);
    }
    return calculateResidual(s, p, rhs, level, imaxLevel, jmaxLevel, kmaxLevel);
  }

  // pre-smoothing
  for (int i = 0; i < s->presmooth; i++) {
    smooth(s, p, rhs, level, imaxLevel, jmaxLevel, kmaxLevel);
    if (level == FINEST_LEVEL)
      setBoundaryCondition(s, p, imaxLevel, jmaxLevel, kmaxLevel);
  }

  commExchangeLevel(s->comm, p, level);
  res = calculateResidual(s, p, rhs, level, imaxLevel, jmaxLevel, kmaxLevel);

  // restrict
  restrictMG(s, level, imaxLevel, jmaxLevel, kmaxLevel);

  // MGSolver on residual and error.
  multiGrid(s, s->e[level + 1], s->r[level + 1], level + 1, imaxLevel / 2,
            jmaxLevel / 2, kmaxLevel / 2);

  // prolongate
  prolongate(s, level, imaxLevel, jmaxLevel, kmaxLevel);

  // correct p on finer level using residual
  correct(s, p, level, imaxLevel, jmaxLevel, kmaxLevel);
  if (level == FINEST_LEVEL)
    setBoundaryCondition(s, p, imaxLevel, jmaxLevel, kmaxLevel);

  // post-smoothing
  for (int i = 0; i < s->postsmooth; i++) {
    smooth(s, p, rhs, level, imaxLevel, jmaxLevel, kmaxLevel);
    if (level == FINEST_LEVEL)
      setBoundaryCondition(s, p, imaxLevel, jmaxLevel, kmaxLevel);
  }

  return res;
}

void initSolver(Solver *s, Discretization *d, Parameter *p) {
  s->eps = p->eps;
  s->omega = p->omg;
  s->itermax = p->itermax;
  s->levels = p->levels;
  s->grid = &d->grid;
  s->comm = &d->comm;
  s->presmooth = p->presmooth;
  s->postsmooth = p->postsmooth;

  int imaxLocal = s->comm->imaxLocal;
  int jmaxLocal = s->comm->jmaxLocal;
  int kmaxLocal = s->comm->kmaxLocal;
  int levels = s->levels;

  int minLocal = imaxLocal;
  if (jmaxLocal < minLocal)
    minLocal = jmaxLocal;
  if (kmaxLocal < minLocal)
    minLocal = kmaxLocal;

  int maxLevels = 1;
  int tmp = minLocal;
  while (tmp >= 4 && tmp % 2 == 0) {
    maxLevels++;
    tmp /= 2;
  }

  if (levels > maxLevels) {
    if (commIsMaster(s->comm)) {
      printf("WARNING: Reducing multigrid levels from %d to %d\n", levels,
             maxLevels);
      printf("  Smallest local dimension: %d (local grid: %dx%dx%d)\n",
             minLocal, imaxLocal, jmaxLocal, kmaxLocal);
      printf(
          "  For %d levels, smallest local dimension must be divisible by %d\n",
          levels, 1 << (levels - 1));
    }
    // levels = maxLevels;
    // s->levels = levels;
  }

  if (commIsMaster(s->comm)) {
    printf("Using Multigrid solver with %d levels\n", levels);
    printf("  Local grid: %dx%dx%d\n", imaxLocal, jmaxLocal, kmaxLocal);
    printf("  Coarsest local grid: %dx%dx%d\n", imaxLocal >> (levels - 1),
           jmaxLocal >> (levels - 1), kmaxLocal >> (levels - 1));
  }

  commSetupMG(s->comm, levels);

  // WARN: EVERY LEVEL IS THE SAME SIZE
  s->r = malloc(levels * sizeof(double *));
  s->e = malloc(levels * sizeof(double *));

  size_t size = (imaxLocal + 2) * (jmaxLocal + 2) * (kmaxLocal + 2);

  for (int j = 0; j < levels; j++) {
    s->r[j] = allocate(64, size * sizeof(double));
    s->e[j] = allocate(64, size * sizeof(double));

    for (size_t i = 0; i < size; i++) {
      s->r[j][i] = 0.0;
      s->e[j][i] = 0.0;
    }
  }
}

double solve(Solver *s, double *p, double *rhs) {
  double res = multiGrid(s, p, rhs, 0, s->comm->imaxLocal, s->comm->jmaxLocal,
                         s->comm->kmaxLocal);
  return res;

#ifdef VERBOSE
  if (commIsMaster(s->comm)) {
    printf("Residuum: %.6f\n", res);
  }
#endif
}
