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
#include "solver.h"
#include "util.h"

#define FINEST_LEVEL   0
#define COARSEST_LEVEL (s->levels - 1)
#define S(i, j, k)     s[(k) * (imax + 2) * (jmax + 2) + (j) * (imax + 2) + (i)]
#define E(i, j, k)     e[(k) * (imax + 2) * (jmax + 2) + (j) * (imax + 2) + (i)]
#define R(i, j, k)     r[(k) * (imax + 2) * (jmax + 2) + (j) * (imax + 2) + (i)]
#define OLD(i, j, k)   old[(k) * (imax + 2) * (jmax + 2) + (j) * (imax + 2) + (i)]

static void restrictMG(Solver* s, int level, int imax, int jmax, int kmax)
{
    double* r   = s->r[level + 1];
    double* old = s->r[level];

    for (int k = 1; k < (kmax + 1) / 2; k++) {
        for (int j = 1; j < (jmax + 1) / 2; j++) {
            for (int i = 1; i < (imax + 1) / 2; ++i) {
                R(i, j, k) = (OLD(2 * i - 1, 2 * j - 1, 2 * k) +
                                 OLD(2 * i, 2 * j - 1, 2 * k) * 2 +
                                 OLD(2 * i + 1, 2 * j - 1, 2 * k) +
                                 OLD(2 * i - 1, 2 * j, 2 * k) * 2 +
                                 OLD(2 * i, 2 * j, 2 * k) * 8 +
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

static void prolongate(Solver* s, int level, int imax, int jmax, int kmax)
{
    double* old = s->r[level + 1];
    double* e   = s->r[level];

    for (int k = 2; k < kmax + 1; k++) {
        for (int j = 2; j < jmax + 1; j++) {
            for (int i = 2; i < imax + 1; i++) {
                E(i, j, k) = OLD((i + 1) / 2, (j + 1) / 2, (k + 1) / 2);
            }
        }
    }
}

static void correct(Solver* s, double* p, int level, int imax, int jmax, int kmax)
{
    double* e = s->e[level];

    for (int k = 1; k < kmax + 1; ++k) {
        for (int j = 1; j < jmax + 1; ++j) {
            for (int i = 1; i < imax + 1; ++i) {
                P(i, j, k) += E(i, j, k);
            }
        }
    }
}

static void setBoundaryCondition(double* p, int imax, int jmax, int kmax)
{
    for (int j = 1; j < jmax + 1; j++) {
        for (int i = 1; i < imax + 1; i++) {
            P(i, j, 0)        = P(i, j, 1);
            P(i, j, kmax + 1) = P(i, j, kmax);
        }
    }

    for (int k = 1; k < kmax + 1; k++) {
        for (int i = 1; i < imax + 1; i++) {
            P(i, 0, k)        = P(i, 1, k);
            P(i, jmax + 1, k) = P(i, jmax, k);
        }
    }

    for (int k = 1; k < kmax + 1; k++) {
        for (int j = 1; j < jmax + 1; j++) {
            P(0, j, k)        = P(1, j, k);
            P(imax + 1, j, k) = P(imax, j, k);
        }
    }
}

static void smooth(
    Solver* s, double* p, double* rhs, int level, int imax, int jmax, int kmax)
{
    double eps    = s->eps;
    int itermax   = s->itermax;
    double dx2    = s->grid->dx * s->grid->dx;
    double dy2    = s->grid->dy * s->grid->dy;
    double dz2    = s->grid->dz * s->grid->dz;
    double idx2   = 1.0 / dx2;
    double idy2   = 1.0 / dy2;
    double idz2   = 1.0 / dz2;
    double factor = s->omega * 0.5 * (dx2 * dy2 * dz2) /
                    (dy2 * dz2 + dx2 * dz2 + dx2 * dy2);
    double* r    = s->r[level];
    double epssq = eps * eps;
    int it       = 0;
    int pass, ksw, jsw, isw;
    double res = 1.0;

    ksw = 1;

    for (pass = 0; pass < 2; pass++) {
        jsw = ksw;

        for (int k = 1; k < kmax + 1; k++) {
            isw = jsw;
            for (int j = 1; j < jmax + 1; j++) {
                for (int i = isw; i < imax + 1; i += 2) {

                    P(i, j, k) -=
                        factor *
                        (RHS(i, j, k) -
                            ((P(i + 1, j, k) - 2.0 * P(i, j, k) + P(i - 1, j, k)) * idx2 +
                                (P(i, j + 1, k) - 2.0 * P(i, j, k) + P(i, j - 1, k)) *
                                    idy2 +
                                (P(i, j, k + 1) - 2.0 * P(i, j, k) + P(i, j, k - 1)) *
                                    idz2));
                }
                isw = 3 - isw;
            }
            jsw = 3 - jsw;
        }
        ksw = 3 - ksw;
    }
}

static double calculateResidual(
    Solver* s, double* p, double* rhs, int level, int imax, int jmax, int kmax)
{
    double eps    = s->eps;
    int itermax   = s->itermax;
    double dx2    = s->grid->dx * s->grid->dx;
    double dy2    = s->grid->dy * s->grid->dy;
    double dz2    = s->grid->dz * s->grid->dz;
    double idx2   = 1.0 / dx2;
    double idy2   = 1.0 / dy2;
    double idz2   = 1.0 / dz2;
    double factor = s->omega * 0.5 * (dx2 * dy2 * dz2) /
                    (dy2 * dz2 + dx2 * dz2 + dx2 * dy2);
    double* r    = s->r[level];
    double epssq = eps * eps;
    int it       = 0;
    int pass, ksw, jsw, isw;
    double res = 1.0;

    ksw = 1;

    for (pass = 0; pass < 2; pass++) {
        jsw = ksw;

        for (int k = 1; k < kmax + 1; k++) {
            isw = jsw;
            for (int j = 1; j < jmax + 1; j++) {
                for (int i = isw; i < imax + 1; i += 2) {

                    R(i,
                        j,
                        k) = (RHS(i, j, k) -
                              ((P(i + 1, j, k) - 2.0 * P(i, j, k) + P(i - 1, j, k)) *
                                      idx2 +
                                  (P(i, j + 1, k) - 2.0 * P(i, j, k) + P(i, j - 1, k)) *
                                      idy2 +
                                  (P(i, j, k + 1) - 2.0 * P(i, j, k) + P(i, j, k - 1)) *
                                      idz2));

                    res += (R(i, j, k) * R(i, j, k));
                }
                isw = 3 - isw;
            }
            jsw = 3 - jsw;
        }
        ksw = 3 - ksw;
    }

    res = res / (double)(imax * jmax * kmax);

    return res;
}

static bool checkDimension(int imax, int jmax, int kmax)
{
    if (imax < 2 || jmax < 2 || kmax < 2) {
        return true;
    }
    return false;
}

static double multiGrid(
    Solver* s, double* p, double* rhs, int level, int imax, int jmax, int kmax)
{
    double res = 0.0;

    // coarsest level
    if (level == COARSEST_LEVEL || checkDimension(imax, jmax, kmax)) {
        for (int i = 0; i < s->presmooth; i++) {
            smooth(s, p, rhs, level, imax, jmax, kmax);
        }
        return calculateResidual(s, p, rhs, level, imax, jmax, kmax);
    }

    // pre-smoothing
    for (int i = 0; i < s->presmooth; i++) {
        smooth(s, p, rhs, level, imax, jmax, kmax);
        if (level == FINEST_LEVEL) setBoundaryCondition(p, imax, jmax, kmax);
    }

    res = calculateResidual(s, p, rhs, level, imax, jmax, kmax);

    // restrict
    restrictMG(s, level, imax, jmax, kmax);

    // MGSolver on residual and error.
    multiGrid(s,
        s->e[level + 1],
        s->r[level + 1],
        level + 1,
        imax / 2,
        jmax / 2,
        kmax / 2);

    // prolongate
    prolongate(s, level, imax, jmax, kmax);

    // correct p on finer level using residual
    correct(s, p, level, imax, jmax, kmax);
    if (level == FINEST_LEVEL) setBoundaryCondition(p, imax, jmax, kmax);

    // post-smoothing
    for (int i = 0; i < s->postsmooth; i++) {
        smooth(s, p, rhs, level, imax, jmax, kmax);
        if (level == FINEST_LEVEL) setBoundaryCondition(p, imax, jmax, kmax);
    }

    return res;
}

void initSolver(Solver* s, Discretization* d, Parameter* p)
{
    s->eps        = p->eps;
    s->omega      = p->omg;
    s->itermax    = p->itermax;
    s->levels     = p->levels;
    s->grid       = &d->grid;
    s->presmooth  = p->presmooth;
    s->postsmooth = p->postsmooth;

    int imax   = s->grid->imax;
    int jmax   = s->grid->jmax;
    int kmax   = s->grid->kmax;
    int levels = s->levels;
    printf("Using Multigrid solver with %d levels\n", levels);

    s->r = malloc(levels * sizeof(double*));
    s->e = malloc(levels * sizeof(double*));

    size_t size = (imax + 2) * (jmax + 2) * (kmax + 2);

    for (int j = 0; j < levels; j++) {
        s->r[j] = allocate(64, size * sizeof(double));
        s->e[j] = allocate(64, size * sizeof(double));

        for (size_t i = 0; i < size; i++) {
            s->r[j][i] = 0.0;
            s->e[j][i] = 0.0;
        }
    }
}

double solve(Solver* s, double* p, double* rhs)
{
    double res = multiGrid(s, p, rhs, 0, s->grid->imax, s->grid->jmax, s->grid->kmax);

#ifdef VERBOSE
    printf("Residuum: %.6f\n", res);
#endif

    return res;
}
