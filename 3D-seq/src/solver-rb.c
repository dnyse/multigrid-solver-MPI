/*
 * Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of nusif-solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "solver.h"
#include "timing.h"
#include "util.h"

void initSolver(Solver* s, Discretization* d, Parameter* p)
{
    s->grid    = &d->grid;
    s->itermax = p->itermax;
    s->eps     = p->eps;
    s->omega   = p->omg;
}

double solve(Solver* s, double* p, double* rhs)
{
    int imax      = s->grid->imax;
    int jmax      = s->grid->jmax;
    int kmax      = s->grid->kmax;
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
    double epssq = eps * eps;
    int it       = 0;
    double res   = 1.0;
    int pass, ksw, jsw, isw;

    double timeStart = getTimeStamp();
    while ((res >= epssq) && (it < itermax)) {
        res = 0.0;
        ksw = 1;

        for (pass = 0; pass < 2; pass++) {
            jsw = ksw;

            for (int k = 1; k < kmax + 1; k++) {
                isw = jsw;
                for (int j = 1; j < jmax + 1; j++) {
                    for (int i = isw; i < imax + 1; i += 2) {

                        double r =
                            RHS(i, j, k) -
                            ((P(i + 1, j, k) - 2.0 * P(i, j, k) + P(i - 1, j, k)) * idx2 +
                                (P(i, j + 1, k) - 2.0 * P(i, j, k) + P(i, j - 1, k)) *
                                    idy2 +
                                (P(i, j, k + 1) - 2.0 * P(i, j, k) + P(i, j, k - 1)) *
                                    idz2);

                        P(i, j, k) -= (factor * r);
                        res += (r * r);
                    }
                    isw = 3 - isw;
                }
                jsw = 3 - jsw;
            }
            ksw = 3 - ksw;
        }

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

        res = res / (double)(imax * jmax * kmax);
#ifdef DEBUG
        printf("%d Residuum: %e\n", it, res);
#endif
        it++;
    }
    double timeStop = getTimeStamp();

#ifdef VERBOSE
    printf("Solver took %d iterations to reach %f\n", it, sqrt(res));
    double datapoints = (kmax + 2) * (jmax + 2) * (imax + 2);
    double points     = kmax * jmax * imax;
    double mlups      = (it * points) * 1.E-6;
    printf("Dataset size: %f MB Performance: %f MLUPS/s\n",
        datapoints * 16.0 * 1.E-6,
        mlups / (timeStop - timeStart));
#endif

    return res;
}
