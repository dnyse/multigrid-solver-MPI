/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of nusif-solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 */
#ifndef __SOLVER_H_
#define __SOLVER_H_
#include "discretization.h"
#include "grid.h"
#include "parameter.h"

typedef struct {
    /* geometry and grid information */
    Grid* grid;
    /* parameters */
    double eps, omega, rho;
    int itermax;
    int levels;
    double **r, **e;
    int presmooth, postsmooth;
} Solver;

extern void initSolver(Solver*, Discretization*, Parameter*);
extern double solve(Solver*, double*, double*);

#endif
