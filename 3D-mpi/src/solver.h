/*
 * Copyright (C) 2022 NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of nusif-solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 */
#ifndef __SOLVER_H_
#define __SOLVER_H_
#include "comm.h"
#include "discretization.h"
#include "grid.h"
#include "parameter.h"

typedef struct {
  /* geometry and grid information */
  Grid *grid;
  /* parameters */
  double eps, omega, rho;
  /* time stepping */
  int itermax;
  int levels;
  /* multigrid paramters */
  double **r, **e;
  int presmooth, postsmooth;

  Comm *comm;

} Solver;

extern void initSolver(Solver *, Discretization *, Parameter *);
extern void solve(Solver *, double *, double *);
#endif
