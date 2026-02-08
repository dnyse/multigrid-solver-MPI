/*
 * Copyright (C) 2024 NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file.
 */
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "allocate.h"
#include "comm.h"
#include "discretization.h"
#include "parameter.h"
#include "progress.h"
#include "solver.h"
#include "timing.h"
#include "vtkWriter.h"
#include <math.h>

int main(int argc, char **argv) {
  double timeStart, timeStop;
  Parameter p;
  Discretization d;
  Solver s;
  initParameter(&p);
  commInit(&d.comm, argc, argv);

  if (argc != 2) {
    printf("Usage: %s <configFile>\n", argv[0]);
    exit(EXIT_SUCCESS);
  }

  readParameter(&p, argv[1]);
  commPartition(&d.comm, p.kmax, p.jmax, p.imax);
  if (commIsMaster(&d.comm)) {
    printParameter(&p);
  }
  initDiscretization(&d, &p);
  initSolver(&s, &d, &p);
#ifndef VERBOSE
  if (commIsMaster(&d.comm)) {
    initProgress(d.te);
  }
#endif

  double tau = d.tau;
  double te = d.te;
  double t = 0.0;
  int nt = 0;
  double conv = DBL_MAX;
  double eps2 = p.eps * p.eps;

  timeStart = getTimeStamp();
  while (t <= te && conv > eps2) {
    if (tau > 0.0)
      computeTimestep(&d);
    setBoundaryConditions(&d);
    setSpecialBoundaryCondition(&d);
    computeFG(&d);
    computeRHS(&d);
    // if (nt % 10 == 0) {
      // normalizePressure(&d);
      // commExchange(&d.comm, d.p);
    // }

    conv = solve(&s, d.p, d.rhs);
    adaptUV(&d);
    t += d.dt;
    nt++;

#ifdef VERBOSE
    if (commIsMaster(&d.comm)) {
      printf("TIME %f, TIMESTEP %f, Residual %f \n", t, d.dt, sqrt(conv));
    }
#else
    if (commIsMaster(&d.comm)) {
      printProgress(t);
    }
#endif
  }
  timeStop = getTimeStamp();
#ifndef VERBOSE
  if (commIsMaster(&d.comm)) {
    stopProgress();
  }
#endif
  if (commIsMaster(&d.comm)) {
    printf("Solution took %.2fs in %d iterations.\n", timeStop - timeStart, nt);
    printf("Final residiuum %.2fs.\n", sqrt(conv));
  }

  double *pg, *ug, *vg, *wg;

  if (commIsMaster(&d.comm)) {
    size_t bytesize = d.grid.imax * d.grid.jmax * d.grid.kmax * sizeof(double);

    pg = allocate(64, bytesize);
    ug = allocate(64, bytesize);
    vg = allocate(64, bytesize);
    wg = allocate(64, bytesize);
  }

  commCollectResult(&d.comm, ug, vg, wg, pg, d.u, d.v, d.w, d.p, d.grid.kmax,
                    d.grid.jmax, d.grid.imax);

  if (commIsMaster(&d.comm)) {
    VtkOptions opts = {.grid = d.grid};
    vtkOpen(&opts, d.problem);
    vtkScalar(&opts, "pressure", pg);
    vtkVector(&opts, "velocity", (VtkVector){ug, vg, wg});
    vtkClose(&opts);
  }

  commFinalize(&d.comm);
  return EXIT_SUCCESS;
}
