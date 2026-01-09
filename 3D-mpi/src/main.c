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
  initProgress(d.te);
#endif

  double tau = d.tau;
  double te = d.te;
  double t = 0.0;
  int nt = 0;

  timeStart = getTimeStamp();
  while (t <= te) {
    if (tau > 0.0)
      computeTimestep(&d);
    setBoundaryConditions(&d);
    setSpecialBoundaryCondition(&d);
    computeFG(&d);
    computeRHS(&d);
    // TODO: Check Correctness -> MPI
    if (nt % 100 == 0)
      normalizePressure(&d);

    solve(&s, d.p, d.rhs);
    adaptUV(&d);
    t += d.dt;
    nt++;

#ifdef VERBOSE
    if (commIsMaster(&s.comm)) {
      printf("TIME %f , TIMESTEP %f\n", t, s.dt);
    }
#else
    printProgress(t);
#endif
  }
  timeStop = getTimeStamp();
#ifndef VERBOSE
  stopProgress();
#endif
  if (commIsMaster(&d.comm)) {
    printf("Solution took %.2fs\n", timeStop - timeStart);
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
