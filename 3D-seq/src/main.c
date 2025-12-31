/*
 * Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "allocate.h"
#include "discretization.h"
#include "parameter.h"
#include "progress.h"
#include "solver.h"
#include "timing.h"
#include "vtkWriter.h"

#define G(v, i, j, k) v[(k) * (imax + 2) * (jmax + 2) + (j) * (imax + 2) + (i)]

static FILE* initResidualWriter()
{
    FILE* fp;
    fp = fopen("residual.dat", "w");

    if (fp == NULL) {
        printf("Error!\n");
        exit(EXIT_FAILURE);
    }

    return fp;

}

static void writeResidual(FILE* fp, double ts, double res)
{
    fprintf(fp, "%f, %f\n", ts, res);
}

static void createBulkArrays(
    Discretization* s, double* pg, double* ug, double* vg, double* wg)
{
    int imax = s->grid.imax;
    int jmax = s->grid.jmax;
    int kmax = s->grid.kmax;
    int idx  = 0;

    for (int k = 1; k < kmax + 1; k++) {
        for (int j = 1; j < jmax + 1; j++) {
            for (int i = 1; i < imax + 1; i++) {
                pg[idx++] = G(s->p, i, j, k);
            }
        }
    }

    idx = 0;

    for (int k = 1; k < kmax + 1; k++) {
        for (int j = 1; j < jmax + 1; j++) {
            for (int i = 1; i < imax + 1; i++) {
                ug[idx++] = (G(s->u, i, j, k) + G(s->u, i - 1, j, k)) / 2.0;
            }
        }
    }

    idx = 0;

    for (int k = 1; k < kmax + 1; k++) {
        for (int j = 1; j < jmax + 1; j++) {
            for (int i = 1; i < imax + 1; i++) {
                vg[idx++] = (G(s->v, i, j, k) + G(s->v, i, j - 1, k)) / 2.0;
            }
        }
    }

    idx = 0;

    for (int k = 1; k < kmax + 1; k++) {
        for (int j = 1; j < jmax + 1; j++) {
            for (int i = 1; i < imax + 1; i++) {
                wg[idx++] = (G(s->w, i, j, k) + G(s->w, i, j, k - 1)) / 2.0;
            }
        }
    }
}

int main(int argc, char** argv)
{
    double timeStart, timeStop;
    Parameter p;
    Discretization d;
    Solver s;
    initParameter(&p);

    FILE* fp;
    fp = initResidualWriter();

    if (argc != 2) {
        printf("Usage: %s <configFile>\n", argv[0]);
        exit(EXIT_SUCCESS);
    }

    readParameter(&p, argv[1]);
    printParameter(&p);
    initDiscretization(&d, &p);
    initSolver(&s, &d, &p);
#ifndef VERBOSE
    initProgress(d.te);
#endif

    double tau = d.tau;
    double te  = d.te;
    double t   = 0.0;
    int nt     = 0;
    double res = 0.0;

    timeStart = getTimeStamp();
    while (t <= te) {
        if (tau > 0.0) computeTimestep(&d);
        setBoundaryConditions(&d);
        setSpecialBoundaryCondition(&d);
        computeFG(&d);
        computeRHS(&d);
        if (nt % 100 == 0) normalizePressure(&d);
        res = solve(&s, d.p, d.rhs);
        adaptUV(&d);

        writeResidual(fp, t, res);

        t += d.dt;
        nt++;

#ifdef VERBOSE
        printf("TIME %f , TIMESTEP %f\n", t, d.dt);
#else
        printProgress(t);
#endif
    }
    timeStop = getTimeStamp();
    stopProgress();
    printf("Solution took %.2fs\n", timeStop - timeStart);

    timeStart = getTimeStamp();
    double *pg, *ug, *vg, *wg;

    size_t bytesize = (size_t)(d.grid.imax * d.grid.jmax * d.grid.kmax) * sizeof(double);

    pg = allocate(64, bytesize);
    ug = allocate(64, bytesize);
    vg = allocate(64, bytesize);
    wg = allocate(64, bytesize);
    
    fclose(fp);
    createBulkArrays(&d, pg, ug, vg, wg);
    VtkOptions opts = { .grid = d.grid };
    vtkOpen(&opts, d.problem);
    vtkScalar(&opts, "pressure", pg);
    vtkVector(&opts, "velocity", (VtkVector) { ug, vg, wg });
    vtkClose(&opts);

    timeStop = getTimeStamp();
    printf("Result output took %.2fs\n", timeStop - timeStart);

    return EXIT_SUCCESS;
}
