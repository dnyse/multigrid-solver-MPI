/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of nusif-solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vtkWriter.h"
#define G(v, i, j, k) v[(k)*imax * jmax + (j)*imax + (i)]

static float floatSwap(float f)
{
    union {
        float f;
        char b[4];
    } dat1, dat2;

    dat1.f    = f;
    dat2.b[0] = dat1.b[3];
    dat2.b[1] = dat1.b[2];
    dat2.b[2] = dat1.b[1];
    dat2.b[3] = dat1.b[0];
    return dat2.f;
}

static void writeHeader(VtkOptions* o)
{
    fprintf(o->fh, "# vtk DataFile Version 3.0\n");
    fprintf(o->fh, "PAMPI cfd solver output\n");
    if (o->fmt == ASCII) {
        fprintf(o->fh, "ASCII\n");
    } else if (o->fmt == BINARY) {
        fprintf(o->fh, "BINARY\n");
    }

    fprintf(o->fh, "DATASET STRUCTURED_POINTS\n");
    fprintf(o->fh, "DIMENSIONS %d %d %d\n", o->grid.imax, o->grid.jmax, o->grid.kmax);
    fprintf(o->fh,
        "ORIGIN %f %f %f\n",
        o->grid.dx * 0.5,
        o->grid.dy * 0.5,
        o->grid.dz * 0.5);
    fprintf(o->fh, "SPACING %f %f %f\n", o->grid.dx, o->grid.dy, o->grid.dz);
    fprintf(o->fh, "POINT_DATA %d\n", o->grid.imax * o->grid.jmax * o->grid.kmax);
}

void vtkOpen(VtkOptions* o, char* problem)
{
    char filename[50];
    snprintf(filename, 50, "%s.vtk", problem);
    o->fh = fopen(filename, "w");
    writeHeader(o);

    printf("Writing VTK output for %s\n", problem);
}

void vtkScalar(VtkOptions* o, char* name, double* s)
{
    int imax = o->grid.imax;
    int jmax = o->grid.jmax;
    int kmax = o->grid.kmax;

    printf("Register scalar %s\n", name);

    if (o->fh == NULL) {
        printf("vtkWriter not initialize! Call vtkOpen first!\n");
        exit(EXIT_FAILURE);
    }
    fprintf(o->fh, "SCALARS %s float\n", name);
    fprintf(o->fh, "LOOKUP_TABLE default\n");

    for (int k = 0; k < kmax; k++) {
        for (int j = 0; j < jmax; j++) {
            for (int i = 0; i < imax; i++) {
                if (o->fmt == ASCII) {
                    fprintf(o->fh, "%f\n", G(s, i, j, k));
                } else if (o->fmt == BINARY) {
                    fwrite((float[1]) { floatSwap(G(s, i, j, k)) },
                        sizeof(float),
                        1,
                        o->fh);
                }
            }
        }
    }
    if (o->fmt == BINARY) fprintf(o->fh, "\n");
}

void vtkVector(VtkOptions* o, char* name, VtkVector vec)
{
    int imax = o->grid.imax;
    int jmax = o->grid.jmax;
    int kmax = o->grid.kmax;

    if (o->fh == NULL) {
        printf("vtkWriter not initialize! Call vtkOpen first!\n");
        exit(EXIT_FAILURE);
    }

    fprintf(o->fh, "VECTORS %s float\n", name);

    for (int k = 0; k < kmax; k++) {
        for (int j = 0; j < jmax; j++) {
            for (int i = 0; i < imax; i++) {
                if (o->fmt == ASCII) {
                    fprintf(o->fh,
                        "%f %f %f\n",
                        G(vec.u, i, j, k),
                        G(vec.v, i, j, k),
                        G(vec.w, i, j, k));
                } else if (o->fmt == BINARY) {
                    fwrite((float[3]) { floatSwap(G(vec.u, i, j, k)),
                               floatSwap(G(vec.v, i, j, k)),
                               floatSwap(G(vec.w, i, j, k)) },
                        sizeof(float),
                        3,
                        o->fh);
                }
            }
        }
    }
    if (o->fmt == BINARY) fprintf(o->fh, "\n");
}

void vtkClose(VtkOptions* o)
{
    fclose(o->fh);
    o->fh = NULL;
}
