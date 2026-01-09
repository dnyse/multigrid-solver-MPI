#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "comm.h"
#include "parameter.h"
#include "util.h"
#include "discretization.h"
#include "allocate.h"

#define P(i, j, k)                                                                       \
    p[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]
#define F(i, j, k)                                                                       \
    f[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]
#define G(i, j, k)                                                                       \
    g[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]
#define H(i, j, k)                                                                       \
    h[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]
#define U(i, j, k)                                                                       \
    u[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]
#define V(i, j, k)                                                                       \
    v[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]
#define W(i, j, k)                                                                       \
    w[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]
#define RHS(i, j, k)                                                                     \
    rhs[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]


static void printConfig(Discretization* d)
{
    if (commIsMaster(&d->comm)) {
        printf("Parameters for #%s#\n", d->problem);
        printf("BC Left:%d Right:%d Bottom:%d Top:%d Front:%d Back:%d\n",
            d->bcLeft,
            d->bcRight,
            d->bcBottom,
            d->bcTop,
            d->bcFront,
            d->bcBack);
        printf("\tReynolds number: %.2f\n", d->re);
        printf("\tGx Gy: %.2f %.2f %.2f\n", d->gx, d->gy, d->gz);
        printf("Geometry data:\n");
        printf("\tDomain box size (x, y, z): %.2f, %.2f, %.2f\n",
            d->grid.xlength,
            d->grid.ylength,
            d->grid.zlength);
        printf("\tCells (x, y, z): %d, %d, %d\n",
            d->grid.imax,
            d->grid.jmax,
            d->grid.kmax);
        printf("\tCell size (dx, dy, dz): %f, %f, %f\n",
            d->grid.dx,
            d->grid.dy,
            d->grid.dz);
        printf("Timestep parameters:\n");
        printf("\tDefault stepsize: %.2f, Final time %.2f\n", d->dt, d->te);
        printf("\tdt bound: %.6f\n", d->dtBound);
        printf("\tTau factor: %.2f\n", d->tau);
        printf("Iterative parameters:\n");
        // printf("\tMax iterations: %d\n", d->itermax);
        printf("\tepsilon (stopping tolerance) : %f\n", d->eps);
        printf("\tgamma factor: %f\n", d->gamma);
        printf("\tomega (SOR relaxation): %f\n", d->omega);
    }
    commPrintConfig(&d->comm);
}

void initDiscretization(Discretization* d, Parameter* p){
d->problem = p->name;
  d->bcLeft = p->bcLeft;
  d->bcRight = p->bcRight;
  d->bcBottom = p->bcBottom;
  d->bcTop = p->bcTop;
  d->bcFront = p->bcFront;
  d->bcBack = p->bcBack;

  d->grid.imax = p->imax;
  d->grid.jmax = p->jmax;
  d->grid.kmax = p->kmax;
  d->grid.xlength = p->xlength;
  d->grid.ylength = p->ylength;
  d->grid.zlength = p->zlength;
  d->grid.dx = p->xlength / p->imax;
  d->grid.dy = p->ylength / p->jmax;
  d->grid.dz = p->zlength / p->kmax;
 d->re = p->re;
  d->gx = p->gx;
  d->gy = p->gy;
  d->gz = p->gz;
  d->dt = p->dt;
  d->te = p->te;
  d->tau = p->tau;
  d->gamma = p->gamma;

  /* allocate arrays */
  int imaxLocal = d->comm.imaxLocal;
  int jmaxLocal = d->comm.jmaxLocal;
  int kmaxLocal = d->comm.kmaxLocal;
  size_t size = (imaxLocal + 2) * (jmaxLocal + 2) * (kmaxLocal + 2);

  d->u = allocate(64, size * sizeof(double));
  d->v = allocate(64, size * sizeof(double));
  d->w = allocate(64, size * sizeof(double));
  d->p = allocate(64, size * sizeof(double));
  d->rhs = allocate(64, size * sizeof(double));
  d->f = allocate(64, size * sizeof(double));
  d->g = allocate(64, size * sizeof(double));
  d->h = allocate(64, size * sizeof(double));

  for (int i = 0; i < size; i++) {
    d->u[i] = p->u_init;
    d->v[i] = p->v_init;
    d->w[i] = p->w_init;
    d->p[i] = p->p_init;
    d->rhs[i] = 0.0;
    d->f[i] = 0.0;
    d->g[i] = 0.0;
    d->h[i] = 0.0;
  }

  double dx = d->grid.dx;
  double dy = d->grid.dy;
  double dz = d->grid.dz;

  double invSqrSum = 1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz);
  d->dtBound = 0.5 * d->re * 1.0 / invSqrSum;

}

void computeRHS(Discretization* d)
{
    int imaxLocal = d->comm.imaxLocal;
    int jmaxLocal = d->comm.jmaxLocal;
    int kmaxLocal = d->comm.kmaxLocal;

    double idx = 1.0 / d->grid.dx;
    double idy = 1.0 / d->grid.dy;
    double idz = 1.0 / d->grid.dz;
    double idt = 1.0 / d->dt;

    double* rhs = d->rhs;
    double* f   = d->f;
    double* g   = d->g;
    double* h   = d->h;

    commShift(&d->comm, f, g, h);

    for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
            for (int i = 1; i < imaxLocal + 1; i++) {
                RHS(i, j, k) = ((F(i, j, k) - F(i - 1, j, k)) * idx +
                                   (G(i, j, k) - G(i, j - 1, k)) * idy +
                                   (H(i, j, k) - H(i, j, k - 1)) * idz) *
                               idt;
            }
        }
    }
}

static double maxElement(Discretization* d, double* m)
{
    int size = (d->comm.imaxLocal + 2) * (d->comm.jmaxLocal + 2) *
               (d->comm.kmaxLocal + 2);
    double maxval = DBL_MIN;

    for (int i = 0; i < size; i++) {
        maxval = MAX(maxval, fabs(m[i]));
    }
    commReduction(&maxval, MAX);
    return maxval;
}

void normalizePressure(Discretization* d)
{
    int imaxLocal = d->comm.imaxLocal;
    int jmaxLocal = d->comm.jmaxLocal;
    int kmaxLocal = d->comm.kmaxLocal;

    double* p   = d->p;
    double avgP = 0.0;

    for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
            for (int i = 1; i < imaxLocal + 1; i++) {
                avgP += P(i, j, k);
            }
        }
    }
    commReduction(&avgP, SUM);
    avgP /= (d->grid.imax * d->grid.jmax * d->grid.kmax);

    for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
            for (int i = 1; i < imaxLocal + 1; i++) {
                P(i, j, k) = P(i, j, k) - avgP;
            }
        }
    }
}

void computeTimestep(Discretization* d)
{
    double dt = d->dtBound;
    double dx = d->grid.dx;
    double dy = d->grid.dy;
    double dz = d->grid.dz;

    double umax = maxElement(d, d->u);
    double vmax = maxElement(d, d->v);
    double wmax = maxElement(d, d->w);

    if (umax > 0) {
        dt = (dt > dx / umax) ? dx / umax : dt;
    }
    if (vmax > 0) {
        dt = (dt > dy / vmax) ? dy / vmax : dt;
    }
    if (wmax > 0) {
        dt = (dt > dz / wmax) ? dz / wmax : dt;
    }

    d->dt = dt * d->tau;
}

void setBoundaryConditions(Discretization* d)
{
    int imaxLocal = d->comm.imaxLocal;
    int jmaxLocal = d->comm.jmaxLocal;
    int kmaxLocal = d->comm.kmaxLocal;

    double* u = d->u;
    double* v = d->v;
    double* w = d->w;

    if (commIsBoundary(&d->comm, TOP)) {
        switch (d->bcTop) {
        case NOSLIP:
            for (int k = 1; k < kmaxLocal + 1; k++) {
                for (int i = 1; i < imaxLocal + 1; i++) {
                    U(i, jmaxLocal + 1, k) = -U(i, jmaxLocal, k);
                    V(i, jmaxLocal, k)     = 0.0;
                    W(i, jmaxLocal + 1, k) = -W(i, jmaxLocal, k);
                }
            }
            break;
        case SLIP:
            for (int k = 1; k < kmaxLocal + 1; k++) {
                for (int i = 1; i < imaxLocal + 1; i++) {
                    U(i, jmaxLocal + 1, k) = U(i, jmaxLocal, k);
                    V(i, jmaxLocal, k)     = 0.0;
                    W(i, jmaxLocal + 1, k) = W(i, jmaxLocal, k);
                }
            }
            break;
        case OUTFLOW:
            for (int k = 1; k < kmaxLocal + 1; k++) {
                for (int i = 1; i < imaxLocal + 1; i++) {
                    U(i, jmaxLocal + 1, k) = U(i, jmaxLocal, k);
                    V(i, jmaxLocal, k)     = V(i, jmaxLocal - 1, k);
                    W(i, jmaxLocal + 1, k) = W(i, jmaxLocal, k);
                }
            }
            break;
        case PERIODIC:
            break;
        }
    }

    if (commIsBoundary(&d->comm, BOTTOM)) {
        switch (d->bcBottom) {
        case NOSLIP:
            for (int k = 1; k < kmaxLocal + 1; k++) {
                for (int i = 1; i < imaxLocal + 1; i++) {
                    U(i, 0, k) = -U(i, 1, k);
                    V(i, 0, k) = 0.0;
                    W(i, 0, k) = -W(i, 1, k);
                }
            }
            break;
        case SLIP:
            for (int k = 1; k < kmaxLocal + 1; k++) {
                for (int i = 1; i < imaxLocal + 1; i++) {
                    U(i, 0, k) = U(i, 1, k);
                    V(i, 0, k) = 0.0;
                    W(i, 0, k) = W(i, 1, k);
                }
            }
            break;
        case OUTFLOW:
            for (int k = 1; k < kmaxLocal + 1; k++) {
                for (int i = 1; i < imaxLocal + 1; i++) {
                    U(i, 0, k) = U(i, 1, k);
                    V(i, 0, k) = V(i, 1, k);
                    W(i, 0, k) = W(i, 1, k);
                }
            }
            break;
        case PERIODIC:
            break;
        }
    }

    if (commIsBoundary(&d->comm, LEFT)) {
        switch (d->bcLeft) {
        case NOSLIP:
            for (int k = 1; k < kmaxLocal + 1; k++) {
                for (int j = 1; j < jmaxLocal + 1; j++) {
                    U(0, j, k) = 0.0;
                    V(0, j, k) = -V(1, j, k);
                    W(0, j, k) = -W(1, j, k);
                }
            }
            break;
        case SLIP:
            for (int k = 1; k < kmaxLocal + 1; k++) {
                for (int j = 1; j < jmaxLocal + 1; j++) {
                    U(0, j, k) = 0.0;
                    V(0, j, k) = V(1, j, k);
                    W(0, j, k) = W(1, j, k);
                }
            }
            break;
        case OUTFLOW:
            for (int k = 1; k < kmaxLocal + 1; k++) {
                for (int j = 1; j < jmaxLocal + 1; j++) {
                    U(0, j, k) = U(1, j, k);
                    V(0, j, k) = V(1, j, k);
                    W(0, j, k) = W(1, j, k);
                }
            }
            break;
        case PERIODIC:
            break;
        }
    }

    if (commIsBoundary(&d->comm, RIGHT)) {
        switch (d->bcRight) {
        case NOSLIP:
            for (int k = 1; k < kmaxLocal + 1; k++) {
                for (int j = 1; j < jmaxLocal + 1; j++) {
                    U(imaxLocal, j, k)     = 0.0;
                    V(imaxLocal + 1, j, k) = -V(imaxLocal, j, k);
                    W(imaxLocal + 1, j, k) = -W(imaxLocal, j, k);
                }
            }
            break;
        case SLIP:
            for (int k = 1; k < kmaxLocal + 1; k++) {
                for (int j = 1; j < jmaxLocal + 1; j++) {
                    U(imaxLocal, j, k)     = 0.0;
                    V(imaxLocal + 1, j, k) = V(imaxLocal, j, k);
                    W(imaxLocal + 1, j, k) = W(imaxLocal, j, k);
                }
            }
            break;
        case OUTFLOW:
            for (int k = 1; k < kmaxLocal + 1; k++) {
                for (int j = 1; j < jmaxLocal + 1; j++) {
                    U(imaxLocal, j, k)     = U(imaxLocal - 1, j, k);
                    V(imaxLocal + 1, j, k) = V(imaxLocal, j, k);
                    W(imaxLocal + 1, j, k) = W(imaxLocal, j, k);
                }
            }
            break;
        case PERIODIC:
            break;
        }
    }

    if (commIsBoundary(&d->comm, FRONT)) {
        switch (d->bcFront) {
        case NOSLIP:
            for (int j = 1; j < jmaxLocal + 1; j++) {
                for (int i = 1; i < imaxLocal + 1; i++) {
                    U(i, j, 0) = -U(i, j, 1);
                    V(i, j, 0) = -V(i, j, 1);
                    W(i, j, 0) = 0.0;
                }
            }
            break;
        case SLIP:
            for (int j = 1; j < jmaxLocal + 1; j++) {
                for (int i = 1; i < imaxLocal + 1; i++) {
                    U(i, j, 0) = U(i, j, 1);
                    V(i, j, 0) = V(i, j, 1);
                    W(i, j, 0) = 0.0;
                }
            }
            break;
        case OUTFLOW:
            for (int j = 1; j < jmaxLocal + 1; j++) {
                for (int i = 1; i < imaxLocal + 1; i++) {
                    U(i, j, 0) = U(i, j, 1);
                    V(i, j, 0) = V(i, j, 1);
                    W(i, j, 0) = W(i, j, 1);
                }
            }
            break;
        case PERIODIC:
            break;
        }
    }

    if (commIsBoundary(&d->comm, BACK)) {
        switch (d->bcBack) {
        case NOSLIP:
            for (int j = 1; j < jmaxLocal + 1; j++) {
                for (int i = 1; i < imaxLocal + 1; i++) {
                    U(i, j, kmaxLocal + 1) = -U(i, j, kmaxLocal);
                    V(i, j, kmaxLocal + 1) = -V(i, j, kmaxLocal);
                    W(i, j, kmaxLocal)     = 0.0;
                }
            }
            break;
        case SLIP:
            for (int j = 1; j < jmaxLocal + 1; j++) {
                for (int i = 1; i < imaxLocal + 1; i++) {
                    U(i, j, kmaxLocal + 1) = U(i, j, kmaxLocal);
                    V(i, j, kmaxLocal + 1) = V(i, j, kmaxLocal);
                    W(i, j, kmaxLocal)     = 0.0;
                }
            }
            break;
        case OUTFLOW:
            for (int j = 1; j < jmaxLocal + 1; j++) {
                for (int i = 1; i < imaxLocal + 1; i++) {
                    U(i, j, kmaxLocal + 1) = U(i, j, kmaxLocal);
                    V(i, j, kmaxLocal + 1) = V(i, j, kmaxLocal);
                    W(i, j, kmaxLocal)     = W(i, j, kmaxLocal - 1);
                }
            }
            break;
        case PERIODIC:
            break;
        }
    }
}

void setSpecialBoundaryCondition(Discretization* d)
{
    int imaxLocal = d->comm.imaxLocal;
    int jmaxLocal = d->comm.jmaxLocal;
    int kmaxLocal = d->comm.kmaxLocal;

    double* u = d->u;

    if (strcmp(d->problem, "dcavity") == 0) {
        if (commIsBoundary(&d->comm, TOP)) {
            for (int k = 1; k < kmaxLocal; k++) {
                for (int i = 1; i < imaxLocal; i++) {
                    U(i, jmaxLocal + 1, k) = 2.0 - U(i, jmaxLocal, k);
                }
            }
        }
    } else if (strcmp(d->problem, "canal") == 0) {
        if (commIsBoundary(&d->comm, LEFT)) {
            for (int k = 1; k < kmaxLocal + 1; k++) {
                for (int j = 1; j < jmaxLocal + 1; j++) {
                    U(0, j, k) = 2.0;
                }
            }
        }
    }
}

void computeFG(Discretization* d)
{
    int imaxLocal = d->comm.imaxLocal;
    int jmaxLocal = d->comm.jmaxLocal;
    int kmaxLocal = d->comm.kmaxLocal;

    double* u = d->u;
    double* v = d->v;
    double* w = d->w;
    double* f = d->f;
    double* g = d->g;
    double* h = d->h;

    double gx = d->gx;
    double gy = d->gy;
    double gz = d->gz;
    double dt = d->dt;

    double gamma     = d->gamma;
    double inverseRe = 1.0 / d->re;
    double inverseDx = 1.0 / d->grid.dx;
    double inverseDy = 1.0 / d->grid.dy;
    double inverseDz = 1.0 / d->grid.dz;
    double du2dx, dv2dy, dw2dz;
    double duvdx, duwdx, duvdy, dvwdy, duwdz, dvwdz;
    double du2dx2, du2dy2, du2dz2;
    double dv2dx2, dv2dy2, dv2dz2;
    double dw2dx2, dw2dy2, dw2dz2;

    commExchange(&d->comm, u);
    commExchange(&d->comm, v);
    commExchange(&d->comm, w);

    for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
            for (int i = 1; i < imaxLocal + 1; i++) {
                du2dx = inverseDx * 0.25 *
                            ((U(i, j, k) + U(i + 1, j, k)) *
                                    (U(i, j, k) + U(i + 1, j, k)) -
                                (U(i, j, k) + U(i - 1, j, k)) *
                                    (U(i, j, k) + U(i - 1, j, k))) +
                        gamma * inverseDx * 0.25 *
                            (fabs(U(i, j, k) + U(i + 1, j, k)) *
                                    (U(i, j, k) - U(i + 1, j, k)) +
                                fabs(U(i, j, k) + U(i - 1, j, k)) *
                                    (U(i, j, k) - U(i - 1, j, k)));

                duvdy = inverseDy * 0.25 *
                            ((V(i, j, k) + V(i + 1, j, k)) *
                                    (U(i, j, k) + U(i, j + 1, k)) -
                                (V(i, j - 1, k) + V(i + 1, j - 1, k)) *
                                    (U(i, j, k) + U(i, j - 1, k))) +
                        gamma * inverseDy * 0.25 *
                            (fabs(V(i, j, k) + V(i + 1, j, k)) *
                                    (U(i, j, k) - U(i, j + 1, k)) +
                                fabs(V(i, j - 1, k) + V(i + 1, j - 1, k)) *
                                    (U(i, j, k) - U(i, j - 1, k)));

                duwdz = inverseDz * 0.25 *
                            ((W(i, j, k) + W(i + 1, j, k)) *
                                    (U(i, j, k) + U(i, j, k + 1)) -
                                (W(i, j, k - 1) + W(i + 1, j, k - 1)) *
                                    (U(i, j, k) + U(i, j, k - 1))) +
                        gamma * inverseDz * 0.25 *
                            (fabs(W(i, j, k) + W(i + 1, j, k)) *
                                    (U(i, j, k) - U(i, j, k + 1)) +
                                fabs(W(i, j, k - 1) + W(i + 1, j, k - 1)) *
                                    (U(i, j, k) - U(i, j, k - 1)));

                du2dx2 = inverseDx * inverseDx *
                         (U(i + 1, j, k) - 2.0 * U(i, j, k) + U(i - 1, j, k));
                du2dy2 = inverseDy * inverseDy *
                         (U(i, j + 1, k) - 2.0 * U(i, j, k) + U(i, j - 1, k));
                du2dz2 = inverseDz * inverseDz *
                         (U(i, j, k + 1) - 2.0 * U(i, j, k) + U(i, j, k - 1));
                F(i, j, k) = U(i, j, k) + dt * (inverseRe * (du2dx2 + du2dy2 + du2dz2) -
                                                   du2dx - duvdy - duwdz + gx);

                duvdx = inverseDx * 0.25 *
                            ((U(i, j, k) + U(i, j + 1, k)) *
                                    (V(i, j, k) + V(i + 1, j, k)) -
                                (U(i - 1, j, k) + U(i - 1, j + 1, k)) *
                                    (V(i, j, k) + V(i - 1, j, k))) +
                        gamma * inverseDx * 0.25 *
                            (fabs(U(i, j, k) + U(i, j + 1, k)) *
                                    (V(i, j, k) - V(i + 1, j, k)) +
                                fabs(U(i - 1, j, k) + U(i - 1, j + 1, k)) *
                                    (V(i, j, k) - V(i - 1, j, k)));

                dv2dy = inverseDy * 0.25 *
                            ((V(i, j, k) + V(i, j + 1, k)) *
                                    (V(i, j, k) + V(i, j + 1, k)) -
                                (V(i, j, k) + V(i, j - 1, k)) *
                                    (V(i, j, k) + V(i, j - 1, k))) +
                        gamma * inverseDy * 0.25 *
                            (fabs(V(i, j, k) + V(i, j + 1, k)) *
                                    (V(i, j, k) - V(i, j + 1, k)) +
                                fabs(V(i, j, k) + V(i, j - 1, k)) *
                                    (V(i, j, k) - V(i, j - 1, k)));

                dvwdz = inverseDz * 0.25 *
                            ((W(i, j, k) + W(i, j + 1, k)) *
                                    (V(i, j, k) + V(i, j, k + 1)) -
                                (W(i, j, k - 1) + W(i, j + 1, k - 1)) *
                                    (V(i, j, k) + V(i, j, k + 1))) +
                        gamma * inverseDz * 0.25 *
                            (fabs(W(i, j, k) + W(i, j + 1, k)) *
                                    (V(i, j, k) - V(i, j, k + 1)) +
                                fabs(W(i, j, k - 1) + W(i, j + 1, k - 1)) *
                                    (V(i, j, k) - V(i, j, k + 1)));

                dv2dx2 = inverseDx * inverseDx *
                         (V(i + 1, j, k) - 2.0 * V(i, j, k) + V(i - 1, j, k));
                dv2dy2 = inverseDy * inverseDy *
                         (V(i, j + 1, k) - 2.0 * V(i, j, k) + V(i, j - 1, k));
                dv2dz2 = inverseDz * inverseDz *
                         (V(i, j, k + 1) - 2.0 * V(i, j, k) + V(i, j, k - 1));
                G(i, j, k) = V(i, j, k) + dt * (inverseRe * (dv2dx2 + dv2dy2 + dv2dz2) -
                                                   duvdx - dv2dy - dvwdz + gy);

                duwdx = inverseDx * 0.25 *
                            ((U(i, j, k) + U(i, j, k + 1)) *
                                    (W(i, j, k) + W(i + 1, j, k)) -
                                (U(i - 1, j, k) + U(i - 1, j, k + 1)) *
                                    (W(i, j, k) + W(i - 1, j, k))) +
                        gamma * inverseDx * 0.25 *
                            (fabs(U(i, j, k) + U(i, j, k + 1)) *
                                    (W(i, j, k) - W(i + 1, j, k)) +
                                fabs(U(i - 1, j, k) + U(i - 1, j, k + 1)) *
                                    (W(i, j, k) - W(i - 1, j, k)));

                dvwdy = inverseDy * 0.25 *
                            ((V(i, j, k) + V(i, j, k + 1)) *
                                    (W(i, j, k) + W(i, j + 1, k)) -
                                (V(i, j - 1, k + 1) + V(i, j - 1, k)) *
                                    (W(i, j, k) + W(i, j - 1, k))) +
                        gamma * inverseDy * 0.25 *
                            (fabs(V(i, j, k) + V(i, j, k + 1)) *
                                    (W(i, j, k) - W(i, j + 1, k)) +
                                fabs(V(i, j - 1, k + 1) + V(i, j - 1, k)) *
                                    (W(i, j, k) - W(i, j - 1, k)));

                dw2dz = inverseDz * 0.25 *
                            ((W(i, j, k) + W(i, j, k + 1)) *
                                    (W(i, j, k) + W(i, j, k + 1)) -
                                (W(i, j, k) + W(i, j, k - 1)) *
                                    (W(i, j, k) + W(i, j, k - 1))) +
                        gamma * inverseDz * 0.25 *
                            (fabs(W(i, j, k) + W(i, j, k + 1)) *
                                    (W(i, j, k) - W(i, j, k + 1)) +
                                fabs(W(i, j, k) + W(i, j, k - 1)) *
                                    (W(i, j, k) - W(i, j, k - 1)));

                dw2dx2 = inverseDx * inverseDx *
                         (W(i + 1, j, k) - 2.0 * W(i, j, k) + W(i - 1, j, k));
                dw2dy2 = inverseDy * inverseDy *
                         (W(i, j + 1, k) - 2.0 * W(i, j, k) + W(i, j - 1, k));
                dw2dz2 = inverseDz * inverseDz *
                         (W(i, j, k + 1) - 2.0 * W(i, j, k) + W(i, j, k - 1));
                H(i, j, k) = W(i, j, k) + dt * (inverseRe * (dw2dx2 + dw2dy2 + dw2dz2) -
                                                   duwdx - dvwdy - dw2dz + gz);
            }
        }
    }

    /* ----------------------------- boundary of F ---------------------------
     */
    if (commIsBoundary(&d->comm, LEFT)) {
        for (int k = 1; k < kmaxLocal + 1; k++) {
            for (int j = 1; j < jmaxLocal + 1; j++) {
                F(0, j, k) = U(0, j, k);
            }
        }
    }

    if (commIsBoundary(&d->comm, RIGHT)) {
        for (int k = 1; k < kmaxLocal + 1; k++) {
            for (int j = 1; j < jmaxLocal + 1; j++) {
                F(imaxLocal, j, k) = U(imaxLocal, j, k);
            }
        }
    }

    /* ----------------------------- boundary of G ---------------------------
     */
    if (commIsBoundary(&d->comm, BOTTOM)) {
        for (int k = 1; k < kmaxLocal + 1; k++) {
            for (int i = 1; i < imaxLocal + 1; i++) {
                G(i, 0, k) = V(i, 0, k);
            }
        }
    }

    if (commIsBoundary(&d->comm, TOP)) {
        for (int k = 1; k < kmaxLocal + 1; k++) {
            for (int i = 1; i < imaxLocal + 1; i++) {
                G(i, jmaxLocal, k) = V(i, jmaxLocal, k);
            }
        }
    }

    /* ----------------------------- boundary of H ---------------------------
     */
    if (commIsBoundary(&d->comm, FRONT)) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
            for (int i = 1; i < imaxLocal + 1; i++) {
                H(i, j, 0) = W(i, j, 0);
            }
        }
    }

    if (commIsBoundary(&d->comm, BACK)) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
            for (int i = 1; i < imaxLocal + 1; i++) {
                H(i, j, kmaxLocal) = W(i, j, kmaxLocal);
            }
        }
    }
}

void adaptUV(Discretization* d)
{
    int imaxLocal = d->comm.imaxLocal;
    int jmaxLocal = d->comm.jmaxLocal;
    int kmaxLocal = d->comm.kmaxLocal;

    double* p = d->p;
    double* u = d->u;
    double* v = d->v;
    double* w = d->w;
    double* f = d->f;
    double* g = d->g;
    double* h = d->h;

    double factorX = d->dt / d->grid.dx;
    double factorY = d->dt / d->grid.dy;
    double factorZ = d->dt / d->grid.dz;

    for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
            for (int i = 1; i < imaxLocal + 1; i++) {
                U(i, j, k) = F(i, j, k) - (P(i + 1, j, k) - P(i, j, k)) * factorX;
                V(i, j, k) = G(i, j, k) - (P(i, j + 1, k) - P(i, j, k)) * factorY;
                W(i, j, k) = H(i, j, k) - (P(i, j, k + 1) - P(i, j, k)) * factorZ;
            }
        }
    }
}




