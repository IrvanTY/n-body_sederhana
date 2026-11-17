#include "physics.h"
#include <math.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif
const double PHYSICS_G = 0.01720209895 * 0.01720209895;
const double PHYSICS_EPS2 = 1e-12;
State *state_alloc(int n, int force_n) {
    State *s = malloc(sizeof(State));
    if (!s) {
        return 0;
    }
    int fn = force_n > 0 ? force_n : 1;
    s->x = malloc(sizeof(double) * (size_t)n);
    s->y = malloc(sizeof(double) * (size_t)n);
    s->z = malloc(sizeof(double) * (size_t)n);
    s->vx = malloc(sizeof(double) * (size_t)n);
    s->vy = malloc(sizeof(double) * (size_t)n);
    s->vz = malloc(sizeof(double) * (size_t)n);
    s->m = malloc(sizeof(double) * (size_t)n);
    s->fx = malloc(sizeof(double) * (size_t)fn);
    s->fy = malloc(sizeof(double) * (size_t)fn);
    s->fz = malloc(sizeof(double) * (size_t)fn);
    if (!s->x || !s->y || !s->z || !s->vx || !s->vy || !s->vz || !s->m || !s->fx || !s->fy || !s->fz) {
        state_free(s);
        return 0;
    }
    return s;
}
void state_free(State *s) {
    if (!s) {
        return;
    }
    free(s->x);
    free(s->y);
    free(s->z);
    free(s->vx);
    free(s->vy);
    free(s->vz);
    free(s->m);
    free(s->fx);
    free(s->fy);
    free(s->fz);
    free(s);
}
void physics_forces(const double *x, const double *y, const double *z, const double *m, int n, int start, int count, double *fx, double *fy, double *fz) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int k = 0; k < count; k++) {
        int i = start + k;
        double xi = x[i];
        double yi = y[i];
        double zi = z[i];
        double fxi = 0.0;
        double fyi = 0.0;
        double fzi = 0.0;
        for (int j = 0; j < n; j++) {
            if (i == j) {
                continue;
            }
            double dx = x[j] - xi;
            double dy = y[j] - yi;
            double dz = z[j] - zi;
            double r2 = dx * dx + dy * dy + dz * dz + PHYSICS_EPS2;
            double inv = 1.0 / (sqrt(r2) * r2);
            double s = PHYSICS_G * m[j] * inv;
            fxi += s * dx;
            fyi += s * dy;
            fzi += s * dz;
        }
        fx[k] = fxi;
        fy[k] = fyi;
        fz[k] = fzi;
    }
}
void physics_kick(double *vx, double *vy, double *vz, const double *fx, const double *fy, const double *fz, double dt, int start, int count) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int k = 0; k < count; k++) {
        int i = start + k;
        vx[i] += 0.5 * fx[k] * dt;
        vy[i] += 0.5 * fy[k] * dt;
        vz[i] += 0.5 * fz[k] * dt;
    }
}
void physics_drift(double *x, double *y, double *z, const double *vx, const double *vy, const double *vz, double dt, int start, int count) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int k = 0; k < count; k++) {
        int i = start + k;
        x[i] += vx[i] * dt;
        y[i] += vy[i] * dt;
        z[i] += vz[i] * dt;
    }
}
void physics_energy_split(const double *x, const double *y, const double *z, const double *vx, const double *vy, const double *vz, const double *m, int n, int start, int count, double *kinetic, double *potential) {
    double ke = 0.0;
    double pe = 0.0;
    for (int k = 0; k < count; k++) {
        int i = start + k;
        ke += 0.5 * m[i] * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
        for (int j = i + 1; j < n; j++) {
            double dx = x[j] - x[i];
            double dy = y[j] - y[i];
            double dz = z[j] - z[i];
            double r = sqrt(dx * dx + dy * dy + dz * dz + PHYSICS_EPS2);
            pe -= PHYSICS_G * m[i] * m[j] / r;
        }
    }
    *kinetic = ke;
    *potential = pe;
}
