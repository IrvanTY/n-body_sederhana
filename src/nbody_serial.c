#include <stdio.h>
#include <stdlib.h>
#include "loader.h"
#include "physics.h"
int main(int argc, char **argv) {
    if (argc != 6) {
        fprintf(stderr, "usage: %s data steps dt snap_out energy_out\n", argv[0]);
        return 1;
    }
    Body *bodies = 0;
    int n = 0;
    if (loader_read(argv[1], &bodies, &n) != 0) {
        return 2;
    }
    int steps = atoi(argv[2]);
    double dt = atof(argv[3]);
    State *s = state_alloc(n, n);
    if (!s) {
        loader_free(bodies);
        return 3;
    }
    for (int i = 0; i < n; i++) {
        s->x[i] = bodies[i].x;
        s->y[i] = bodies[i].y;
        s->z[i] = bodies[i].z;
        s->vx[i] = bodies[i].vx;
        s->vy[i] = bodies[i].vy;
        s->vz[i] = bodies[i].vz;
        s->m[i] = bodies[i].m;
    }
    for (int step = 0; step < steps; step++) {
        physics_forces(s->x, s->y, s->z, s->m, n, 0, n, s->fx, s->fy, s->fz);
        physics_kick(s->vx, s->vy, s->vz, s->fx, s->fy, s->fz, dt, 0, n);
        physics_drift(s->x, s->y, s->z, s->vx, s->vy, s->vz, dt, 0, n);
        physics_forces(s->x, s->y, s->z, s->m, n, 0, n, s->fx, s->fy, s->fz);
        physics_kick(s->vx, s->vy, s->vz, s->fx, s->fy, s->fz, dt, 0, n);
    }
    double kinetic = 0.0;
    double potential = 0.0;
    physics_energy_split(s->x, s->y, s->z, s->vx, s->vy, s->vz, s->m, n, 0, n, &kinetic, &potential);
    int status = 0;
    if (loader_write_snapshot(argv[4], s->x, s->y, s->z, s->vx, s->vy, s->vz, s->m, n) != 0) {
        status = 4;
    } else if (loader_write_energy(argv[5], kinetic + potential) != 0) {
        status = 5;
    }
    state_free(s);
    loader_free(bodies);
    return status;
}
