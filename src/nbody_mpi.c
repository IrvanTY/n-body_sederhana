#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <limits.h>
#include "loader.h"
#include "physics.h"

int block_counts(int n, int p, int rank) {
    return n / p + (rank < n % p ? 1 : 0);
}

int block_start(int n, int p, int rank) {
    int base = n / p;
    int rem = n % p;
    return rank * base + (rank < rem ? rank : rem);
}

void pack_triple(double *buf, const double *x, const double *y, const double *z, int start, int count) {
    for (int k = 0; k < count; k++) {
        int i = start + k;
        buf[3 * k] = x[i];
        buf[3 * k + 1] = y[i];
        buf[3 * k + 2] = z[i];
    }
}

void unpack_triple(const double *buf, double *x, double *y, double *z, int n) {
    for (int i = 0; i < n; i++) {
        x[i] = buf[3 * i];
        y[i] = buf[3 * i + 1];
        z[i] = buf[3 * i + 2];
    }
}

int main(int argc, char **argv) {
    if (argc != 6) {
        fprintf(stderr, "usage: %s data steps dt snap_out energy_out\n", argv[0]);
        return 1;
    }
    int init_rc = MPI_Init(&argc, &argv);
    if (init_rc != MPI_SUCCESS) {
        fprintf(stderr, "%s: MPI_Init failed\n", argv[0]);
        return 1;
    }
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#ifdef _OPENMP
    if (rank == 0) {
        fprintf(stderr, "%s: threads=%d\n", argv[0], omp_get_max_threads());
    }
#endif
    Body *all = 0;
    int n = 0;
    if (rank == 0) {
        if (loader_read(argv[1], &all, &n) != 0) {
            MPI_Abort(MPI_COMM_WORLD, 2);
        }
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    size_t body_bytes = sizeof(Body) * (size_t)n;
    if (body_bytes > (size_t)INT_MAX) {
        fprintf(stderr, "%s: body count exceeds MPI int range\n", argv[0]);
        MPI_Abort(MPI_COMM_WORLD, 2);
    }
    if (rank != 0) {
        all = malloc(sizeof(Body) * (size_t)n);
        if (!all) {
            fprintf(stderr, "%s: allocation failed\n", argv[0]);
            MPI_Abort(MPI_COMM_WORLD, 3);
        }
    }
    MPI_Bcast(all, (int)body_bytes, MPI_BYTE, 0, MPI_COMM_WORLD);
    int steps = atoi(argv[2]);
    double dt = atof(argv[3]);
    int local_n = block_counts(n, size, rank);
    int start = block_start(n, size, rank);
    double t_comp = 0.0;
    double t_comm = 0.0;
    double tick = 0.0;
    State *s = state_alloc(n, local_n);
    double *pos = malloc(sizeof(double) * (size_t)3 * (size_t)n);
    double *vel = malloc(sizeof(double) * (size_t)3 * (size_t)n);
    int *counts = malloc(sizeof(int) * (size_t)size);
    int *displs = malloc(sizeof(int) * (size_t)size);
    double *send = malloc(sizeof(double) * (size_t)3 * (size_t)(local_n > 0 ? local_n : 1));
    if (!s || !pos || !vel || !counts || !displs || !send) {
        fprintf(stderr, "%s: allocation failed\n", argv[0]);
        state_free(s);
        free(pos);
        free(vel);
        free(counts);
        free(displs);
        free(send);
        if (rank == 0) {
            loader_free(all);
        } else {
            free(all);
        }
        MPI_Abort(MPI_COMM_WORLD, 3);
    }
    for (int i = 0; i < n; i++) {
        s->x[i] = all[i].x;
        s->y[i] = all[i].y;
        s->z[i] = all[i].z;
        s->vx[i] = all[i].vx;
        s->vy[i] = all[i].vy;
        s->vz[i] = all[i].vz;
        s->m[i] = all[i].m;
    }
    for (int r = 0; r < size; r++) {
        counts[r] = 3 * block_counts(n, size, r);
        displs[r] = 3 * block_start(n, size, r);
    }
    for (int step = 0; step < steps; step++) {
        tick = MPI_Wtime();
        physics_forces(s->x, s->y, s->z, s->m, n, start, local_n, s->fx, s->fy, s->fz);
        physics_kick(s->vx, s->vy, s->vz, s->fx, s->fy, s->fz, dt, start, local_n);
        physics_drift(s->x, s->y, s->z, s->vx, s->vy, s->vz, dt, start, local_n);
        pack_triple(send, s->x, s->y, s->z, start, local_n);
        t_comp += MPI_Wtime() - tick;
        tick = MPI_Wtime();
        MPI_Allgatherv(send, 3 * local_n, MPI_DOUBLE, pos, counts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
        t_comm += MPI_Wtime() - tick;
        tick = MPI_Wtime();
        unpack_triple(pos, s->x, s->y, s->z, n);
        physics_forces(s->x, s->y, s->z, s->m, n, start, local_n, s->fx, s->fy, s->fz);
        physics_kick(s->vx, s->vy, s->vz, s->fx, s->fy, s->fz, dt, start, local_n);
        pack_triple(send, s->vx, s->vy, s->vz, start, local_n);
        t_comp += MPI_Wtime() - tick;
        tick = MPI_Wtime();
        MPI_Allgatherv(send, 3 * local_n, MPI_DOUBLE, vel, counts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
        t_comm += MPI_Wtime() - tick;
        tick = MPI_Wtime();
        unpack_triple(vel, s->vx, s->vy, s->vz, n);
        t_comp += MPI_Wtime() - tick;
    }
    double epart[2];
    double eglob[2];
    epart[0] = 0.0;
    epart[1] = 0.0;
    eglob[0] = 0.0;
    eglob[1] = 0.0;
    tick = MPI_Wtime();
    physics_energy_split(s->x, s->y, s->z, s->vx, s->vy, s->vz, s->m, n, start, local_n, &epart[0], &epart[1]);
    t_comp += MPI_Wtime() - tick;
    tick = MPI_Wtime();
    MPI_Allreduce(epart, eglob, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    t_comm += MPI_Wtime() - tick;
    double rank_times[3];
    double max_times[3];
    rank_times[0] = t_comp;
    rank_times[1] = t_comm;
    rank_times[2] = t_comp + t_comm;
    MPI_Reduce(rank_times, max_times, 3, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        printf("TIMING t_comp=%.6f t_comm=%.6f t_sum=%.6f\n", max_times[0], max_times[1], max_times[2]);
        fflush(stdout);
    }
    int status = 0;
    if (rank == 0) {
        if (loader_write_snapshot(argv[4], s->x, s->y, s->z, s->vx, s->vy, s->vz, s->m, n) != 0) {
            status = 4;
        } else if (loader_write_energy(argv[5], eglob[0] + eglob[1]) != 0) {
            status = 5;
        }
    }
    MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    free(send);
    free(counts);
    free(displs);
    free(pos);
    free(vel);
    state_free(s);
    if (rank == 0) {
        loader_free(all);
    } else {
        free(all);
    }
    MPI_Finalize();
    return status;
}
