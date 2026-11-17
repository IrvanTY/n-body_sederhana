#include "loader.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
void loader_free(Body *b) {
    free(b);
}
int loader_read(const char *path, Body **out, int *nout) {
    FILE *f = fopen(path, "r");
    if (!f) return 1;
    char line[1024];
    if (!fgets(line, sizeof(line), f)) {
        fclose(f);
        return 2;
    }
    int cap = 2048;
    int n = 0;
    Body *b = malloc(sizeof(Body) * (size_t)cap);
    if (!b) {
        fclose(f);
        return 3;
    }
    while (fgets(line, sizeof(line), f)) {
        double x, y, z, vx, vy, vz, m;
        if (sscanf(line, "%lf,%lf,%lf,%lf,%lf,%lf,%lf", &x, &y, &z, &vx, &vy, &vz, &m) != 7) continue;
        if (n >= cap) {
            cap *= 2;
            Body *nb = realloc(b, sizeof(Body) * (size_t)cap);
            if (!nb) {
                free(b);
                fclose(f);
                return 4;
            }
            b = nb;
        }
        b[n].x = x;
        b[n].y = y;
        b[n].z = z;
        b[n].vx = vx;
        b[n].vy = vy;
        b[n].vz = vz;
        b[n].m = m;
        n++;
    }
    fclose(f);
    *out = b;
    *nout = n;
    return 0;
}
int loader_write_snapshot(const char *path, const double *x, const double *y, const double *z, const double *vx, const double *vy, const double *vz, const double *m, int n) {
    FILE *f = fopen(path, "w");
    if (!f) {
        return 1;
    }
    for (int i = 0; i < n; i++) {
        fprintf(f, "%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n", x[i], y[i], z[i], vx[i], vy[i], vz[i], m[i]);
    }
    fclose(f);
    return 0;
}
int loader_write_energy(const char *path, double energy) {
    FILE *f = fopen(path, "w");
    if (!f) {
        return 1;
    }
    fprintf(f, "%.12e\n", energy);
    fclose(f);
    return 0;
}
