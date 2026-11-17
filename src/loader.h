#ifndef LOADER_H
#define LOADER_H
typedef struct {
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
    double m;
} Body;
int loader_read(const char *path, Body **out, int *nout);
void loader_free(Body *b);
int loader_write_snapshot(const char *path, const double *x, const double *y, const double *z, const double *vx, const double *vy, const double *vz, const double *m, int n);
int loader_write_energy(const char *path, double energy);
#endif
