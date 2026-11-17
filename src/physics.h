#ifndef PHYSICS_H
#define PHYSICS_H
typedef struct {
    double *x;
    double *y;
    double *z;
    double *vx;
    double *vy;
    double *vz;
    double *m;
    double *fx;
    double *fy;
    double *fz;
} State;
State *state_alloc(int n, int force_n);
void state_free(State *s);
extern const double PHYSICS_G;
extern const double PHYSICS_EPS2;
void physics_forces(const double *x, const double *y, const double *z, const double *m, int n, int start, int count, double *fx, double *fy, double *fz);
void physics_kick(double *vx, double *vy, double *vz, const double *fx, const double *fy, const double *fz, double dt, int start, int count);
void physics_drift(double *x, double *y, double *z, const double *vx, const double *vy, const double *vz, double dt, int start, int count);
void physics_energy_split(const double *x, const double *y, const double *z, const double *vx, const double *vy, const double *vz, const double *m, int n, int start, int count, double *kinetic, double *potential);
#endif
