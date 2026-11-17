CC=cc
MPICC=mpicc
CFLAGS=-O3 -march=native -std=c11 -Wall
LDLIBS=-lm
BINDIR=bin
all: $(BINDIR)/nbody_serial $(BINDIR)/nbody_mpi $(BINDIR)/nbody_mpi_omp
$(BINDIR)/nbody_serial: nbody_serial.c loader.c physics.c | $(BINDIR)
	$(CC) $(CFLAGS) -o $@ nbody_serial.c loader.c physics.c $(LDLIBS)
$(BINDIR)/nbody_mpi: nbody_mpi.c loader.c physics.c | $(BINDIR)
	$(MPICC) $(CFLAGS) -o $@ nbody_mpi.c loader.c physics.c $(LDLIBS)
$(BINDIR)/nbody_mpi_omp: nbody_mpi.c loader.c physics.c | $(BINDIR)
	$(MPICC) $(CFLAGS) -fopenmp -o $@ nbody_mpi.c loader.c physics.c $(LDLIBS)
$(BINDIR):
	mkdir -p $(BINDIR)
clean:
	rm -rf $(BINDIR)
