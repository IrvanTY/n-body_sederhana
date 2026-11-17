#!/bin/bash
#SBATCH -J nbody
#SBATCH -N 2
#SBATCH -n 16
#SBATCH --ntasks-per-node=8
#SBATCH -c 1
#SBATCH -o nbody_%j.out
set -e
ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}"
SIF="${SIF:-$ROOT/container/nbody.sif}"
DATA="${DATA:-$ROOT/data/asteroids.csv}"
STEPS="${STEPS:-100}"
DT="${DT:-0.5}"
SNAP_OUT="${SNAP_OUT:-$ROOT/results/snap.csv}"
ENERGY_OUT="${ENERGY_OUT:-$ROOT/results/energy.txt}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
srun apptainer exec --bind /home:/home "$SIF" /opt/nbody/nbody_mpi_omp "$DATA" "$STEPS" "$DT" "$SNAP_OUT" "$ENERGY_OUT"
