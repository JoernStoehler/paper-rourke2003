# Convex Partitions of the Square

This project is a numerical exploration of convex partitions of the square, with the goal of finding a partition that minimizes the maximum aspect ratio of its constituent pieces. The aspect ratio is defined as the circumradius divided by the inradius.

This work is an attempt to replicate and improve upon the results of the 2003 paper by Rourke et al.

## Quickstart

```bash
make
make test
make jupyter   # open JupyterLab
```

This will create a virtual environment in `.venv` and install the required dependencies.

## Layout
src/optproj/numerics.c # C source for numerical computations
src/optproj/numerics.h # C header for numerical computations
src/optproj/         # Python package
scripts/             # build script
tests/               # pytest tests
build/               # intermediate files (ignored)
