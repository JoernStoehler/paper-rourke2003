#!/usr/bin/env python3
import os
import glob
import shutil
from cffi import FFI  # type: ignore[import]

HERE = os.path.dirname(__file__)
ROOT = os.path.abspath(os.path.join(HERE, ".."))
SRC_ROOT = os.path.join(ROOT, "src")
PKG = "optproj"
C_DIR = os.path.join(SRC_ROOT, PKG, "c")
BUILD_DIR = os.path.join(ROOT, "build", PKG)

ffibuilder = FFI()
ffibuilder.cdef("""
    double l2_norm(const double *arr, int n);
""")

ffibuilder.set_source(
    f"{PKG}._cext",
    '#include "norm.h"\n',
    sources=[os.path.join(C_DIR, "norm.c")],
    include_dirs=[C_DIR],
)

if __name__ == "__main__":
    os.makedirs(BUILD_DIR, exist_ok=True)
    ffibuilder.compile(tmpdir=BUILD_DIR, verbose=True)
    # cffi puts the .so in build/optproj/optproj/ (due to module name)
    so_files = glob.glob(os.path.join(BUILD_DIR, PKG, "_cext.*.so"))
    if so_files:
        shutil.move(so_files[0], os.path.join(SRC_ROOT, PKG, os.path.basename(so_files[0])))
