# optproj (minimal scaffold)

Minimal Python+C scaffold using cffi.  
No hidden magic, all build products isolated in `build/`.

## Quickstart
pip install -r requirements.txt
make build
make test
make jupyter   # open JupyterLab

## Layout
src/optproj/c/       # C sources and headers
src/optproj/         # Python package
scripts/             # build script
tests/               # pytest tests
build/               # intermediate files (ignored)

## Notes
- The `build/` directory may contain compiler/intermediate clutter due to cffi/gcc quirks. This is normal for minimal cffi builds.
- Run `make clean` to remove all build products, intermediates, and Python bytecode caches (`__pycache__`, `.pyc`).