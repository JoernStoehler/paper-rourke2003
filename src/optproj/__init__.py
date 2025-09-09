try:
    from . import _cext as _cext_module
except Exception as exc:
    raise ImportError(
        "C extension not built. Run `make build` to compile."
    ) from exc

ffi = _cext_module.ffi
lib = _cext_module.lib

def l2_norm(arr):
    """Compute L2 norm using the fast C hotspot."""
    import numpy as _np
    a = _np.ascontiguousarray(arr, dtype=_np.float64)
    ptr = ffi.cast("double *", ffi.from_buffer(a))
    return lib.l2_norm(ptr, a.size)
