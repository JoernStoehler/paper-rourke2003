# GEMINI Project Guide

This document provides a guide for AI developers to understand and contribute to this project.

## 1. Project Overview

This project aims to replicate and improve upon the 2003 paper by Rourke et al. on convex partitions of the square. The primary goal is to find a convex partition with a minimal maximum aspect ratio (circumradius divided by inradius) for its constituent pieces.

The project employs a hybrid approach:
- **C**: For high-performance numerical computations (numerical lifting).
- **Python**: For rapid development, orchestration, analysis, and testing.
- **Jupyter**: For interactive experiments and visualization.

The methodology involves numerically optimizing the positions of vertices within the square to discover partitions with low aspect ratios. These numerical examples are intended to inspire an explicit mathematical construction and, ultimately, a proof of optimality.

## 2. Folder Structure

The project follows a canonical structure to ensure predictability and ease of navigation:

- **`/src/optproj/`**: The main Python package source code.
- **`/src/optproj/numerics.c`**: C source file for numerical computations.
- **`/src/optproj/numerics.h`**: C header file for numerical computations.
- **`/src/optproj/numerics.py`**: Python wrapper for the C functions.
- **`/tests/`**: Unit tests for the Python and C components. All test files should be prefixed with `test_`.
- **`/scripts/`**: Build and utility scripts.
- **`/build/`**: Intermediate and final build artifacts. This directory is ignored by Git.
- **`/Makefile`**: Defines common project commands.
- **`/requirements.txt`**: A list of Python dependencies.

## 3. Common Commands

The `Makefile` provides a set of commands for common development tasks:

- **`make install`**: Creates a virtual environment and installs dependencies from `requirements.txt`.
- **`make build`**: Compiles the C extension using `cffi`.
- **`make test`**: Runs the `pytest` test suite.
- **`make clean`**: Removes all build artifacts and temporary files.

## 4. Style Conventions

To maintain code quality and consistency, please adhere to the following style conventions:

### 4.1. Python

- **Functional Programming**: Where appropriate, favor a functional programming style. Use pure functions, avoid side effects, and leverage functions from libraries like `itertools` and `functools`.
- **Type Annotations**: All new Python code **must** include type annotations.
- **Shape Annotations**: For all code involving `numpy` arrays or other tensor-like objects, use `jaxtyping` to provide shape annotations. This is critical for ensuring the correctness of numerical computations.

  ```python
  from jaxtyping import Float, Array
  import numpy as np

  def my_function(x: Float[Array, "batch channels"]) -> Float[Array, "batch"]:
      # function implementation
      ...
  ```

### 4.2. Unit Tests

- All new functionality must be accompanied by unit tests.
- Tests should be placed in the `/tests` directory.
- Use `pytest` for writing and running tests.

## 5. Onboarding Wishlist

This section is a living document where new developers can add what they would have liked to see during their onboarding.

- *Add your suggestions here.*

## 6. Further Reading

For a deeper understanding of the project's mathematical background, goals, and current progress, please refer to the following documents:

- **`README.md`**: High-level project description.
- **`docs/`**: Detailed mathematical formulations, roadmaps, and progress reports. (Coming soon)
