# BCNNeuroModel

A C++ neuromodeling project with Eigen and gflags dependencies.

## Building

```bash
mkdir Build && cd Build
cmake ..
make
```

## Testing

Tests are built automatically. Run them with:

```bash
cd Build
ctest
```

## Dependencies

- **Eigen3**: Linear algebra library (automatically fetched if not found)
- **gflags**: Command-line flag parsing library (must be installed)
- **GoogleTest**: Unit testing framework (automatically fetched)

### Dependency Resolution

- **Eigen**: The build system checks in order:
  1. Local copy in `src/Eigen/`
  2. System-installed Eigen3 via CMake's find_package
  3. Automatic fetch from GitLab (version 3.4.0) if not found

- **gflags**: Must be installed on your system. Install via:
  - macOS: `brew install gflags`
  - Ubuntu/Debian: `sudo apt-get install libgflags-dev`
  - Or build from source

- **GoogleTest**: Automatically fetched during build, no installation needed
