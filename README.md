![Magnoom](https://github.com/n-s-kiselev/magnoom/blob/master/MagnoomWiki/TitleImage.png)

# Magnoom

Magnoom is software for atomistic spin-dynamics simulations with real-time OpenGL visualization and interactive parameter control. See the [project wiki](https://github.com/n-s-kiselev/magnoom/wiki) for additional documentation and examples.

## Main features

- Portable C/C++ code for Linux, macOS, and Windows (MinGW)
- Multithreaded spin-dynamics calculations
- Real-time parameter control through AntTweakBar
- OpenGL visualization and post-processing tools, including slicing and filtering
- Import and export support for formats used by
  [OOMMF](https://math.nist.gov/oommf/) and
  [MuMax3](https://mumax.github.io/), including OVF

## Building

Magnoom uses a single cross-platform [`nob.c`](https://github.com/tsoding/nob.h) build program. Bootstrap it once from the repository root:

```sh
cc nob.c -o nob
```

Then run:

```sh
./nob          # build Magnoom and its vendored dependencies
./nob -clean   # remove all generated build output
./nob -help    # list the supported options
```

On Windows with MinGW, use:

```sh
gcc nob.c -o nob.exe
./nob.exe
```

The build produces:

- `build/magnoom` (`build/magnoom.exe` on Windows) — the application
- `build/libAntTweakBar.a` — the statically linked AntTweakBar library
- intermediate objects under `build/`

The `nob` executable automatically rebuilds itself when `nob.c` or its vendored `nob.h` changes. All normal build products are contained in `build/`, except for the bootstrapped `nob` executable itself.

`nob.c` is the only build system for this repository. No CMake, Makefile, Visual Studio project, package-manager build, or dependency download is required.

## Dependencies

The complete source for the application dependencies is stored under [`vendor`](vendor):

- [AntTweakBar Legacy](https://github.com/n-s-kiselev/AntTweakBar-Legacy) — the real-time OpenGL parameter interface, trimmed to the source needed by Magnoom and built as a static library from `vendor/AntTweakBar-Legacy`
- GLFW 2.7.9 — window creation, input, and the event loop, built from `vendor/glfw2`
- [GLAD](https://glad.dav1d.de/) — OpenGL function loading, built from `vendor/glad`
- [`nob.h`](https://github.com/tsoding/nob.h) — build-system implementation, stored in `vendor/nob`

AntTweakBar's example programs are intentionally not included. GLFW and GLAD are built directly from the vendored source, so no system installation of GLFW, GLAD, or AntTweakBar is needed. Magnoom uses an OpenGL compatibility context because its renderer uses the fixed-function API.

The remaining platform dependencies provide the compiler, system OpenGL libraries, native window-system headers, and threading support:

### Linux

A C/C++ compiler, `ar`, OpenGL, X11, and XRandR development files are required. On Debian-based Linux:

```sh
sudo apt update
sudo apt install build-essential libgl1-mesa-dev libx11-dev libxrandr-dev
```

### macOS

Install the Xcode Command Line Tools:

```sh
xcode-select --install
```

The build uses the system OpenGL, Cocoa, AppKit, Foundation, IOKit, and
CoreVideo frameworks. No Homebrew packages are required.

### Windows

Use a MinGW-w64 toolchain that provides `gcc`, `g++`, and `ar`. The build
links against the Windows OpenGL, GDI, and multimedia system libraries.
Run the bootstrap and build commands from a MinGW-compatible terminal.

Supported platforms are Linux, macOS, and Windows with MinGW. The current
build has been verified directly on macOS ARM64; the Linux and Windows paths
are implemented but have not yet been verified in this repository.

## Authors

Developer: [Nikolai S. Kiselev](https://github.com/n-s-kiselev)

Contributors:

1. [Vladyslav M. Kuchkin](https://www.researchgate.net/profile/Vladyslav-Kuchkin)
2. [Andriy S. Savchenko](https://www.researchgate.net/profile/A_Savchenko)
3. [Filipp N. Rybakov](http://www.quantumandclassical.com/excalibur/)

## Citation

If you use Magnoom in scientific research, please cite the following publication:

A. S. Savchenko, V. M. Kuchkin, F. N. Rybakov, S. Blügel, and N. S. Kiselev,
"Chiral standing spin waves in skyrmion lattice,"
*APL Materials* **10**, 071111 (2022).
[https://doi.org/10.1063/5.0097651](https://doi.org/10.1063/5.0097651)

## License

See [`LICENSE`](LICENSE). AntTweakBar and the other vendored components retain
their own license files in their respective directories under [`vendor`](vendor).
