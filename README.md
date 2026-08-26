![Magnoom](https://github.com/n-s-kiselev/magnoom/blob/master/MagnoomWiki/TitleImage.png)

# Magnoom

Magnoom is cross-platform software for atomistic spin-dynamics simulations with real-time OpenGL visualization and interactive parameter control, running natively on Linux, macOS, and Windows. See the [project wiki](https://github.com/n-s-kiselev/magnoom/wiki) for additional documentation and examples.

## Main features

- Real-time parameter control through [AntTweakBar](https://github.com/n-s-kiselev/AntTweakBar-Legacy)
- Cross-platform C99 application code for Linux, macOS, and Windows (MinGW), built with the same single cross-platform build script on every platform
- Multithreaded solvers
- OpenGL visualization and post-processing tools, including slicing and filtering
- Import and export support for formats used by
  [OOMMF](https://math.nist.gov/oommf/) and
  [MuMax3](https://mumax.github.io/), including [OVF](https://math.nist.gov/oommf/doc/userguide20b0/userguide/Vector_Field_File_Format_OV.html)

## External magnetic field

Magnoom separates the external magnetic field into two components:

- `BextDC` is the static component, controlled by its magnitude and direction angles.
- `BextAC` is the time-dependent component, controlled by its amplitude, direction, waveform, and timing parameters.

The AC component supports sinusoidal, Gaussian pulse, sinc pulse, and circular waveforms. The F5 `BextAC` control bar also provides phase-resolved dynamical-mode recording controls. The arrow shown in the 3D view represents `BextDC`; its visual scale does not change the field used by the solvers.

## Visualization lighting

Press `l` or use the View bar's `Lighting` control to cycle between three modes:

- `Off` disables lighting for all 3D geometry.
- `Fixed` uses the configurable directional light.
- `Adaptive` places a point light at the camera so it follows the view.

## Building

Clone Magnoom with its AntTweakBar submodule, then enter the repository:

```sh
git clone --recurse-submodules https://github.com/n-s-kiselev/magnoom.git
cd magnoom
```

If you already cloned without `--recurse-submodules`, initialize AntTweakBar from
the repository root:

```sh
git submodule update --init
```

`vendor/AntTweakBar-Legacy` is a git submodule. The `nob` build tool also
initializes it automatically on the first build if it is missing.

Magnoom uses a single cross-platform [nob.h](https://github.com/tsoding/nob.h)
build system. Bootstrap it once from the repository root:

```sh
cc nob.c -o nob
```

Then run:

```sh
./nob          # build Magnoom and its vendored dependencies
./nob -clean   # remove all generated build output
./nob -test    # build and run the automated tests
./nob -help    # list the supported options
```

On Windows with [MinGW](https://www.mingw-w64.org/), use:

```sh
gcc nob.c -o nob.exe
./nob.exe
```

The build produces:

- `build/magnoom` (`build/magnoom.exe` on Windows) — the application
- `build/Magnoom.app` on macOS — an application bundle with the Magnoom icon
- `build/share/` on Linux — staged desktop-launcher and hicolor icon files
- `build/libAntTweakBar.a` — the statically linked AntTweakBar library
- intermediate objects under `build/`

On macOS, launch `build/Magnoom.app` to use the application name and Dock icon;
the raw `build/magnoom` executable remains available for terminal use. On
Windows, the icon is embedded directly in `build/magnoom.exe`. On Linux, the
window publishes its icon and `Magnoom` application class directly, while the
files under `build/share/` can be copied into a standard installation prefix
to provide a desktop launcher.

## Simulation file locations

The Initial State panel provides separate `Input directory` and `Output
directory` fields. Both start as the directory containing the running Magnoom
executable and can be edited for the current session. Relative input and output
file names are resolved against these fields; an absolute file name overrides
the corresponding directory.

The panel has one `Import` button and one `Export` button. Magnoom selects the
format from the final lowercase file extension. Import supports `.csv`, `.ovf`,
`.vtk`, and `.bin`; export additionally supports `.png`. Extensions are
case-sensitive. Before importing, Magnoom checks that CSV contains exactly the
expected bounded six-value spin rows, OVF has the expected header, dimensions,
and data marker, VTK contains a complete binary-float vector field for the
current grid, and BIN has the exact payload size required by the current grid.
Missing, unsupported, export-only, or content-mismatched extensions produce a
modal error message. Manual import stops the simulation and leaves it stopped;
export writes the current spin state without changing whether the simulation is
running.

The output directory applies to manual CSV, OVF, VTK, BIN, and PNG exports as
well as `table.csv`, `phase*.vtk`, and `dTdF*.vtk`. Changing it starts a new
`table.csv` after flushing buffered records to the previous file. Magnoom does
not create missing directories, and reports the resolved path when opening a
file fails. Paths may contain spaces and dots and are limited to 4095 bytes.
Windows paths use the MinGW C runtime's narrow-character encoding. If an
installed executable directory is read-only, select a writable output
directory before recording or exporting data. On Windows, use fully qualified
drive or UNC paths; drive-relative paths such as `C:file.csv` and current-drive
rooted paths such as `\file.csv` are rejected as ambiguous.

Magnoom's application sources are compiled as C99. The vendored AntTweakBar implementation remains C++, so the final executable is linked with the C++ compiler driver.

The active crystal basis is selected in `magnoom.c` immediately after
`readConfigFile()`. The one-atom simple-cubic basis is enabled by default;
commented B20, EuSi, FCC2, and FCC3 alternatives are provided beside it. Keep each
alternative's lattice vectors and atom positions together when selecting it.
`magnoom_ctx_set_block()` copies up to 100 Cartesian atom positions, validates
them against the half-open unit cell, and must run before neighbor maps and
spin or drawing arrays are created.

The `nob` executable automatically rebuilds itself when `nob.c` or its vendored `nob.h` changes. All normal build products are contained in `build/`, except for the bootstrapped `nob` executable itself.

`nob.c` is the only build system for this repository. No CMake, Makefile, Visual Studio project, or package-manager build is required. The only network access a build may trigger is the one-time git submodule fetch described above; there is no package-manager dependency resolution.

The optional `nob_icons.c` tool regenerates platform icon assets from
`assets/icon/magnoom.svg` using Inkscape. It is separate from the normal build;
see `assets/icon/icon_instructions.md` for usage.

## Dependencies

The complete source for the application dependencies is stored under [`vendor`](vendor):

- [AntTweakBar](https://github.com/n-s-kiselev/AntTweakBar-Legacy) — the real-time OpenGL parameter interface, included as a git submodule at `vendor/AntTweakBar-Legacy` and built as a static library from the subset of its sources Magnoom needs (see `nob.c`'s `atb_common_sources`)
- [GLFW](https://www.glfw.org/) [v*2.7.9*](https://sourceforge.net/projects/glfw/files/glfw/2.7.9/) — window creation, input, and the event loop, built from `vendor/glfw2`
- [GLAD](https://glad.dav1d.de/) — OpenGL function loading, built from `vendor/glad`
- [stb_image](https://github.com/nothings/stb) and [stb_image_write](https://github.com/nothings/stb) — PNG decoding for the optional icon exporter and PNG image export for Magnoom, stored in `vendor/stb`
- [nob.h](https://github.com/tsoding/nob.h) — build-system implementation, stored in `vendor/nob`

Magnoom uses an OpenGL compatibility context because its renderer uses the fixed-function API (no shaders).

The remaining platform dependencies provide the compiler, system OpenGL libraries, native window-system headers, and threading support:

### Linux

A C99 compiler, a C++ compiler for AntTweakBar, `ar`, OpenGL, X11, and XRandR development files are required. On Debian-based Linux:

```sh
sudo apt update
sudo apt install build-essential libgl1-mesa-dev libx11-dev libxrandr-dev
```

### macOS

Install the Xcode Command Line Tools:

```sh
xcode-select --install
```

The build uses the system OpenGL, Cocoa, AppKit, Foundation, IOKit, and CoreVideo frameworks. No Homebrew packages are required.

### Windows

Use a MinGW-w64 toolchain that provides `gcc`, `g++`, and `ar`. The build links against the Windows OpenGL, GDI, and multimedia system libraries. Run the bootstrap and build commands from a MinGW-compatible terminal.

Supported platforms are Linux, macOS, and Windows with MinGW, using the same `nob.c` build script on all three — no platform-specific build files or tooling are required. The current build has been verified directly on macOS Intel and ARM64, and on Debian Linux. On Windows/MinGW, known MinGW-w64 toolchain issues have been fixed (a `MAKEINTRESOURCE` pointer-narrowing error, an unused Direct3D vtable reference pulled in by the OpenGL-only AntTweakBar build, and a GLAD loader failure caused by GLFW2's `glfwGetProcAddress()` lacking a `GetProcAddress()` fallback), using solutions already verified in the standalone [AntTweakBar-Legacy](https://github.com/n-s-kiselev/AntTweakBar-Legacy) project; running the build directly on a Windows machine has not yet been re-verified in this repository.

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

See [`LICENSE`](LICENSE). AntTweakBar and the other vendored components retain their own license files in their respective directories under [`vendor`](vendor).
