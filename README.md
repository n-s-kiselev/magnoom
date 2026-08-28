![Magnoom](https://github.com/n-s-kiselev/magnoom/blob/master/MagnoomWiki/TitleImage.png)

# Magnoom

Magnoom is cross-platform software for atomistic spin-dynamics simulations with real-time OpenGL visualization and interactive parameter control, running natively on Linux, macOS, and Windows. See the [project wiki](https://github.com/n-s-kiselev/magnoom/wiki) for additional documentation and examples.

## Main features

- Cross-platform C99 application code for Linux, macOS, and Windows (MinGW), built with the same single cross-platform build script on every platform
- A few multithreaded solvers for Landau-Lifshitz-Gilbert equation, including finite-temperature simulation via stochastic LLG
- A wide range of pairwise interactions (exchange, Dzyaloshinskii-Moriya, biquadratic) and single-ion anisotropy up to 4th order in general tensor form
- Support for crystal lattices of arbitrary complexity
- Import and export support for formats used by
  [OOMMF](https://math.nist.gov/oommf/) and
  [MuMax3](https://mumax.github.io/), including [OVF](https://math.nist.gov/oommf/doc/userguide20b0/userguide/Vector_Field_File_Format_OV.html)
- OpenGL visualization and post-processing tools, including slicing and filtering
- Real-time parameter control


## Building

Clone Magnoom with its AntTweakBar submodule, then enter the repository:

```sh
git clone --recurse-submodules https://github.com/n-s-kiselev/magnoom.git
cd magnoom
```

If you by have already cloned without `--recurse-submodules`, initialize AntTweakBar from
the repository root:

```sh
git submodule update --init
```

Magnoom uses a cross-platform [nob.h](https://github.com/tsoding/nob.h) build system. 
Bootstrap it once from the repository root:

```sh
cc nob.c -o nob
```
and then 

```sh
./nob
```

You can use `./nob -help` to see the list of supported optionscan

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

Magnoom's application sources are compiled as C99. The vendored AntTweakBar
implementation remains C++, so the final executable is linked with the C++
compiler driver.

The `nob` executable automatically rebuilds itself when `nob.c` or its
vendored `nob.h` changes. All normal build products are contained in `build/`,
except for the bootstrapped `nob` executable itself. `nob.c` is the only build
system for this repository — no CMake, Makefile, Visual Studio project, or
package-manager build is required. The only network access a build may
trigger is the one-time git submodule fetch described above; there is no
package-manager dependency resolution.

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

Developer: [Nikolai S. Kiselev](https://www.researchgate.net/profile/Nikolai-Kiselev-3)

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
