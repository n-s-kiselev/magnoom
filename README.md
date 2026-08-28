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

`magnoom.cfg` (repo root, next to the built executable) configures Magnoom at
startup. It uses plain `key = value` lines; `#` starts a comment that runs to
the end of the line, blank lines are ignored, and whitespace around keys,
`=`, and values is ignored. Unknown keys and malformed or out-of-range values
are rejected with a `magnoom.cfg:<line>: ...` error and the program exits — a
typo in the file is never silently ignored. If `magnoom.cfg` is missing, a
default file is written in its place and the program starts from built-in
defaults; if it can't be written (e.g. a read-only directory), the program
just starts from the built-in defaults. Values not mentioned in the file keep
their default. The old `# begin magnoom config` / `# key: value` syntax is no
longer supported; a file starting with the old header produces a clear error
telling you to convert it.

`magnoom.cfg.example` (repo root) is a reference listing of every key
`readConfigFile()` understands — lattice vectors, cell counts, shell count,
boundary conditions, the default atom position, and the full set of per-atom
anisotropy keys — each set to the same value `magnoom_ctx_init()` already
uses, so loading it changes nothing. Copy it to `magnoom.cfg` and edit the
values you actually want to change; any key you omit keeps its default.

```text
# Lattice vectors
ax = 1.0
ay = 0.0
az = 0.0
bx = 0.0
by = 1.0
bz = 0.0
cx = 0.0
cy = 0.0
cz = 1.0

# Number of cells
Na = 50
Nb = 50
Nc = 50

Shells = 4

# Boundary conditions (optional; defaults shown)
BCa = 1
BCb = 1
BCc = 0
```

**Atom positions.** The crystal basis (the atoms inside one unit cell,
`magnoom_ctx_set_block()`'s target) is configured with `Atom<N>x`/`Atom<N>y`/
`Atom<N>z` keys, `N` starting at `0` and counting up with no gaps — the number
of atoms is inferred from how many are given, up to `MAX_ATOMS_PER_BLOCK = 5`.
If no `AtomNx/y/z` keys are present, the built-in one-atom simple-cubic basis
(`{0.5, 0.5, 0.5}`) is kept. Positions are given in the same Cartesian units
as the lattice vectors above, and must map to a fractional position inside
the half-open unit cell `[0, 1)` once expressed in that lattice basis. For
example, a two-atom basis:

```text
Atom0x = 0.0
Atom0y = 0.0
Atom0z = 0.0

Atom1x = 0.5
Atom1y = 0.5
Atom1z = 0.5
```

A few other crystal structures, for reference (set the matching lattice
vectors above alongside the atom positions):

```text
# B20 (u = 0.138), cubic unit cell (identity lattice vectors, as shown above):
Atom0x = 0.0     Atom0y = 0.0     Atom0z = 0.0
Atom1x = 0.5     Atom1y = 0.224   Atom1z = 0.724
Atom2x = 0.724   Atom2y = 0.5     Atom2z = 0.224
Atom3x = 0.224   Atom3y = 0.724   Atom3z = 0.5

# FCC2, orthogonal unit cell (Cartesian positions; by=1.4142136, cz=1.4142136):
Atom0x = 0.0     Atom0y = 0.0        Atom0z = 0.0
Atom1x = 0.0     Atom1y = 0.7071068  Atom1z = 0.0
Atom2x = 0.5     Atom2y = 0.3535534  Atom2z = 0.3535534
Atom3x = 0.5     Atom3y = 1.0606602  Atom3z = 0.3535534
Atom4x = 0.0     Atom4y = 0.7071068  Atom4z = 0.7071068
```

**Anisotropy.** Per-atom magnetocrystalline anisotropy tensor and rotation
components are set with `Atom<N>_...` keys (or `AtomAll_...` to assign every
atom at once, matching the atom `-1` wildcard used internally):

```text
# Rank-2 and rank-4 components use one-based tensor indices 1..3
Atom0_K11 = 0.5
Atom0_K1111 = 0.1

# Local-to-global rotation quaternion
Atom0_Qx = 0
Atom0_Qy = 0
Atom0_Qz = 0
Atom0_Qs = 1

# AtomAll applies to every atom in the active basis
AtomAll_K22 = 0.2
```

The `K11`..`K33` (6) and `K1111`..`K3333` (15) component names match the F6
Anisotropy bar's control labels exactly. Any subset may be given — omitted
components stay at their default (normally zero); this is why the syntax is
sparse. The four quaternion components (`Qx`, `Qy`, `Qz`, `Qs`, AntTweakBar's
`{x, y, z, scalar}` order) must be given together if any one of them is; they
are normalized automatically and must have a nonzero finite norm. Keys are
applied in file order, so a later key overrides an earlier one for the same
atom and component.

Press `F6` to open the Anisotropy bar. It always exposes all six independent K2
and fifteen independent K4 components, including components initialized to
zero. Global mode applies atom 0's tensor to every spin; Individual mode
provides an atom selector and uses each basis atom's tensor. The copy button
copies atom 0's local K2/K4 tensors to every atom without changing their
quaternions. The Rotation widget exposes all four `{qx, qy, qz, qs}` components
and a 3D rotation ball; editing any component or dragging the ball immediately
refreshes the active atom's rotated tensor. The Axis and Angle fields define an
additional rotation; pressing **compose axis-angle** multiplies the corresponding
quaternion by the current one so rotations can be built up one after another.
**reset to identity** restores the identity quaternion `{0, 0, 0, 1}`.

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
