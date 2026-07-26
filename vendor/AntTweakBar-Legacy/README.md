# AntTweakBar (Legacy)

[AntTweakBar](https://anttweakbar.sourceforge.io/doc) (**ATB**) is a small and
easy-to-use C/C++ library developed by
[Philippe Decaudin](https://phildec.users.sourceforge.net/) that adds a
lightweight, intuitive GUI to OpenGL-based graphics programs for real-time
parameter tweaking. Representative screenshots can be found, for example, on
the [magnoom](https://github.com/n-s-kiselev/magnoom) page.

This is a maintained checkout of the legacy AntTweakBar development version,
built with a single cross-platform [`nob.c`](https://github.com/tsoding/nob.h)
build script instead of the original per-platform Makefiles/Visual Studio
project. Only the OpenGL backend is targeted (Direct3D9/10/11 and the OpenGL
Core Profile are not supported by this build).

## Changes from the upstream development version

Compared with the development version published on the
[official AntTweakBar website](https://anttweakbar.sourceforge.io/doc/), this
repository:

- corrects a few small typos inherited from the upstream source and
  documentation;
- replaces the original platform-specific build files with the single
  cross-platform `nob.c` build;
- builds the supported examples from vendored GLFW2, FreeGLUT, and GLAD
  sources;
- fixes the AntTweakBar library's custom cursor creation on modern macOS; and
- includes cross-platform example and integration fixes, including HiDPI
  sizing, resizing, and safe FreeGLUT shutdown.

## GLFW and OpenGL compatibility

The examples in this repository are compiled with the vendored GLFW v2
implementation and use AntTweakBar's legacy OpenGL compatibility-profile
backend.

For AntTweakBar with GLFW v3 and the OpenGL Core Profile (OpenGL 3.0 and
later), see [n-s-kiselev/AntTweakBarGLFW3](https://github.com/n-s-kiselev/AntTweakBarGLFW3).

## macOS custom cursors

AntTweakBar stores its custom 32×32 cursor shapes as data in the library's
source code. The original macOS code used this data to create a compact 2-bit
image, which modern macOS versions display as fully transparent. As a result,
the pointer disappeared whenever AntTweakBar selected a custom cursor. The
problem affected both GLFW and FreeGLUT applications because the faulty image
was created by AntTweakBar itself.

This repository fixes the source code in `CTwMgr::PixmapCursor`. While an
application is running, the function now creates a standard 32-bit RGBA image
in memory from the built-in cursor shape and passes that image to macOS. The
cursor mask supplies the transparent and opaque areas. This is a library code
fix; it does not generate, modify, or post-process any bitmap files during the
build.

## Building

Bootstrap the build tool once, from the repository root:

```sh
gcc nob.c -o nob
```

Then:

```sh
./nob            # build the library (lib/libAntTweakBar.a + the platform dynamic library)
./nob -clean     # remove all generated build output
./nob -examples  # build the example programs (requires ./nob to have run first)
./nob -help      # list all flags
```

`./nob` produces:

- `lib/libAntTweakBar.a` — static library
- `lib/libAntTweakBar.so` (Linux) / `lib/libAntTweakBar.dylib` (macOS) /
  `lib/libAntTweakBar.dll` + `lib/libAntTweakBar.dll.a` (Windows/MinGW) —
  dynamic library

`./nob -examples` compiles the supported examples listed in `nob.c` statically against
`lib/libAntTweakBar.a` into `build/examples/`, and fails with a clear message
if `lib/libAntTweakBar.a` doesn't exist yet (i.e. if `./nob` hasn't been run).
GLFW2 and FreeGLUT are vendored under `vendor/` and built
from source as part of this step, so no external GLFW/GLUT install is
needed on Linux, macOS, or Windows; see
[`examples/Readme.txt`](examples/Readme.txt) for details. The deprecated
macOS system `GLUT.framework` is not used.

Supported platforms: Linux, macOS, and Windows (MinGW).

## License

See [`License.txt`](License.txt).
