# Changing the Magnoom Icon

## 1. Edit the master artwork

Edit `assets/icon/magnoom.svg` in Inkscape. Keep a square 1024x1024 view box,
transparent corners, and details that remain readable at 16x16. Do not edit the
generated PNG, ICO, ICNS, or X11 header manually.

## 2. Bootstrap the exporter

On Linux or macOS:

```sh
cc -std=c99 -O2 nob_icons.c -o nob_icons
```

On Windows with MinGW:

```sh
gcc -std=c99 -O2 nob_icons.c -o nob_icons.exe
```

The exporter requires Inkscape. It uses `INKSCAPE` when set, otherwise the
standard macOS application location when present, and finally `inkscape` from
`PATH`:

```sh
INKSCAPE=/path/to/inkscape ./nob_icons
```

The approved PNG decoder is vendored in `vendor/stb/stb_image.h` and is compiled
only into the exporter. Python, Pillow, and ImageMagick are not required.

## 3. Generate the assets

Run:

```sh
./nob_icons
```

The exporter creates 16, 32, 48, 64, 128, 256, 512, and 1024 pixel PNGs in
`build/icons/`, validates their dimensions and transparency, then generates:

- `assets/icon/magnoom-256.png`
- `assets/icon/magnoom-512.png`
- `assets/icon/magnoom.ico`
- `assets/icon/magnoom_x11_icon.h`
- `assets/icon/magnoom.icns` on macOS

ICO data is assembled directly from Inkscape's PNG files. The X11 header is
generated from the decoded 64x64 RGBA image. On macOS, the exporter calls the
system `iconutil` command. On Linux and Windows it preserves the existing ICNS
file because `iconutil` is macOS-only.

Available commands:

```sh
./nob_icons          # regenerate only when inputs are newer
./nob_icons -force   # regenerate every icon
./nob_icons -clean   # remove build/icons temporary files
./nob_icons -help    # show usage
```

## 4. Rebuild Magnoom

The regular build detects the regenerated assets:

```sh
./nob
```

For complete validation:

```sh
./nob -clean
./nob
```

## 5. Check each platform

On macOS:

```sh
plutil -lint build/Magnoom.app/Contents/Info.plist
codesign --verify --deep --strict build/Magnoom.app
open build/Magnoom.app
```

On Windows, inspect `build/magnoom.exe`, its window, and its taskbar icon.

On Linux, inspect the running window and the launcher assets under
`build/share/icons/hicolor/` and `build/share/applications/`.
