# Magnoom Icon Assets

`magnoom.svg` is the editable master artwork. The PNG, ICO, ICNS, and X11
header files are platform-specific renderings of that violet-and-gold source.
Regenerate them with the standalone `nob_icons.c` exporter described in
`icon_instructions.md`; Python and ImageMagick are not required.

- `magnoom.ico` is embedded in the Windows executable through `magnoom.rc`.
- `magnoom.icns` is copied into the macOS application bundle.
- `magnoom-256.png` and `magnoom-512.png` are staged for Linux launchers.
- `magnoom_x11_icon.h` contains the 64x64 ARGB pixels published through
  `_NET_WM_ICON` by the vendored GLFW2 X11 backend.
