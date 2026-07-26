// glfw2_unity.c - single-file compilation of the vendored GLFW 2.7.9 sources
// under lib/, so the TwSimpleGLFW2.c example links against a built-in GLFW2
// instead of requiring one to be installed system-wide. Every GLFW example
// uses this version because it matches the event and native-cursor integration
// AntTweakBar originally supported.
//
// Platform backend selection (must be defined by the build before this file
// is compiled - see append_glfw2_flags() in nob.c):
//   _GLFW2_WIN32   Win32 API (Windows/MinGW)
//   _GLFW2_COCOA   Cocoa frameworks (macOS)
//   _GLFW2_X11     X Window System (Linux)

// lib/init.c must be included first: internal.h's GLFWGLOBAL macro (empty
// vs. "extern") is decided by whether _init_c_ is already defined at
// internal.h's *first* inclusion (its own #include guard makes every later
// #include "internal.h" a no-op) - only init.c defines _init_c_, so it must
// be the file that triggers that first inclusion, or the global variable
// definitions never get emitted anywhere in this translation unit.
#include "lib/init.c"

// Common modules shared by every platform
#include "lib/enable.c"
#include "lib/fullscreen.c"
#include "lib/glext.c"
#include "lib/image.c"
#include "lib/input.c"
#include "lib/joystick.c"
#include "lib/stream.c"
#include "lib/tga.c"
#include "lib/thread.c"
#include "lib/time.c"
#include "lib/window.c"

#if defined(_GLFW2_WIN32)
#   include "lib/win32/win32_enable.c"
#   include "lib/win32/win32_fullscreen.c"
#   include "lib/win32/win32_glext.c"
#   include "lib/win32/win32_init.c"
#   include "lib/win32/win32_joystick.c"
#   include "lib/win32/win32_thread.c"
#   include "lib/win32/win32_time.c"
#   include "lib/win32/win32_window.c"
#endif

#if defined(_GLFW2_X11)
#   include "lib/x11/x11_enable.c"
#   include "lib/x11/x11_fullscreen.c"
#   include "lib/x11/x11_glext.c"
#   include "lib/x11/x11_init.c"
#   include "lib/x11/x11_joystick.c"
#   include "lib/x11/x11_keysym2unicode.c"
#   include "lib/x11/x11_thread.c"
#   include "lib/x11/x11_time.c"
#   include "lib/x11/x11_window.c"
#endif

#if defined(_GLFW2_COCOA)
#   include "lib/cocoa/cocoa_enable.m"
#   include "lib/cocoa/cocoa_fullscreen.m"
#   include "lib/cocoa/cocoa_glext.m"
#   include "lib/cocoa/cocoa_init.m"
#   include "lib/cocoa/cocoa_joystick.m"
#   include "lib/cocoa/cocoa_thread.c"
#   include "lib/cocoa/cocoa_time.m"
#   include "lib/cocoa/cocoa_window.m"
#endif
