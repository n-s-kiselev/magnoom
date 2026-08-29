#define NOB_IMPLEMENTATION
#include "vendor/nob/nob.h"

#define ATB_ROOT "vendor/AntTweakBar-Legacy/"
#define ATB_SRC ATB_ROOT "src/"
#define ATB_INCLUDE ATB_ROOT "include/"
#define BUILD_DIR "build/"

#if defined(_WIN32)
#define EXE_EXT ".exe"
#define PLATFORM_SUFFIX "-windows"
#elif defined(__APPLE__)
#define EXE_EXT ""
#define PLATFORM_SUFFIX "-macos"
#else
#define EXE_EXT ""
#define PLATFORM_SUFFIX "-linux"
#endif

// Every output below is named per-platform (PLATFORM_SUFFIX) so that this
// repository can live in a shared, Dropbox-synced folder and be built on
// Linux, macOS, and Windows in turn without one platform's build overwriting
// or being cleaned out by another's.
#define ATB_BUILD_DIR BUILD_DIR "anttweakbar" PLATFORM_SUFFIX "/"
#define ATB_LIB BUILD_DIR "libAntTweakBar" PLATFORM_SUFFIX ".a"
#define GLAD_INCLUDE "vendor/glad/include/"
#define GLAD_SOURCE "vendor/glad/src/glad.c"
#define GLAD_OBJECT BUILD_DIR "glad" PLATFORM_SUFFIX ".o"
#define GLFW_INCLUDE "vendor/glfw2/include/"
#define GLFW_SOURCE "vendor/glfw2/glfw2_unity.c"
#define GLFW_OBJECT BUILD_DIR "glfw2" PLATFORM_SUFFIX ".o"
#define MAGNOOM_OBJECT BUILD_DIR "magnoom" PLATFORM_SUFFIX ".o"
#define BLOCK_SETTER_TEST_SOURCE "tests/block_setter_test.c"
#define BLOCK_SETTER_TEST_OBJECT BUILD_DIR "block_setter_test" PLATFORM_SUFFIX ".o"
#define BLOCK_SETTER_TEST_OUTPUT BUILD_DIR "block_setter_test" PLATFORM_SUFFIX
#define OUTPUT BUILD_DIR "magnoom" PLATFORM_SUFFIX
#define ICON_DIR "assets/icon/"
// Windows- and macOS-only outputs are not suffixed: each is only ever
// produced under its own platform's #ifdef branch, so there is no
// cross-platform collision to guard against.
#define WINDOWS_RESOURCE_OBJECT BUILD_DIR "magnoom-resource.o"
#define MACOS_APP BUILD_DIR "Magnoom.app/"
#define NOB_HEADER "vendor/nob/nob.h"

static const char *atb_common_sources[] = {
    ATB_SRC "TwColors.cpp", ATB_SRC "TwFonts.cpp", ATB_SRC "TwOpenGL.cpp",
    ATB_SRC "TwOpenGLCore.cpp", ATB_SRC "TwBar.cpp", ATB_SRC "TwMgr.cpp",
    ATB_SRC "LoadOGL.cpp", ATB_SRC "LoadOGLCore.cpp",
    ATB_SRC "TwEventGLFW.c",
};

static bool build_needed(const char *output, const char **inputs, size_t count)
{
    int result = nob_needs_rebuild(output, inputs, count);
    if (result < 0) exit(1);
    return result > 0;
}

static const char *object_path(const char *source)
{
    char *name = nob_temp_strdup(nob_path_name(source));
    char *dot = strrchr(name, '.');
    if (dot) *dot = '\0';
    return nob_temp_sprintf(ATB_BUILD_DIR "%s.o", name);
}

static bool build_atb_object(const char *source)
{
    const char *output = object_path(source);
    const char *inputs[] = { source, "nob.c", NOB_HEADER };
    if (!build_needed(output, inputs, NOB_ARRAY_LEN(inputs))) return true;

    Nob_Cmd cmd = {0};
#if defined(__APPLE__)
    nob_cmd_append(&cmd, "c++", "-x", "objective-c++", "-D_MACOSX", "-DGL_SILENCE_DEPRECATION",
                   "-Wno-deprecated-declarations");
#elif defined(_WIN32)
    nob_cmd_append(&cmd, nob_sv_ends_with_cstr(nob_sv_from_cstr(source), ".cpp") ? "c++" : "cc", "-D_WINDOWS");
    // TwMgr.cpp's CreateCursors() calls MAKEINTRESOURCE() on values like
    // IDC_ARROW, which MinGW-w64's own <winuser.h> already defines via
    // MAKEINTRESOURCE (i.e. as a pointer, not a raw integer id) - so this is
    // MAKEINTRESOURCE applied twice. That's always been a redundant no-op at
    // runtime (the "fake pointer" round-trips through the WORD truncation
    // unchanged), but recent MinGW-w64 GCC (confirmed with 14.2.0) now
    // treats the implicit pointer-to-WORD narrowing inside that expansion as
    // a hard error in C++ instead of a warning. -fpermissive downgrades it
    // back to a warning without touching this stock upstream AntTweakBar
    // source (same fix already applied and confirmed working in the
    // standalone AntTweakBar-Legacy project this is vendored from).
    nob_cmd_append(&cmd, "-fpermissive");
    // This build targets the OpenGL backend only and does not compile/link
    // TwDirect3D9/10/11.cpp (see atb_common_sources above) - without this,
    // TwMgr.cpp's TwCreateGraph() unconditionally does "new
    // CTwGraphDirect3D9/10/11" under plain ANT_WINDOWS, requiring those
    // classes' vtables at link time even though magnoom never requests
    // TW_DIRECT3D9/10/11.
    nob_cmd_append(&cmd, "-DTW_NO_DIRECT3D");
#else
    nob_cmd_append(&cmd, nob_sv_ends_with_cstr(nob_sv_from_cstr(source), ".cpp") ? "c++" : "cc", "-D_UNIX");
#endif
    nob_cmd_append(&cmd, "-Wall", "-O2", "-fno-strict-aliasing", "-DTW_STATIC",
                   "-I" ATB_INCLUDE, "-c", source, "-o", output);
    return nob_cmd_run(&cmd);
}

static bool build_ant_tweak_bar(void)
{
    Nob_File_Paths objects = {0};
    for (size_t i = 0; i < NOB_ARRAY_LEN(atb_common_sources); ++i) {
        if (!build_atb_object(atb_common_sources[i])) return false;
        nob_da_append(&objects, object_path(atb_common_sources[i]));
    }
    if (!build_needed(ATB_LIB, objects.items, objects.count)) return true;
    if (nob_file_exists(ATB_LIB) && !nob_delete_file(ATB_LIB)) return false;
    Nob_Cmd cmd = {0};
    nob_cmd_append(&cmd, "ar", "rcs", ATB_LIB);
    for (size_t i = 0; i < objects.count; ++i) nob_cmd_append(&cmd, objects.items[i]);
    return nob_cmd_run(&cmd);
}

static bool build_glad(void)
{
    const char *inputs[] = { GLAD_SOURCE, "nob.c", NOB_HEADER };
    if (!build_needed(GLAD_OBJECT, inputs, NOB_ARRAY_LEN(inputs))) return true;
    Nob_Cmd cmd = {0};
    nob_cmd_append(&cmd, "cc", "-O2", "-I" GLAD_INCLUDE, "-c", GLAD_SOURCE, "-o", GLAD_OBJECT);
    return nob_cmd_run(&cmd);
}

static bool build_glfw(void)
{
    const char *inputs[] = {
        GLFW_SOURCE, "vendor/glfw2/lib/x11/x11_window.c",
        ICON_DIR "magnoom_x11_icon.h", "nob.c", NOB_HEADER,
    };
    if (!build_needed(GLFW_OBJECT, inputs, NOB_ARRAY_LEN(inputs))) return true;
    Nob_Cmd cmd = {0};
    nob_cmd_append(&cmd, "cc", "-O2", "-I" GLFW_INCLUDE, "-Ivendor/glfw2/lib",
                   "-I" ICON_DIR, "-DGLFW_NO_GLU");
#if defined(_WIN32)
    nob_cmd_append(&cmd, "-D_GLFW2_WIN32", "-Ivendor/glfw2/lib/win32");
#elif defined(__APPLE__)
    nob_cmd_append(&cmd, "-D_GLFW2_COCOA", "-Ivendor/glfw2/lib/cocoa", "-x", "objective-c",
                   "-Wno-deprecated-declarations");
#else
    nob_cmd_append(&cmd, "-D_GLFW2_X11", "-Ivendor/glfw2/lib/x11",
                   "-D_GLFW_HAS_XRANDR", "-D_GLFW_HAS_PTHREAD", "-D_GLFW_HAS_SCHED_YIELD",
                   "-D_GLFW_HAS_GLXGETPROCADDRESS", "-D_GLFW_HAS_SYSCONF", "-D_GLFW_USE_LINUX_JOYSTICKS");
#endif
    nob_cmd_append(&cmd, "-c", GLFW_SOURCE, "-o", GLFW_OBJECT);
    return nob_cmd_run(&cmd);
}

static bool build_windows_resource(void)
{
#if defined(_WIN32)
    const char *inputs[] = {
        ICON_DIR "magnoom.rc", ICON_DIR "magnoom.ico", "nob.c", NOB_HEADER,
    };
    if (!build_needed(WINDOWS_RESOURCE_OBJECT, inputs, NOB_ARRAY_LEN(inputs))) return true;

    Nob_Cmd cmd = {0};
    nob_cmd_append(&cmd, "windres", ICON_DIR "magnoom.rc", WINDOWS_RESOURCE_OBJECT);
    return nob_cmd_run(&cmd);
#else
    return true;
#endif
}

static bool build_magnoom_object(void)
{
    const char *inputs[] = {
        "magnoom.c", "anisotropy.c", "solvers.c", "lattice_geometry.c", "initial_states.c",
        "math_utils.c", "visualization.c", "linmath.h",
        ATB_LIB, GLAD_OBJECT, GLFW_OBJECT,
        "vendor/stb/stb_image_write.h",
        "vendor/glfw2/TwGLFW2.h", "nob.c", NOB_HEADER,
    };
    if (!build_needed(MAGNOOM_OBJECT, inputs, NOB_ARRAY_LEN(inputs))) return true;

    Nob_Cmd cmd = {0};
    nob_cmd_append(&cmd, "cc", "-std=c99", "-O3", "-Wall", "-fno-strict-aliasing", "-DTW_STATIC",
                   "-I" ATB_INCLUDE, "-I" GLAD_INCLUDE, "-I" GLFW_INCLUDE,
                   "-c", "magnoom.c", "-o", MAGNOOM_OBJECT);
    return nob_cmd_run(&cmd);
}

static bool build_block_setter_test_object(void)
{
    const char *inputs[] = {
        BLOCK_SETTER_TEST_SOURCE, "magnoom.c", "anisotropy.c", "solvers.c", "lattice_geometry.c",
        "initial_states.c", "math_utils.c", "visualization.c", "linmath.h",
        ATB_LIB, GLAD_OBJECT, GLFW_OBJECT, "vendor/stb/stb_image_write.h",
        "vendor/glfw2/TwGLFW2.h", "nob.c", NOB_HEADER,
    };
    if (!build_needed(BLOCK_SETTER_TEST_OBJECT, inputs, NOB_ARRAY_LEN(inputs))) return true;

    Nob_Cmd cmd = {0};
    nob_cmd_append(&cmd, "cc", "-std=c99", "-O2", "-Wall", "-fno-strict-aliasing", "-DTW_STATIC",
                   "-I" ATB_INCLUDE, "-I" GLAD_INCLUDE, "-I" GLFW_INCLUDE,
                   "-c", BLOCK_SETTER_TEST_SOURCE, "-o", BLOCK_SETTER_TEST_OBJECT);
    return nob_cmd_run(&cmd);
}

static bool build_magnoom(void)
{
#if defined(_WIN32)
    const char *inputs[] = {
        MAGNOOM_OBJECT, ATB_LIB, GLAD_OBJECT, GLFW_OBJECT, WINDOWS_RESOURCE_OBJECT,
    };
#else
    const char *inputs[] = { MAGNOOM_OBJECT, ATB_LIB, GLAD_OBJECT, GLFW_OBJECT };
#endif
    const char *output = OUTPUT EXE_EXT;
    if (!build_needed(output, inputs, NOB_ARRAY_LEN(inputs))) return true;

    Nob_Cmd cmd = {0};
    nob_cmd_append(&cmd, "c++", MAGNOOM_OBJECT, GLAD_OBJECT, ATB_LIB, GLFW_OBJECT);
#if defined(_WIN32)
    nob_cmd_append(&cmd, WINDOWS_RESOURCE_OBJECT);
#endif
    nob_cmd_append(&cmd, "-o", output);
#if defined(__APPLE__)
    nob_cmd_append(&cmd, "-framework", "OpenGL", "-framework", "Cocoa",
                   "-framework", "AppKit", "-framework", "Foundation", "-framework", "IOKit",
                   "-framework", "CoreVideo", "-pthread",
                   "-Wno-deprecated-declarations", "-lobjc");
#elif defined(_WIN32)
    // Statically link the MinGW runtime (libstdc++, libgcc, and libwinpthread if pulled in)
    // so the .exe doesn't depend on DLLs that may not be present on a machine without the
    // MinGW toolchain installed (this links via c++, which dynamically links libstdc++ by default).
    nob_cmd_append(&cmd, "-static", "-lopengl32", "-lgdi32", "-lwinmm");
#else
    nob_cmd_append(&cmd, "-lGL", "-lX11", "-lXrandr", "-lpthread", "-ldl", "-lm");
#endif
    return nob_cmd_run(&cmd);
}

static bool build_block_setter_test(void)
{
    const char *inputs[] = { BLOCK_SETTER_TEST_OBJECT, ATB_LIB, GLAD_OBJECT, GLFW_OBJECT };
    const char *output = BLOCK_SETTER_TEST_OUTPUT EXE_EXT;
    if (!build_needed(output, inputs, NOB_ARRAY_LEN(inputs))) return true;

    Nob_Cmd cmd = {0};
    nob_cmd_append(&cmd, "c++", BLOCK_SETTER_TEST_OBJECT, GLAD_OBJECT, ATB_LIB, GLFW_OBJECT,
                   "-o", output);
#if defined(__APPLE__)
    nob_cmd_append(&cmd, "-framework", "OpenGL", "-framework", "Cocoa",
                   "-framework", "AppKit", "-framework", "Foundation", "-framework", "IOKit",
                   "-framework", "CoreVideo", "-pthread",
                   "-Wno-deprecated-declarations", "-lobjc");
#elif defined(_WIN32)
    nob_cmd_append(&cmd, "-static", "-lopengl32", "-lgdi32", "-lwinmm");
#else
    nob_cmd_append(&cmd, "-lGL", "-lX11", "-lXrandr", "-lpthread", "-ldl", "-lm");
#endif
    return nob_cmd_run(&cmd);
}

static bool run_block_setter_test(void)
{
    Nob_Cmd cmd = {0};
    nob_cmd_append(&cmd, BLOCK_SETTER_TEST_OUTPUT EXE_EXT);
    return nob_cmd_run(&cmd);
}

static bool copy_if_needed(const char *source, const char *destination, bool *copied)
{
    const char *inputs[] = { source };
    if (!build_needed(destination, inputs, NOB_ARRAY_LEN(inputs))) return true;
    if (!nob_copy_file(source, destination)) return false;
    if (copied) *copied = true;
    return true;
}

static bool package_platform_assets(void)
{
#if defined(__APPLE__)
    bool bundle_changed = false;
    if (!nob_mkdir_if_not_exists(MACOS_APP) ||
        !nob_mkdir_if_not_exists(MACOS_APP "Contents/") ||
        !nob_mkdir_if_not_exists(MACOS_APP "Contents/MacOS/") ||
        !nob_mkdir_if_not_exists(MACOS_APP "Contents/Resources/")) return false;

    if (!copy_if_needed(OUTPUT, MACOS_APP "Contents/MacOS/magnoom", &bundle_changed) ||
        !copy_if_needed("assets/macos/Info.plist", MACOS_APP "Contents/Info.plist",
                        &bundle_changed) ||
        !copy_if_needed(ICON_DIR "magnoom.icns",
                        MACOS_APP "Contents/Resources/magnoom.icns",
                        &bundle_changed)) return false;

    if (!bundle_changed &&
        nob_file_exists(MACOS_APP "Contents/_CodeSignature/CodeResources")) return true;

    Nob_Cmd cmd = {0};
    nob_cmd_append(&cmd, "codesign", "--force", "--deep", "--sign", "-", MACOS_APP);
    return nob_cmd_run(&cmd);
#elif defined(_WIN32)
    return true;
#else
    if (!nob_mkdir_if_not_exists(BUILD_DIR "share/") ||
        !nob_mkdir_if_not_exists(BUILD_DIR "share/applications/") ||
        !nob_mkdir_if_not_exists(BUILD_DIR "share/icons/") ||
        !nob_mkdir_if_not_exists(BUILD_DIR "share/icons/hicolor/") ||
        !nob_mkdir_if_not_exists(BUILD_DIR "share/icons/hicolor/256x256/") ||
        !nob_mkdir_if_not_exists(BUILD_DIR "share/icons/hicolor/256x256/apps/") ||
        !nob_mkdir_if_not_exists(BUILD_DIR "share/icons/hicolor/512x512/") ||
        !nob_mkdir_if_not_exists(BUILD_DIR "share/icons/hicolor/512x512/apps/")) return false;

    return copy_if_needed("assets/linux/magnoom.desktop",
                          BUILD_DIR "share/applications/magnoom.desktop", NULL) &&
           copy_if_needed(ICON_DIR "magnoom-256.png",
                          BUILD_DIR "share/icons/hicolor/256x256/apps/magnoom.png", NULL) &&
           copy_if_needed(ICON_DIR "magnoom-512.png",
                          BUILD_DIR "share/icons/hicolor/512x512/apps/magnoom.png", NULL);
#endif
}

static bool clear_directory(const char *folder)
{
    if (!nob_file_exists(folder)) return true;
    Nob_File_Paths children = {0};
    if (!nob_read_entire_dir(folder, &children)) return false;
    for (size_t i = 0; i < children.count; ++i) {
        if (strcmp(children.items[i], ".") == 0 || strcmp(children.items[i], "..") == 0) continue;
        const char *path = nob_temp_sprintf("%s%s", folder, children.items[i]);
        if (nob_get_file_type(path) == NOB_FILE_DIRECTORY) {
            if (!clear_directory(nob_temp_sprintf("%s/", path))) return false;
        } else if (!nob_delete_file(path)) {
            return false;
        }
    }
    return nob_delete_file(folder);
}

static bool delete_if_exists(const char *path)
{
    return !nob_file_exists(path) || nob_delete_file(path);
}

static bool clean(void)
{
    bool ok = true;
    ok = clear_directory(ATB_BUILD_DIR) && ok;
    ok = delete_if_exists(ATB_LIB) && ok;
    ok = delete_if_exists(GLAD_OBJECT) && ok;
    ok = delete_if_exists(GLFW_OBJECT) && ok;
    ok = delete_if_exists(MAGNOOM_OBJECT) && ok;
    ok = delete_if_exists(BLOCK_SETTER_TEST_OBJECT) && ok;
    ok = delete_if_exists(BLOCK_SETTER_TEST_OUTPUT EXE_EXT) && ok;
    ok = delete_if_exists(OUTPUT EXE_EXT) && ok;
#if defined(_WIN32)
    ok = delete_if_exists(WINDOWS_RESOURCE_OBJECT) && ok;
#elif defined(__APPLE__)
    ok = clear_directory(MACOS_APP) && ok;
#else
    ok = clear_directory(BUILD_DIR "share/") && ok;
#endif

    // One-time sweep of un-suffixed outputs from before per-platform build
    // isolation was introduced, so stale orphans left by an old ./nob build
    // don't linger forever in the shared build/ directory.
    ok = clear_directory(BUILD_DIR "anttweakbar/") && ok;
    ok = delete_if_exists(BUILD_DIR "libAntTweakBar.a") && ok;
    ok = delete_if_exists(BUILD_DIR "glad.o") && ok;
    ok = delete_if_exists(BUILD_DIR "glfw2.o") && ok;
    ok = delete_if_exists(BUILD_DIR "magnoom.o") && ok;
    ok = delete_if_exists(BUILD_DIR "block_setter_test.o") && ok;
    ok = delete_if_exists(BUILD_DIR "block_setter_test") && ok;
    ok = delete_if_exists(BUILD_DIR "magnoom") && ok;

    return ok;
}

static bool ensure_submodules(void)
{
    const char *marker = ATB_SRC "TwMgr.cpp";
    if (nob_file_exists(marker)) return true;

    if (!nob_file_exists(".git")) {
        nob_log(NOB_ERROR, "%s is missing and this is not a git checkout (no .git found) - "
                "fetch the vendor/AntTweakBar-Legacy submodule manually", marker);
        return false;
    }

    nob_log(NOB_INFO, "vendor/AntTweakBar-Legacy submodule is not initialized, running `git submodule update --init`");
    Nob_Cmd cmd = {0};
    nob_cmd_append(&cmd, "git", "submodule", "update", "--init", ATB_ROOT);
    if (!nob_cmd_run(&cmd)) return false;

    if (!nob_file_exists(marker)) {
        nob_log(NOB_ERROR, "%s is still missing after `git submodule update --init`", marker);
        return false;
    }
    return true;
}

static void usage(const char *program)
{
    printf("usage: %s [-clean] [-test] [-help]\n", program);
    printf("  -clean  remove generated build files and exit\n");
    printf("  -test   build and run automated tests\n");
    printf("  -help   print this help and exit\n");
}

int main(int argc, char **argv)
{
    NOB_GO_REBUILD_URSELF_PLUS(argc, argv, NOB_HEADER);
    bool run_tests = false;
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-clean") == 0) return clean() ? 0 : 1;
        if (strcmp(argv[i], "-test") == 0) {
            run_tests = true;
            continue;
        }
        if (strcmp(argv[i], "-help") == 0 || strcmp(argv[i], "--help") == 0) {
            usage(argv[0]);
            return 0;
        }
        nob_log(NOB_ERROR, "unknown argument: %s", argv[i]);
        usage(argv[0]);
        return 1;
    }

    if (!ensure_submodules()) return 1;
    if (!nob_mkdir_if_not_exists(BUILD_DIR) || !nob_mkdir_if_not_exists(ATB_BUILD_DIR)) return 1;
    if (run_tests) {
        if (!build_ant_tweak_bar() || !build_glad() || !build_glfw() ||
            !build_block_setter_test_object() || !build_block_setter_test()) return 1;
        return run_block_setter_test() ? 0 : 1;
    }
    if (!build_ant_tweak_bar() || !build_glad() || !build_glfw() ||
        !build_windows_resource() || !build_magnoom_object() ||
        !build_magnoom() || !package_platform_assets()) return 1;
    nob_log(NOB_INFO, "built %s", OUTPUT EXE_EXT);
    return 0;
}
