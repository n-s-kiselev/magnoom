#define NOB_IMPLEMENTATION
#include "vendor/nob/nob.h"

#define ATB_ROOT "vendor/AntTweakBar-Legacy/"
#define ATB_SRC ATB_ROOT "src/"
#define ATB_INCLUDE ATB_ROOT "include/"
#define BUILD_DIR "build/"
#define ATB_BUILD_DIR BUILD_DIR "anttweakbar/"
#define ATB_LIB BUILD_DIR "libAntTweakBar.a"
#define GLAD_INCLUDE "vendor/glad/include/"
#define GLAD_SOURCE "vendor/glad/src/glad.c"
#define GLAD_OBJECT BUILD_DIR "glad.o"
#define GLFW_INCLUDE "vendor/glfw2/include/"
#define GLFW_SOURCE "vendor/glfw2/glfw2_unity.c"
#define GLFW_OBJECT BUILD_DIR "glfw2.o"
#define MAGNOOM_OBJECT BUILD_DIR "magnoom.o"
#define OUTPUT BUILD_DIR "magnoom"
#define NOB_HEADER "vendor/nob/nob.h"

#if defined(_WIN32)
#define EXE_EXT ".exe"
#else
#define EXE_EXT ""
#endif

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
    nob_cmd_append(&cmd, "c++", "-x", "objective-c++", "-D_MACOSX");
#elif defined(_WIN32)
    nob_cmd_append(&cmd, nob_sv_ends_with_cstr(nob_sv_from_cstr(source), ".cpp") ? "c++" : "cc", "-D_WINDOWS");
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
    const char *inputs[] = { GLFW_SOURCE, "nob.c", NOB_HEADER };
    if (!build_needed(GLFW_OBJECT, inputs, NOB_ARRAY_LEN(inputs))) return true;
    Nob_Cmd cmd = {0};
    nob_cmd_append(&cmd, "cc", "-O2", "-I" GLFW_INCLUDE, "-Ivendor/glfw2/lib", "-DGLFW_NO_GLU");
#if defined(_WIN32)
    nob_cmd_append(&cmd, "-D_GLFW2_WIN32", "-Ivendor/glfw2/lib/win32");
#elif defined(__APPLE__)
    nob_cmd_append(&cmd, "-D_GLFW2_COCOA", "-Ivendor/glfw2/lib/cocoa", "-x", "objective-c");
#else
    nob_cmd_append(&cmd, "-D_GLFW2_X11", "-Ivendor/glfw2/lib/x11",
                   "-D_GLFW_HAS_XRANDR", "-D_GLFW_HAS_PTHREAD", "-D_GLFW_HAS_SCHED_YIELD",
                   "-D_GLFW_HAS_GLXGETPROCADDRESS", "-D_GLFW_HAS_SYSCONF", "-D_GLFW_USE_LINUX_JOYSTICKS");
#endif
    nob_cmd_append(&cmd, "-c", GLFW_SOURCE, "-o", GLFW_OBJECT);
    return nob_cmd_run(&cmd);
}

static bool build_magnoom_object(void)
{
    const char *inputs[] = {
        "main.c", "ENGINE.c", "GEOM.c", "INITSTATE.c", "MATH.c", "OPGL.c",
        "global.h", "linmath.h", ATB_LIB, GLAD_OBJECT, GLFW_OBJECT,
        "vendor/stb/stb_image_write.h",
        "vendor/glfw2/TwGLFW2.h", "nob.c", NOB_HEADER,
    };
    if (!build_needed(MAGNOOM_OBJECT, inputs, NOB_ARRAY_LEN(inputs))) return true;

    Nob_Cmd cmd = {0};
    nob_cmd_append(&cmd, "cc", "-std=c99", "-O3", "-Wall", "-fno-strict-aliasing", "-DTW_STATIC",
                   "-I" ATB_INCLUDE, "-I" GLAD_INCLUDE, "-I" GLFW_INCLUDE,
                   "-c", "main.c", "-o", MAGNOOM_OBJECT);
    return nob_cmd_run(&cmd);
}

static bool build_magnoom(void)
{
    const char *inputs[] = { MAGNOOM_OBJECT, ATB_LIB, GLAD_OBJECT, GLFW_OBJECT };
    const char *output = OUTPUT EXE_EXT;
    if (!build_needed(output, inputs, NOB_ARRAY_LEN(inputs))) return true;

    Nob_Cmd cmd = {0};
    nob_cmd_append(&cmd, "c++", MAGNOOM_OBJECT, GLAD_OBJECT, ATB_LIB, GLFW_OBJECT, "-o", output);
#if defined(__APPLE__)
    nob_cmd_append(&cmd, "-framework", "OpenGL", "-framework", "Cocoa",
                   "-framework", "AppKit", "-framework", "Foundation", "-framework", "IOKit",
                   "-framework", "CoreVideo", "-pthread",
                   "-Wno-deprecated-declarations", "-lobjc");
#elif defined(_WIN32)
    nob_cmd_append(&cmd, "-lopengl32", "-lgdi32", "-lwinmm");
#else
    nob_cmd_append(&cmd, "-lGL", "-lX11", "-lXrandr", "-lpthread", "-ldl", "-lm");
#endif
    return nob_cmd_run(&cmd);
}

static bool clear_directory(const char *folder)
{
    if (!nob_file_exists(folder)) return true;
    Nob_File_Paths children = {0};
    if (!nob_read_entire_dir(folder, &children)) return false;
    for (size_t i = 0; i < children.count; ++i) {
        if (strcmp(children.items[i], ".") == 0 || strcmp(children.items[i], "..") == 0) continue;
        if (!nob_delete_file(nob_temp_sprintf("%s%s", folder, children.items[i]))) return false;
    }
    return nob_delete_file(folder);
}

static bool clean(void)
{
    bool ok = true;
    ok = clear_directory(ATB_BUILD_DIR) && ok;
    if (nob_file_exists(ATB_LIB)) ok = nob_delete_file(ATB_LIB) && ok;
    if (nob_file_exists(GLAD_OBJECT)) ok = nob_delete_file(GLAD_OBJECT) && ok;
    if (nob_file_exists(GLFW_OBJECT)) ok = nob_delete_file(GLFW_OBJECT) && ok;
    if (nob_file_exists(MAGNOOM_OBJECT)) ok = nob_delete_file(MAGNOOM_OBJECT) && ok;
    if (nob_file_exists(OUTPUT EXE_EXT)) ok = nob_delete_file(OUTPUT EXE_EXT) && ok;
    if (nob_file_exists(BUILD_DIR)) ok = nob_delete_file(BUILD_DIR) && ok;
    return ok;
}

static void usage(const char *program)
{
    printf("usage: %s [-clean] [-help]\n", program);
    printf("  -clean  remove generated build files and exit\n");
    printf("  -help   print this help and exit\n");
}

int main(int argc, char **argv)
{
    NOB_GO_REBUILD_URSELF_PLUS(argc, argv, NOB_HEADER);
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-clean") == 0) return clean() ? 0 : 1;
        if (strcmp(argv[i], "-help") == 0 || strcmp(argv[i], "--help") == 0) {
            usage(argv[0]);
            return 0;
        }
        nob_log(NOB_ERROR, "unknown argument: %s", argv[i]);
        usage(argv[0]);
        return 1;
    }

    if (!nob_mkdir_if_not_exists(BUILD_DIR) || !nob_mkdir_if_not_exists(ATB_BUILD_DIR)) return 1;
    if (!build_ant_tweak_bar() || !build_glad() || !build_glfw() ||
        !build_magnoom_object() || !build_magnoom()) return 1;
    nob_log(NOB_INFO, "built %s", OUTPUT EXE_EXT);
    return 0;
}
