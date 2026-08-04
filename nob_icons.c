#define NOB_IMPLEMENTATION
#include "vendor/nob/nob.h"

#define STBI_ONLY_PNG
#define STBI_NO_HDR
#define STBI_NO_LINEAR
#define STBI_NO_SIMD
#define STB_IMAGE_IMPLEMENTATION
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#endif
#include "vendor/stb/stb_image.h"
#if defined(__clang__)
#pragma clang diagnostic pop
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

#include <stdint.h>

#define NOB_HEADER "vendor/nob/nob.h"
#define STB_IMAGE_HEADER "vendor/stb/stb_image.h"
#define SOURCE_SVG "assets/icon/magnoom.svg"
#define ICON_DIR "assets/icon/"
#define TEMP_DIR "build/icons/"
#define ICONSET_DIR TEMP_DIR "magnoom.iconset/"
#define TEMP_ICO TEMP_DIR "magnoom.ico"
#define TEMP_ICNS TEMP_DIR "magnoom.icns"
#define TEMP_X11_HEADER TEMP_DIR "magnoom_x11_icon.h"
#define PUBLISH_PREPARED ICON_DIR ".magnoom-publish-prepared"
#define PUBLISH_PREPARED_NEW ICON_DIR ".magnoom-publish-prepared.new"
#define PUBLISH_COMMITTED ICON_DIR ".magnoom-publish-committed"
#define PUBLISH_COMMITTED_NEW ICON_DIR ".magnoom-publish-committed.new"

typedef struct {
    int size;
    const char *path;
} icon_export;

static const icon_export exports[] = {
    {  16, TEMP_DIR "magnoom-16.png" },
    {  32, TEMP_DIR "magnoom-32.png" },
    {  48, TEMP_DIR "magnoom-48.png" },
    {  64, TEMP_DIR "magnoom-64.png" },
    { 128, TEMP_DIR "magnoom-128.png" },
    { 256, TEMP_DIR "magnoom-256.png" },
    { 512, TEMP_DIR "magnoom-512.png" },
    {1024, TEMP_DIR "magnoom-1024.png"},
};

static bool delete_if_exists(const char *path)
{
    return !nob_file_exists(path) || nob_delete_file(path);
}

static bool write_complete_file(const char *path, const void *data, size_t size)
{
    FILE *file = fopen(path, "wb");
    if (!file) return false;
    bool ok = fwrite(data, size, 1, file) == 1;
    if (fclose(file) != 0) ok = false;
    if (!ok) delete_if_exists(path);
    return ok;
}

static bool clean_temp_files(void)
{
    static const char *iconset_files[] = {
        "icon_16x16.png", "icon_16x16@2x.png",
        "icon_32x32.png", "icon_32x32@2x.png",
        "icon_128x128.png", "icon_128x128@2x.png",
        "icon_256x256.png", "icon_256x256@2x.png",
        "icon_512x512.png", "icon_512x512@2x.png",
    };

    for (size_t i = 0; i < NOB_ARRAY_LEN(exports); ++i) {
        if (!delete_if_exists(exports[i].path)) return false;
    }
    if (!delete_if_exists(TEMP_ICO) || !delete_if_exists(TEMP_ICNS) ||
        !delete_if_exists(TEMP_X11_HEADER)) return false;
    for (size_t i = 0; i < NOB_ARRAY_LEN(iconset_files); ++i) {
        if (!delete_if_exists(nob_temp_sprintf("%s%s", ICONSET_DIR, iconset_files[i])))
            return false;
    }
    if (!delete_if_exists(ICONSET_DIR) || !delete_if_exists(TEMP_DIR)) return false;

    return true;
}

static bool ensure_directories(void)
{
    return nob_mkdir_if_not_exists("build/") &&
           nob_mkdir_if_not_exists(TEMP_DIR) &&
           nob_mkdir_if_not_exists(ICONSET_DIR);
}

static const char *find_inkscape(void)
{
    const char *configured = getenv("INKSCAPE");
    if (configured && configured[0] != '\0') return configured;
#if defined(__APPLE__)
    if (nob_file_exists("/Applications/Inkscape.app/Contents/MacOS/inkscape"))
        return "/Applications/Inkscape.app/Contents/MacOS/inkscape";
#endif
    return "inkscape";
}

static bool output_needs_rebuild(const char *output)
{
    const char *inputs[] = { SOURCE_SVG, "nob_icons.c", NOB_HEADER, STB_IMAGE_HEADER };
    int result = nob_needs_rebuild(output, inputs, NOB_ARRAY_LEN(inputs));
    if (result < 0) exit(1);
    return result > 0;
}

static bool icons_need_rebuild(void)
{
    const char *publication_files[] = {
        PUBLISH_PREPARED,
        PUBLISH_PREPARED_NEW,
        PUBLISH_COMMITTED,
        PUBLISH_COMMITTED_NEW,
        ICON_DIR ".magnoom-256.png.bak",
        ICON_DIR ".magnoom-512.png.bak",
        ICON_DIR ".magnoom.ico.bak",
        ICON_DIR ".magnoom_x11_icon.h.bak",
#if defined(__APPLE__)
        ICON_DIR ".magnoom.icns.bak",
#endif
    };
    const char *outputs[] = {
        ICON_DIR "magnoom-256.png",
        ICON_DIR "magnoom-512.png",
        ICON_DIR "magnoom.ico",
        ICON_DIR "magnoom_x11_icon.h",
#if defined(__APPLE__)
        ICON_DIR "magnoom.icns",
#endif
    };

    for (size_t i = 0; i < NOB_ARRAY_LEN(publication_files); ++i) {
        if (nob_file_exists(publication_files[i])) return true;
    }
    for (size_t i = 0; i < NOB_ARRAY_LEN(outputs); ++i) {
        if (output_needs_rebuild(outputs[i])) return true;
    }
    return false;
}

static bool export_png(const char *inkscape, const icon_export *icon)
{
    Nob_Cmd cmd = {0};
    nob_cmd_append(&cmd, inkscape, SOURCE_SVG,
                   "--export-type=png", "--export-area-page",
                   "--export-background-opacity=0",
                   nob_temp_sprintf("--export-width=%d", icon->size),
                   nob_temp_sprintf("--export-height=%d", icon->size),
                   nob_temp_sprintf("--export-filename=%s", icon->path));
    return nob_cmd_run(&cmd, .max_procs = 0);
}

static bool validate_png(const icon_export *icon)
{
    int width = 0;
    int height = 0;
    int channels = 0;
    unsigned char *pixels = stbi_load(icon->path, &width, &height, &channels, STBI_rgb_alpha);
    if (!pixels) {
        nob_log(NOB_ERROR, "could not decode %s: %s", icon->path, stbi_failure_reason());
        return false;
    }

    bool transparent = false;
    for (size_t i = 0; i < (size_t)width * (size_t)height; ++i) {
        if (pixels[4*i + 3] < 255) {
            transparent = true;
            break;
        }
    }
    stbi_image_free(pixels);

    if (width != icon->size || height != icon->size) {
        nob_log(NOB_ERROR, "%s is %dx%d; expected %dx%d",
                icon->path, width, height, icon->size, icon->size);
        return false;
    }
    if (!transparent) {
        nob_log(NOB_ERROR, "%s has no transparent pixels", icon->path);
        return false;
    }
    return true;
}

static bool write_u16_le(FILE *file, uint16_t value)
{
    unsigned char bytes[] = {
        (unsigned char)(value & 0xffu),
        (unsigned char)((value >> 8) & 0xffu),
    };
    return fwrite(bytes, sizeof(bytes), 1, file) == 1;
}

static bool write_u32_le(FILE *file, uint32_t value)
{
    unsigned char bytes[] = {
        (unsigned char)(value & 0xffu),
        (unsigned char)((value >> 8) & 0xffu),
        (unsigned char)((value >> 16) & 0xffu),
        (unsigned char)((value >> 24) & 0xffu),
    };
    return fwrite(bytes, sizeof(bytes), 1, file) == 1;
}

static bool write_ico(void)
{
    enum { ICO_IMAGE_COUNT = 6 };
    Nob_String_Builder pngs[ICO_IMAGE_COUNT] = {0};
    uint32_t offsets[ICO_IMAGE_COUNT] = {0};
    uint32_t offset = 6u + ICO_IMAGE_COUNT * 16u;
    bool ok = false;
    FILE *file = NULL;

    for (size_t i = 0; i < ICO_IMAGE_COUNT; ++i) {
        if (!nob_read_entire_file(exports[i].path, &pngs[i])) {
            nob_log(NOB_ERROR, "could not read %s", exports[i].path);
            goto cleanup;
        }
        if (pngs[i].count > UINT32_MAX - offset) {
            nob_log(NOB_ERROR, "%s is too large for an ICO file", exports[i].path);
            goto cleanup;
        }
        offsets[i] = offset;
        offset += (uint32_t)pngs[i].count;
    }

    file = fopen(TEMP_ICO, "wb");
    if (!file) {
        nob_log(NOB_ERROR, "could not open %s", TEMP_ICO);
        goto cleanup;
    }

    ok = write_u16_le(file, 0) && write_u16_le(file, 1) &&
         write_u16_le(file, ICO_IMAGE_COUNT);
    for (size_t i = 0; ok && i < ICO_IMAGE_COUNT; ++i) {
        int size = exports[i].size;
        ok = fputc(size == 256 ? 0 : size, file) != EOF &&
             fputc(size == 256 ? 0 : size, file) != EOF &&
             fputc(0, file) != EOF && fputc(0, file) != EOF &&
             write_u16_le(file, 1) && write_u16_le(file, 32) &&
             write_u32_le(file, (uint32_t)pngs[i].count) &&
             write_u32_le(file, offsets[i]);
    }
    for (size_t i = 0; ok && i < ICO_IMAGE_COUNT; ++i)
        ok = fwrite(pngs[i].items, pngs[i].count, 1, file) == 1;

    if (!ok) nob_log(NOB_ERROR, "could not write %s", TEMP_ICO);

cleanup:
    if (file && fclose(file) != 0) ok = false;
    for (size_t i = 0; i < ICO_IMAGE_COUNT; ++i) free(pngs[i].items);
    return ok;
}

static bool write_x11_header(void)
{
    const icon_export *icon = &exports[3];
    int width = 0;
    int height = 0;
    int channels = 0;
    unsigned char *pixels = stbi_load(icon->path, &width, &height, &channels, STBI_rgb_alpha);
    if (!pixels) {
        nob_log(NOB_ERROR, "could not decode %s: %s", icon->path, stbi_failure_reason());
        return false;
    }
    if (width != 64 || height != 64) {
        nob_log(NOB_ERROR, "%s is %dx%d; expected 64x64", icon->path, width, height);
        stbi_image_free(pixels);
        return false;
    }

    Nob_String_Builder output = {0};
    nob_sb_append_cstr(&output,
        "#ifndef MAGNOOM_X11_ICON_H\n"
        "#define MAGNOOM_X11_ICON_H\n\n"
        "static const unsigned long magnoom_x11_icon[] = {\n"
        "    64UL, 64UL,\n");

    for (size_t i = 0; i < 64u * 64u; ++i) {
        if (i % 8 == 0) nob_sb_append_cstr(&output, "    ");
        uint32_t argb = ((uint32_t)pixels[4*i + 3] << 24) |
                        ((uint32_t)pixels[4*i + 0] << 16) |
                        ((uint32_t)pixels[4*i + 1] << 8) |
                        (uint32_t)pixels[4*i + 2];
        nob_sb_appendf(&output, "0x%08xUL,", argb);
        nob_sb_append_cstr(&output, i % 8 == 7 ? "\n" : " ");
    }
    nob_sb_append_cstr(&output, "};\n\n#endif\n");

    bool ok = nob_write_entire_file(TEMP_X11_HEADER, output.items, output.count);
    if (!ok) nob_log(NOB_ERROR, "could not write %s", TEMP_X11_HEADER);
    free(output.items);
    stbi_image_free(pixels);
    return ok;
}

#if defined(__APPLE__)
static bool copy_iconset_file(size_t export_index, const char *name)
{
    return nob_copy_file(exports[export_index].path,
                         nob_temp_sprintf("%s%s", ICONSET_DIR, name));
}

static bool write_icns(void)
{
    if (!copy_iconset_file(0, "icon_16x16.png") ||
        !copy_iconset_file(1, "icon_16x16@2x.png") ||
        !copy_iconset_file(1, "icon_32x32.png") ||
        !copy_iconset_file(3, "icon_32x32@2x.png") ||
        !copy_iconset_file(4, "icon_128x128.png") ||
        !copy_iconset_file(5, "icon_128x128@2x.png") ||
        !copy_iconset_file(5, "icon_256x256.png") ||
        !copy_iconset_file(6, "icon_256x256@2x.png") ||
        !copy_iconset_file(6, "icon_512x512.png") ||
        !copy_iconset_file(7, "icon_512x512@2x.png")) return false;

    if (nob_file_exists(TEMP_ICNS) && !nob_delete_file(TEMP_ICNS)) return false;
    const char *iconutil = getenv("ICONUTIL");
    if (!iconutil || iconutil[0] == '\0') iconutil = "iconutil";
    Nob_Cmd cmd = {0};
    nob_cmd_append(&cmd, iconutil, "-c", "icns", ICONSET_DIR, "-o", TEMP_ICNS);
    return nob_cmd_run(&cmd, .max_procs = 0);
}
#endif

typedef struct {
    const char *generated;
    const char *destination;
    const char *staged;
    const char *backup;
} published_icon;

enum { PUBLISHED_ICON_CAPACITY = 5 };

static size_t get_published_icons(published_icon icons[PUBLISHED_ICON_CAPACITY])
{
    size_t count = 0;
    icons[count++] = (published_icon){exports[5].path, ICON_DIR "magnoom-256.png",
                                     ICON_DIR ".magnoom-256.png.new",
                                     ICON_DIR ".magnoom-256.png.bak"};
    icons[count++] = (published_icon){exports[6].path, ICON_DIR "magnoom-512.png",
                                     ICON_DIR ".magnoom-512.png.new",
                                     ICON_DIR ".magnoom-512.png.bak"};
    icons[count++] = (published_icon){TEMP_ICO, ICON_DIR "magnoom.ico",
                                     ICON_DIR ".magnoom.ico.new",
                                     ICON_DIR ".magnoom.ico.bak"};
    icons[count++] = (published_icon){TEMP_X11_HEADER, ICON_DIR "magnoom_x11_icon.h",
                                     ICON_DIR ".magnoom_x11_icon.h.new",
                                     ICON_DIR ".magnoom_x11_icon.h.bak"};
#if defined(__APPLE__)
    icons[count++] = (published_icon){TEMP_ICNS, ICON_DIR "magnoom.icns",
                                     ICON_DIR ".magnoom.icns.new",
                                     ICON_DIR ".magnoom.icns.bak"};
#endif
    return count;
}

static bool read_prepared_mask(uint32_t *mask)
{
    FILE *file = fopen(PUBLISH_PREPARED, "rb");
    char buffer[32];
    char expected[32];
    unsigned int value = 0;
    if (!file) return false;
    bool ok = fgets(buffer, sizeof(buffer), file) != NULL && fgetc(file) == EOF && !ferror(file);
    if (fclose(file) != 0) ok = false;
    if (!ok || sscanf(buffer, "prepared %x", &value) != 1) return false;
    int length = snprintf(expected, sizeof(expected), "prepared %x\n", value);
    if (length < 0 || (size_t)length >= sizeof(expected) || strcmp(buffer, expected) != 0)
        return false;
    *mask = (uint32_t)value;
    return ok;
}

static bool rollback_publication(const published_icon *icons, size_t icon_count,
                                 const bool *had_original, const bool *installed)
{
    bool ok = true;
    for (size_t i = 0; i < icon_count; ++i) {
        if (had_original[i] && nob_file_exists(icons[i].backup)) {
            if (!nob_rename(icons[i].backup, icons[i].destination)) ok = false;
        } else if (!had_original[i] && installed[i]) {
            if (!delete_if_exists(icons[i].destination)) ok = false;
        }
        if (!delete_if_exists(icons[i].staged)) ok = false;
    }

    if (ok) {
        if (!delete_if_exists(PUBLISH_PREPARED_NEW) ||
            !delete_if_exists(PUBLISH_COMMITTED_NEW) ||
            !delete_if_exists(PUBLISH_PREPARED) ||
            !delete_if_exists(PUBLISH_COMMITTED)) ok = false;
    }
    return ok;
}

static bool recover_publication(void)
{
    published_icon icons[PUBLISHED_ICON_CAPACITY];
    size_t icon_count = get_published_icons(icons);
    bool interrupted = false;
    for (size_t i = 0; i < icon_count; ++i) {
        if (nob_file_exists(icons[i].backup)) interrupted = true;
    }

    if (nob_file_exists(PUBLISH_COMMITTED)) {
        for (size_t i = 0; i < icon_count; ++i) {
            if (!delete_if_exists(icons[i].backup)) return false;
        }
        if (!delete_if_exists(PUBLISH_PREPARED) ||
            !delete_if_exists(PUBLISH_COMMITTED)) return false;
    } else if (nob_file_exists(PUBLISH_PREPARED)) {
        uint32_t original_mask = 0;
        if (!read_prepared_mask(&original_mask) ||
            (original_mask >> icon_count) != 0) {
            nob_log(NOB_ERROR, "invalid icon publication state in %s", PUBLISH_PREPARED);
            return false;
        }
        nob_log(NOB_WARNING, "restoring icon assets after interrupted publication");
        for (size_t i = 0; i < icon_count; ++i) {
            if ((original_mask & (1u << i)) != 0) {
                if (nob_file_exists(icons[i].backup)) {
                    if (!delete_if_exists(icons[i].destination) ||
                        !nob_rename(icons[i].backup, icons[i].destination)) return false;
                } else if (!nob_file_exists(icons[i].destination)) {
                    nob_log(NOB_ERROR, "lost original icon asset %s", icons[i].destination);
                    return false;
                }
            } else if (!delete_if_exists(icons[i].destination)) {
                return false;
            }
        }
        if (!delete_if_exists(PUBLISH_PREPARED)) return false;
    } else if (interrupted) {
        nob_log(NOB_WARNING, "restoring icon assets after interrupted publication");
        for (size_t i = 0; i < icon_count; ++i) {
            if (!nob_file_exists(icons[i].backup)) continue;
            if (!delete_if_exists(icons[i].destination) ||
                !nob_rename(icons[i].backup, icons[i].destination)) return false;
        }
    }

    if (!delete_if_exists(PUBLISH_PREPARED_NEW) ||
        !delete_if_exists(PUBLISH_COMMITTED_NEW)) return false;
    for (size_t i = 0; i < icon_count; ++i) {
        if (!delete_if_exists(icons[i].staged)) return false;
    }
    return true;
}

static bool publish_icons(void)
{
    published_icon icons[PUBLISHED_ICON_CAPACITY];
    size_t icon_count = get_published_icons(icons);
    bool had_original[PUBLISHED_ICON_CAPACITY];
    bool installed[PUBLISHED_ICON_CAPACITY];
    memset(had_original, 0, sizeof(had_original));
    memset(installed, 0, sizeof(installed));

    if (!recover_publication()) return false;

    // Copying may fail or truncate, so stage every file before moving originals.
    for (size_t i = 0; i < icon_count; ++i) {
        if (!nob_copy_file(icons[i].generated, icons[i].staged)) {
            for (size_t j = 0; j < icon_count; ++j)
                delete_if_exists(icons[j].staged);
            return false;
        }
    }

    uint32_t original_mask = 0;
    for (size_t i = 0; i < icon_count; ++i) {
        had_original[i] = nob_file_exists(icons[i].destination);
        if (had_original[i]) original_mask |= 1u << i;
    }
    char prepared[32];
    int prepared_length = snprintf(prepared, sizeof(prepared),
                                   "prepared %x\n", original_mask);
    if (prepared_length < 0 || (size_t)prepared_length >= sizeof(prepared) ||
        !write_complete_file(PUBLISH_PREPARED_NEW, prepared, (size_t)prepared_length) ||
        !nob_rename(PUBLISH_PREPARED_NEW, PUBLISH_PREPARED)) {
        for (size_t i = 0; i < icon_count; ++i)
            delete_if_exists(icons[i].staged);
        delete_if_exists(PUBLISH_PREPARED_NEW);
        delete_if_exists(PUBLISH_PREPARED);
        return false;
    }

    for (size_t i = 0; i < icon_count; ++i) {
        if (had_original[i] && !nob_rename(icons[i].destination, icons[i].backup)) {
            rollback_publication(icons, icon_count, had_original, installed);
            return false;
        }
    }

    for (size_t i = 0; i < icon_count; ++i) {
        if (!nob_rename(icons[i].staged, icons[i].destination)) {
            rollback_publication(icons, icon_count, had_original, installed);
            return false;
        }
        installed[i] = true;
    }

    static const char committed[] = "committed\n";
    if (!write_complete_file(PUBLISH_COMMITTED_NEW, committed, sizeof(committed) - 1) ||
        !nob_rename(PUBLISH_COMMITTED_NEW, PUBLISH_COMMITTED)) {
        rollback_publication(icons, icon_count, had_original, installed);
        return false;
    }

    for (size_t i = 0; i < icon_count; ++i) {
        if (!delete_if_exists(icons[i].backup)) return false;
    }
    if (!delete_if_exists(PUBLISH_PREPARED) ||
        !delete_if_exists(PUBLISH_COMMITTED) ||
        !delete_if_exists(PUBLISH_PREPARED_NEW) ||
        !delete_if_exists(PUBLISH_COMMITTED_NEW)) return false;

#if !defined(__APPLE__)
    nob_log(NOB_WARNING,
            "preserving existing %smagnoom.icns; regenerate on macOS before release",
            ICON_DIR);
#endif
    return true;
}

static bool generate_icons(void)
{
    const char *inkscape = find_inkscape();
    Nob_Cmd version = {0};
    nob_cmd_append(&version, inkscape, "--version");
    if (!nob_cmd_run(&version, .max_procs = 0)) {
        nob_log(NOB_ERROR, "Inkscape was not found; set INKSCAPE to its executable path");
        return false;
    }
    if (!ensure_directories()) return false;

    for (size_t i = 0; i < NOB_ARRAY_LEN(exports); ++i) {
        if (!export_png(inkscape, &exports[i]) || !validate_png(&exports[i])) return false;
    }
    if (!write_ico() || !write_x11_header()) return false;
#if defined(__APPLE__)
    if (!write_icns()) return false;
#endif
    return publish_icons();
}

static void usage(const char *program)
{
    printf("usage: %s [-force | -clean | -help]\n", program);
    printf("  -force  regenerate every icon asset\n");
    printf("  -clean  remove temporary icon-export files\n");
    printf("  -help   print this help and exit\n");
}

int main(int argc, char **argv)
{
    NOB_GO_REBUILD_URSELF_PLUS(argc, argv, NOB_HEADER, STB_IMAGE_HEADER);

    bool force = false;
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-force") == 0) {
            force = true;
        } else if (strcmp(argv[i], "-clean") == 0) {
            return clean_temp_files() ? 0 : 1;
        } else if (strcmp(argv[i], "-help") == 0 || strcmp(argv[i], "--help") == 0) {
            usage(argv[0]);
            return 0;
        } else {
            nob_log(NOB_ERROR, "unknown argument: %s", argv[i]);
            usage(argv[0]);
            return 1;
        }
    }

    if (!recover_publication()) return 1;
    if (!force && !icons_need_rebuild()) {
        nob_log(NOB_INFO, "icon assets are up to date");
        return 0;
    }
    if (!generate_icons()) return 1;
    nob_log(NOB_INFO, "generated icon assets from %s", SOURCE_SVG);
    return 0;
}
