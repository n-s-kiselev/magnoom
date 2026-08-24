#define MAGNOOM_NO_MAIN
#include "../magnoom.c"

static int failures;

#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "%s:%d: check failed: %s\n", __FILE__, __LINE__, #condition); \
        failures += 1; \
    } \
} while (0)

static void free_geometry(magnoom_ctx *ctx)
{
    free(ctx->Block);
    free(ctx->RadiusOfShell);
    free(ctx->NeighborsPerAtom);
}

static void set_identity_cell(magnoom_ctx *ctx)
{
    const float identity[3][3] = {
        {1.0f, 0.0f, 0.0f},
        {0.0f, 1.0f, 0.0f},
        {0.0f, 0.0f, 1.0f}
    };
    memcpy(ctx->abc, identity, sizeof(identity));
}

static void test_path_helpers(void)
{
    char path[128];
    char executable_directory[MAGNOOM_PATH_CAPACITY];

#if defined(_WIN32)
    CHECK(magnoom_join_path(path, sizeof(path), "C:\\data", "input.csv"));
    CHECK(strcmp(path, "C:\\data\\input.csv") == 0);
    CHECK(magnoom_path_is_absolute("C:\\data\\input.csv"));
    CHECK(magnoom_path_is_absolute("\\\\server\\share\\input.csv"));
    CHECK(!magnoom_path_is_absolute("C:input.csv"));
    CHECK(!magnoom_join_path(path, sizeof(path), "D:\\data", "C:input.csv"));
    CHECK(!magnoom_join_path(path, sizeof(path), "D:\\data", "\\input.csv"));
    CHECK(magnoom_join_path(path, sizeof(path), "D:\\ignored", "C:\\absolute\\input.csv"));
    CHECK(strcmp(path, "C:\\absolute\\input.csv") == 0);
#else
    CHECK(magnoom_join_path(path, sizeof(path), "/tmp/data", "input.csv"));
    CHECK(strcmp(path, "/tmp/data/input.csv") == 0);
    CHECK(magnoom_path_is_absolute("/tmp/data/input.csv"));
    CHECK(!magnoom_path_is_absolute("tmp/data/input.csv"));
    CHECK(magnoom_join_path(path, sizeof(path), "/ignored", "/absolute/input.csv"));
    CHECK(strcmp(path, "/absolute/input.csv") == 0);
#endif
    CHECK(magnoom_join_path(path, sizeof(path), "/tmp/data/", "input.csv"));
    CHECK(strcmp(path, "/tmp/data/input.csv") == 0);
    CHECK(magnoom_replace_extension(path, sizeof(path),
        "/tmp/user.name/output.csv", ".ovf"));
    CHECK(strcmp(path, "/tmp/user.name/output.ovf") == 0);
    CHECK(magnoom_replace_extension(path, sizeof(path), "/tmp/.state", ".vtk"));
    CHECK(strcmp(path, "/tmp/.state.vtk") == 0);
    CHECK(!magnoom_join_path(path, 8, "/tmp/data", "input.csv"));
    CHECK(!magnoom_join_path(path, sizeof(path), "/tmp/data", ""));
    CHECK(magnoom_executable_directory(executable_directory,
        sizeof(executable_directory), "./block_setter_test"));
    CHECK(executable_directory[0] != '\0');
}

static void test_solver_state_reset(void)
{
    magnoom_ctx ctx = {0};
    double spins[9] = {2.0, -1.0, 0.5, 0.0, 1.0, 0.0, -0.25, 0.75, -2.0};
    double original[9];
    double display[9] = {NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN};
    double stage[9] = {NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN};
    double increments[9] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double field_x[3] = {2.0, 2.0, 2.0};
    double field_y[3] = {2.0, 2.0, 2.0};
    double field_z[3] = {2.0, 2.0, 2.0};
    float random_x[3] = {NAN, NAN, NAN};
    float random_y[3] = {NAN, NAN, NAN};
    float random_z[3] = {NAN, NAN, NAN};

    memcpy(original, spins, sizeof(spins));

    ctx.NOS = 3;
    ctx.S = spins;
    ctx.bS = display;
    ctx.tS = stage;
    ctx.t2S = increments;
    ctx.Heffx = field_x;
    ctx.Heffy = field_y;
    ctx.Heffz = field_z;
    ctx.RNx = random_x;
    ctx.RNy = random_y;
    ctx.RNz = random_z;

    magnoom_reset_solver_state(&ctx);
    for (int i = 0; i < 9; ++i) {
        CHECK(spins[i] == original[i]);
        CHECK(display[i] == spins[i]);
        CHECK(stage[i] == spins[i]);
        CHECK(increments[i] == 1.0);
    }
    for (int i = 0; i < 3; ++i) {
        CHECK(field_x[i] == 2.0 && field_y[i] == 2.0 && field_z[i] == 2.0);
        CHECK(random_x[i] == 0.0f && random_y[i] == 0.0f && random_z[i] == 0.0f);
    }
}

static void check_cuboid_geometry(const float *vertices, const float *normals, const GLuint *indices)
{
    float center[3] = {0.0f, 0.0f, 0.0f};
    for (int vertex = 0; vertex < 24; ++vertex) {
        for (int component = 0; component < 3; ++component) {
            center[component] += vertices[3*vertex + component] / 24.0f;
        }
    }

    for (int face = 0; face < 6; ++face) {
        const int base = face * 12;
        float face_center[3] = {0.0f, 0.0f, 0.0f};
        for (int vertex = 0; vertex < 4; ++vertex) {
            for (int component = 0; component < 3; ++component) {
                face_center[component] += vertices[base + 3*vertex + component] / 4.0f;
                CHECK(fabsf(normals[base + 3*vertex + component] - normals[base + component]) < 1e-6f);
            }
        }
        float dot = 0.0f;
        for (int component = 0; component < 3; ++component) {
            dot += normals[base + component] * (face_center[component] - center[component]);
        }
        CHECK(dot > 0.0f);
        float normal_length_squared =
            normals[base + 0]*normals[base + 0] +
            normals[base + 1]*normals[base + 1] +
            normals[base + 2]*normals[base + 2];
        CHECK(fabsf(normal_length_squared - 1.0f) < 1e-5f);
    }

    for (int triangle = 0; triangle < 12; ++triangle) {
        GLuint i0 = indices[3*triangle + 0];
        GLuint i1 = indices[3*triangle + 1];
        GLuint i2 = indices[3*triangle + 2];
        CHECK(i0 < 24 && i1 < 24 && i2 < 24);
        if (i0 < 24 && i1 < 24 && i2 < 24) {
            float edge1[3], edge2[3], triangle_normal[3];
            for (int component = 0; component < 3; ++component) {
                edge1[component] = vertices[3*i1 + component] - vertices[3*i0 + component];
                edge2[component] = vertices[3*i2 + component] - vertices[3*i0 + component];
            }
            Crossf(edge1, edge2, triangle_normal);
            float winding_dot =
                triangle_normal[0]*normals[3*i0 + 0] +
                triangle_normal[1]*normals[3*i0 + 1] +
                triangle_normal[2]*normals[3*i0 + 2];
            CHECK(winding_dot < 0.0f);
        }
    }
}

static void test_cuboid_normals(void)
{
    magnoom_ctx ctx = {0};
    const float abc[3][3] = {
        {1.0f, 0.0f, 0.0f},
        {0.2f, 1.0f, 0.0f},
        {0.1f, 0.3f, 1.0f}
    };
    memcpy(ctx.abc, abc, sizeof(abc));

    float vertices[6*4*3] = {0};
    float normals[6*4*3] = {0};
    float colors[6*4*3] = {0};
    GLuint indices[6*2*3] = {0};
    FillProtoVerNorInd(&ctx, vertices, normals, indices, 4, BOX1, 0);
    check_cuboid_geometry(vertices, normals, indices);

    memset(vertices, 0, sizeof(vertices));
    memset(normals, 0, sizeof(normals));
    const float translation[3] = {0.3f, -0.2f, 0.5f};
    const float color[3] = {1.0f, 1.0f, 1.0f};
    parallelepiped(abc, translation, 1.7f, 0.4f, 0.8f,
        color, 0, vertices, normals, colors, indices);
    check_cuboid_geometry(vertices, normals, indices);

    ctx.abc[2][0] = -ctx.abc[2][0];
    ctx.abc[2][1] = -ctx.abc[2][1];
    ctx.abc[2][2] = -ctx.abc[2][2];
    memset(vertices, 0, sizeof(vertices));
    memset(normals, 0, sizeof(normals));
    memset(indices, 0, sizeof(indices));
    FillProtoVerNorInd(&ctx, vertices, normals, indices, 4, BOX1, 0);
    check_cuboid_geometry(vertices, normals, indices);

    memset(vertices, 0, sizeof(vertices));
    memset(normals, 0, sizeof(normals));
    memset(indices, 0, sizeof(indices));
    parallelepiped(ctx.abc, translation, 1.7f, 0.4f, 0.8f,
        color, 0, vertices, normals, colors, indices);
    check_cuboid_geometry(vertices, normals, indices);
}

static void test_default_block(void)
{
    magnoom_ctx ctx;
    memset(&ctx, 0xA5, sizeof(ctx));
    CHECK(magnoom_ctx_init(&ctx));
    CHECK(ctx.AtomsPerBlock == 1);
    CHECK(ctx.Block != NULL);
    CHECK(ctx.NeighborsPerAtom != NULL);
    CHECK(ctx.Block[0][0] == 0.5f);
    CHECK(ctx.Block[0][1] == 0.5f);
    CHECK(ctx.Block[0][2] == 0.5f);
    CHECK(ctx.NOB == 1000);
    CHECK(ctx.NOS == 1000);
    CHECK(ctx.WhichLightingMode == LIGHT_ADAPTIVE);
    CHECK(ctx.shininess >= 0.0f && ctx.shininess <= 128.0f);
    free_geometry(&ctx);
}

static void test_crystal_examples(void)
{
    magnoom_ctx ctx = {0};
    CHECK(magnoom_ctx_init(&ctx));

    const float uB20 = 0.138f;
    float b20[][3] = {
        {0.0f,             0.0f,             0.0f},
        {0.5f,             0.5f-2.0f*uB20, 1.0f-2.0f*uB20},
        {1.0f-2.0f*uB20, 0.5f,             0.5f-2.0f*uB20},
        {0.5f-2.0f*uB20, 1.0f-2.0f*uB20, 0.5f}
    };
    set_identity_cell(&ctx);
    ctx.RadiusOfShell[0] = 7.0f;
    CHECK(magnoom_ctx_set_block(&ctx, 4, b20));
    CHECK(ctx.AtomsPerBlock == 4);
    CHECK(ctx.NOS == 4000);
    CHECK(ctx.RadiusOfShell[0] == 7.0f);
    CHECK(ctx.Block[3][0] == b20[3][0]);
    for (int atom = 0; atom < ctx.AtomsPerBlock; ++atom) {
        CHECK(ctx.NeighborsPerAtom[atom] == 0);
    }
    b20[0][0] = 0.25f;
    CHECK(ctx.Block[0][0] == 0.0f);

    const float a_EuSi = 3.9845f;
    const float b_EuSi = 4.6955f/a_EuSi;
    const float c_EuSi = 11.1528f/a_EuSi;
    const float u_Eu = 0.3595f;
    const float eusi_abc[3][3] = {
        {1.0f, 0.0f, 0.0f},
        {0.0f, b_EuSi, 0.0f},
        {0.0f, 0.0f, c_EuSi}
    };
    const float eusi[][3] = {
        {0.25f, 0.0f,          u_Eu*c_EuSi},
        {0.75f, 0.0f,          (1.0f-u_Eu)*c_EuSi},
        {0.75f, 0.5f*b_EuSi, (0.5f-u_Eu)*c_EuSi},
        {0.25f, 0.5f*b_EuSi, (0.5f+u_Eu)*c_EuSi}
    };
    memcpy(ctx.abc, eusi_abc, sizeof(eusi_abc));
    CHECK(magnoom_ctx_set_block(&ctx, 4, eusi));
    CHECK(ctx.abc[0][0] == 1.0f);
    CHECK(ctx.AtomsPerBlock == 4);
    CHECK(ctx.NOS == 4000);
    CHECK(ctx.Block[3][2] == eusi[3][2]);

    GetShells(&ctx);
    int neighbor_pairs = GetNeighborsNumber(
        ctx.abc, ctx.Block, ctx.AtomsPerBlock, ctx.ShellNumber,
        ctx.RadiusOfShell, ctx.NeighborsPerAtom);
    CHECK(neighbor_pairs == 8);
    for (int atom = 0; atom < ctx.AtomsPerBlock; ++atom) {
        CHECK(ctx.NeighborsPerAtom[atom] == 2);
    }

    int *atom_indices = (int *)calloc((size_t)neighbor_pairs, sizeof(int));
    int *neighbor_indices = (int *)calloc((size_t)neighbor_pairs, sizeof(int));
    int *grid_a = (int *)calloc((size_t)neighbor_pairs, sizeof(int));
    int *grid_b = (int *)calloc((size_t)neighbor_pairs, sizeof(int));
    int *grid_c = (int *)calloc((size_t)neighbor_pairs, sizeof(int));
    int *shell_indices = (int *)calloc((size_t)neighbor_pairs, sizeof(int));
    CHECK(atom_indices != NULL && neighbor_indices != NULL && grid_a != NULL &&
          grid_b != NULL && grid_c != NULL && shell_indices != NULL);
    if (atom_indices != NULL && neighbor_indices != NULL && grid_a != NULL &&
        grid_b != NULL && grid_c != NULL && shell_indices != NULL) {
        CreateMapOfNeighbors(
            ctx.abc, ctx.Block, ctx.AtomsPerBlock, ctx.ShellNumber,
            ctx.RadiusOfShell, atom_indices, neighbor_indices,
            grid_a, grid_b, grid_c, shell_indices);
        for (int pair = 0; pair < neighbor_pairs; ++pair) {
            bool reciprocal_found = false;
            for (int candidate = 0; candidate < neighbor_pairs; ++candidate) {
                if (atom_indices[candidate] == neighbor_indices[pair] &&
                    neighbor_indices[candidate] == atom_indices[pair] &&
                    grid_a[candidate] == -grid_a[pair] &&
                    grid_b[candidate] == -grid_b[pair] &&
                    grid_c[candidate] == -grid_c[pair]) {
                    reciprocal_found = true;
                    break;
                }
            }
            CHECK(reciprocal_found);
            CHECK(shell_indices[pair] == 0);
        }
    }
    free(atom_indices);
    free(neighbor_indices);
    free(grid_a);
    free(grid_b);
    free(grid_c);
    free(shell_indices);

    const float sqrt2 = sqrtf(2.0f);
    const float fcc2_abc[3][3] = {
        {1.0f, 0.0f, 0.0f},
        {0.0f, sqrt2, 0.0f},
        {0.0f, 0.0f, sqrt2}
    };
    const float fcc2[][3] = {
        {0.0f, 0.0f,             0.0f},
        {0.0f, sqrt2/2.0f,      0.0f},
        {0.5f, sqrt2/4.0f,      sqrt2/4.0f},
        {0.5f, 3.0f*sqrt2/4.0f, sqrt2/4.0f},
        {0.0f, sqrt2/2.0f,      sqrt2/2.0f}
    };
    memcpy(ctx.abc, fcc2_abc, sizeof(fcc2_abc));
    CHECK(magnoom_ctx_set_block(&ctx, 5, fcc2));
    CHECK(ctx.AtomsPerBlock == 5);
    CHECK(ctx.NOS == 5000);
    CHECK(ctx.Block[4][2] == fcc2[4][2]);

    const float sqrt3 = sqrtf(3.0f);
    const float sqrt6 = sqrtf(6.0f);
    const float fcc3_abc[3][3] = {
        {sqrt2,      0.0f,       0.0f},
        {sqrt2/2.0f, sqrt6/2.0f, 0.0f},
        {sqrt2/2.0f, sqrt6/6.0f, sqrt3/sqrt2}
    };
    const float fcc3[][3] = {{0.0f, 0.0f, 0.0f}};
    memcpy(ctx.abc, fcc3_abc, sizeof(fcc3_abc));
    CHECK(magnoom_ctx_set_block(&ctx, 1, fcc3));
    CHECK(ctx.AtomsPerBlock == 1);
    CHECK(ctx.NOS == 1000);
    CHECK(ctx.Block[0][0] == 0.0f);

    const float nonorthogonal_point[][3] = {{
        0.25f*fcc3_abc[0][0] + 0.5f*fcc3_abc[1][0] + 0.75f*fcc3_abc[2][0],
        0.25f*fcc3_abc[0][1] + 0.5f*fcc3_abc[1][1] + 0.75f*fcc3_abc[2][1],
        0.25f*fcc3_abc[0][2] + 0.5f*fcc3_abc[1][2] + 0.75f*fcc3_abc[2][2]
    }};
    CHECK(magnoom_ctx_set_block(&ctx, 1, nonorthogonal_point));
    CHECK(ctx.Block[0][2] == nonorthogonal_point[0][2]);

    free_geometry(&ctx);
}

static void test_invalid_inputs_are_transactional(void)
{
    magnoom_ctx ctx = {0};
    CHECK(magnoom_ctx_init(&ctx));
    float (*old_block)[3] = ctx.Block;
    int *old_neighbors = ctx.NeighborsPerAtom;
    int old_atom_count = ctx.AtomsPerBlock;
    int old_nos = ctx.NOS;
    const float inside[][3] = {{0.25f, 0.25f, 0.25f}};
    const float upper_face[][3] = {{1.0f, 0.25f, 0.25f}};
    const float near_upper_face[][3] = {{0.9999f, 0.25f, 0.25f}};
    const float near_lower_outside[][3] = {{-FLT_EPSILON, 0.25f, 0.25f}};
    const float outside[][3] = {{-0.1f, 0.25f, 0.25f}};
    const float too_many_atoms[MAX_ATOMS_PER_BLOCK + 1][3] = {{0.0f, 0.0f, 0.0f}};

    CHECK(!magnoom_ctx_set_block(NULL, 1, inside));
    CHECK(!magnoom_ctx_set_block(&ctx, 0, inside));
    CHECK(!magnoom_ctx_set_block(&ctx, 1, NULL));
    CHECK(!magnoom_ctx_set_block(&ctx, 1, upper_face));
    CHECK(magnoom_ctx_set_block(&ctx, 1, near_upper_face));
    old_block = ctx.Block;
    old_neighbors = ctx.NeighborsPerAtom;
    old_atom_count = ctx.AtomsPerBlock;
    old_nos = ctx.NOS;
    CHECK(!magnoom_ctx_set_block(&ctx, 1, near_lower_outside));
    CHECK(!magnoom_ctx_set_block(&ctx, 1, outside));
    CHECK(!magnoom_ctx_set_block(&ctx, MAX_ATOMS_PER_BLOCK + 1, too_many_atoms));
    CHECK(ctx.Block == old_block);
    CHECK(ctx.NeighborsPerAtom == old_neighbors);
    CHECK(ctx.AtomsPerBlock == old_atom_count);
    CHECK(ctx.NOS == old_nos);

    float saved_abc[3][3];
    memcpy(saved_abc, ctx.abc, sizeof(saved_abc));
    ctx.abc[2][0] = 0.0f;
    ctx.abc[2][1] = 0.0f;
    ctx.abc[2][2] = 0.0f;
    CHECK(!magnoom_ctx_set_block(&ctx, 1, inside));
    CHECK(ctx.Block == old_block);
    memcpy(ctx.abc, saved_abc, sizeof(saved_abc));

    int saved_dimensions[3] = {ctx.uABC[0], ctx.uABC[1], ctx.uABC[2]};
    ctx.uABC[0] = INT_MAX;
    ctx.uABC[1] = 2;
    ctx.uABC[2] = 1;
    CHECK(!magnoom_ctx_set_block(&ctx, 1, inside));
    CHECK(ctx.Block == old_block);
    memcpy(ctx.uABC, saved_dimensions, sizeof(saved_dimensions));

    ctx.NeighborPairs = 1;
    CHECK(!magnoom_ctx_set_block(&ctx, 1, inside));
    CHECK(ctx.Block == old_block);
    ctx.NeighborPairs = 0;

    free_geometry(&ctx);
}

int main(void)
{
    test_path_helpers();
    test_solver_state_reset();
    test_default_block();
    test_cuboid_normals();
    test_crystal_examples();
    test_invalid_inputs_are_transactional();
    if (failures != 0) {
        fprintf(stderr, "%d Magnoom test(s) failed.\n", failures);
        return 1;
    }
    printf("Magnoom tests passed.\n");
    return 0;
}
