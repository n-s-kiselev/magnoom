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
    test_default_block();
    test_crystal_examples();
    test_invalid_inputs_are_transactional();
    if (failures != 0) {
        fprintf(stderr, "%d block setter test(s) failed.\n", failures);
        return 1;
    }
    printf("Block setter tests passed.\n");
    return 0;
}
