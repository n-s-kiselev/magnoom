/* Regression harness for the Sx/Sy/Sz -> interleaved-S vector-field migration.
 * Not part of the normal build (see MAGNOOM_TEST_BUILD in magnoom.c / nob.c's -test target).
 * Builds a small lattice headlessly (no GUI/threads), exercises allocation, a
 * deterministic initial state, one integrator step, and a file round-trip, then
 * prints fixed-format numeric results to stdout so before/after runs can be diffed. */
#define MAGNOOM_TEST_BUILD
#include "magnoom.c"

static void print_spin(const char *label, int n, double sx, double sy, double sz)
{
    printf("%s[%d] = (%.10f, %.10f, %.10f)\n", label, n, sx, sy, sz);
}

int main(void)
{
    static magnoom_ctx ctx; /* static: zero-initialized, matching the real mag_ctx global */

    if (!magnoom_ctx_init(&ctx)) {
        fprintf(stderr, "magnoom_ctx_init failed\n");
        return 1;
    }

    /* Build the (default 10x10x10) neighbor map exactly like main() does. */
    GetShells(ctx.abc, ctx.Block, ctx.AtomsPerBlock, ctx.ShellNumber, ctx.RadiusOfShell);
    ctx.NeighborPairs = GetNeighborsNumber(ctx.abc, ctx.Block, ctx.AtomsPerBlock, ctx.ShellNumber, ctx.RadiusOfShell, ctx.NeighborsPerAtom);
    ctx.AIdxBlock = (int *)calloc(ctx.NeighborPairs, sizeof(int));
    ctx.NIdxBlock = (int *)calloc(ctx.NeighborPairs, sizeof(int));
    ctx.NIdxGridA = (int *)calloc(ctx.NeighborPairs, sizeof(int));
    ctx.NIdxGridB = (int *)calloc(ctx.NeighborPairs, sizeof(int));
    ctx.NIdxGridC = (int *)calloc(ctx.NeighborPairs, sizeof(int));
    ctx.SIdx      = (int *)calloc(ctx.NeighborPairs, sizeof(int));
    CreateMapOfNeighbors(ctx.abc, ctx.Block, ctx.AtomsPerBlock, ctx.ShellNumber, ctx.RadiusOfShell, ctx.AIdxBlock, ctx.NIdxBlock, ctx.NIdxGridA, ctx.NIdxGridB, ctx.NIdxGridC, ctx.SIdx);
    ctx.Jexc = (float *)calloc(ctx.NeighborPairs, sizeof(float));
    ctx.Bexc = (float *)calloc(ctx.NeighborPairs, sizeof(float));
    ctx.Dexc = (float *)calloc(ctx.NeighborPairs, sizeof(float));
    ctx.VDMx = (float *)calloc(ctx.NeighborPairs, sizeof(float));
    ctx.VDMy = (float *)calloc(ctx.NeighborPairs, sizeof(float));
    ctx.VDMz = (float *)calloc(ctx.NeighborPairs, sizeof(float));
    SetExch1(ctx.abc, ctx.Block, ctx.NeighborPairs, ctx.Jij, ctx.Bij, ctx.Dij, ctx.AIdxBlock, ctx.NIdxBlock, ctx.NIdxGridA, ctx.NIdxGridB, ctx.NIdxGridC, ctx.SIdx,
              ctx.Jexc, ctx.Bexc, ctx.Dexc, ctx.VDMx, ctx.VDMy, ctx.VDMz);

    /* --- allocation/free test --- */
    ReallocateMemoryForSpins(&ctx, ctx.NOS);
    ReallocateMemoryForAllOther(&ctx, ctx.NOS);
    ReallocateMemoryForImages(&ctx, ctx.num_images, ctx.NOS);
    printf("alloc: all non-null = %d\n",
           ctx.S && ctx.bS && ctx.tS && ctx.t2S && ctx.t3S && ctx.Image && ctx.dImage);

    /* --- geometry / kind / deterministic (HOMO) initial state --- */
    GetBox(&ctx, ctx.abc, ctx.uABC, ctx.Box);
    UpdateSpinPositions(&ctx, ctx.abc, ctx.uABC, ctx.Block, ctx.AtomsPerBlock, ctx.Box, ctx.Px, ctx.Py, ctx.Pz);
    UpdateKind(&ctx, ctx.Kind, ctx.Px, ctx.Py, ctx.Pz, ctx.NOS, ctx.NOSK);
    InitSpinComponents(&ctx, ctx.Px, ctx.Py, ctx.Pz, ctx.S, 1 /* HOMO, deterministic */);
    for (int i = 0; i < ctx.NOS; i++) {
        VEC_X(ctx.bS,i) = VEC_X(ctx.S,i);
        VEC_Y(ctx.bS,i) = VEC_Y(ctx.S,i);
        VEC_Z(ctx.bS,i) = VEC_Z(ctx.S,i);
    }

    printf("NOS=%d\n", ctx.NOS);
    print_spin("initial", 0, VEC_X(ctx.S,0), VEC_Y(ctx.S,0), VEC_Z(ctx.S,0));
    print_spin("initial", ctx.NOS / 2, VEC_X(ctx.S,ctx.NOS/2), VEC_Y(ctx.S,ctx.NOS/2), VEC_Z(ctx.S,ctx.NOS/2));
    print_spin("initial", ctx.NOS - 1, VEC_X(ctx.S,ctx.NOS-1), VEC_Y(ctx.S,ctx.NOS-1), VEC_Z(ctx.S,ctx.NOS-1));

    /* --- copy S -> tS (memcpy-equivalent) --- */
    memcpy(ctx.tS, ctx.S, 3*(size_t)ctx.NOS * sizeof(double));
    printf("copy S->tS matches: %d\n",
           memcmp(ctx.tS, ctx.S, 3*(size_t)ctx.NOS*sizeof(double)) == 0);

    /* --- one deterministic integrator step (Temperature defaults to 0: no rand() involved) --- */
    StochasticLLG_Heun(&ctx, ctx.S, ctx.tS,
                        ctx.Heffx, ctx.Heffy, ctx.Heffz, ctx.RNx, ctx.RNy, ctx.RNz,
                        ctx.NOS, ctx.damping, ctx.t_step, ctx.Temperature,
                        0, 0, ctx.uABC[0], 0, ctx.uABC[1], 0, ctx.uABC[2]);
    print_spin("after_heun", 0, VEC_X(ctx.S,0), VEC_Y(ctx.S,0), VEC_Z(ctx.S,0));
    print_spin("after_heun", ctx.NOS / 2, VEC_X(ctx.S,ctx.NOS/2), VEC_Y(ctx.S,ctx.NOS/2), VEC_Z(ctx.S,ctx.NOS/2));
    print_spin("after_heun", ctx.NOS - 1, VEC_X(ctx.S,ctx.NOS-1), VEC_Y(ctx.S,ctx.NOS-1), VEC_Z(ctx.S,ctx.NOS-1));

    /* --- file round-trip: Save_VTK (binary) then Read_VTK into bS, compare --- */
    Save_VTK(&ctx, ctx.S, 0, "/tmp/magnoom_test_roundtrip.vtk");
    Read_VTK(&ctx, ctx.bS, "/tmp/magnoom_test_roundtrip.vtk");
    print_spin("roundtrip", 0, VEC_X(ctx.bS,0), VEC_Y(ctx.bS,0), VEC_Z(ctx.bS,0));
    print_spin("roundtrip", ctx.NOS / 2, VEC_X(ctx.bS,ctx.NOS/2), VEC_Y(ctx.bS,ctx.NOS/2), VEC_Z(ctx.bS,ctx.NOS/2));
    print_spin("roundtrip", ctx.NOS - 1, VEC_X(ctx.bS,ctx.NOS-1), VEC_Y(ctx.bS,ctx.NOS-1), VEC_Z(ctx.bS,ctx.NOS-1));
    remove("/tmp/magnoom_test_roundtrip.vtk");

    printf("TEST DONE\n");
    return 0;
}
