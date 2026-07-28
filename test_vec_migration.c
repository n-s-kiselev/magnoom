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
           ctx.Sx && ctx.Sy && ctx.Sz && ctx.tSx && ctx.tSy && ctx.tSz &&
           ctx.t2Sx && ctx.t2Sy && ctx.t2Sz && ctx.t3Sx && ctx.t3Sy && ctx.t3Sz &&
           ctx.bSx && ctx.bSy && ctx.bSz);

    /* --- geometry / kind / deterministic (HOMO) initial state --- */
    GetBox(&ctx, ctx.abc, ctx.uABC, ctx.Box);
    UpdateSpinPositions(&ctx, ctx.abc, ctx.uABC, ctx.Block, ctx.AtomsPerBlock, ctx.Box, ctx.Px, ctx.Py, ctx.Pz);
    UpdateKind(&ctx, ctx.Kind, ctx.Px, ctx.Py, ctx.Pz, ctx.NOS, ctx.NOSK);
    InitSpinComponents(&ctx, ctx.Px, ctx.Py, ctx.Pz, ctx.Sx, ctx.Sy, ctx.Sz, 1 /* HOMO, deterministic */);
    for (int i = 0; i < ctx.NOS; i++) { ctx.bSx[i] = ctx.Sx[i]; ctx.bSy[i] = ctx.Sy[i]; ctx.bSz[i] = ctx.Sz[i]; }

    printf("NOS=%d\n", ctx.NOS);
    print_spin("initial", 0, ctx.Sx[0], ctx.Sy[0], ctx.Sz[0]);
    print_spin("initial", ctx.NOS / 2, ctx.Sx[ctx.NOS/2], ctx.Sy[ctx.NOS/2], ctx.Sz[ctx.NOS/2]);
    print_spin("initial", ctx.NOS - 1, ctx.Sx[ctx.NOS-1], ctx.Sy[ctx.NOS-1], ctx.Sz[ctx.NOS-1]);

    /* --- copy S -> tS (memcpy-equivalent) --- */
    memcpy(ctx.tSx, ctx.Sx, ctx.NOS * sizeof(double));
    memcpy(ctx.tSy, ctx.Sy, ctx.NOS * sizeof(double));
    memcpy(ctx.tSz, ctx.Sz, ctx.NOS * sizeof(double));
    printf("copy S->tS matches: %d\n",
           memcmp(ctx.tSx, ctx.Sx, ctx.NOS*sizeof(double)) == 0 &&
           memcmp(ctx.tSy, ctx.Sy, ctx.NOS*sizeof(double)) == 0 &&
           memcmp(ctx.tSz, ctx.Sz, ctx.NOS*sizeof(double)) == 0);

    /* --- one deterministic integrator step (Temperature defaults to 0: no rand() involved) --- */
    StochasticLLG_Heun(&ctx, ctx.Sx, ctx.Sy, ctx.Sz, ctx.tSx, ctx.tSy, ctx.tSz,
                        ctx.Heffx, ctx.Heffy, ctx.Heffz, ctx.RNx, ctx.RNy, ctx.RNz,
                        ctx.NOS, ctx.damping, ctx.t_step, ctx.Temperature,
                        0, 0, ctx.uABC[0], 0, ctx.uABC[1], 0, ctx.uABC[2]);
    print_spin("after_heun", 0, ctx.Sx[0], ctx.Sy[0], ctx.Sz[0]);
    print_spin("after_heun", ctx.NOS / 2, ctx.Sx[ctx.NOS/2], ctx.Sy[ctx.NOS/2], ctx.Sz[ctx.NOS/2]);
    print_spin("after_heun", ctx.NOS - 1, ctx.Sx[ctx.NOS-1], ctx.Sy[ctx.NOS-1], ctx.Sz[ctx.NOS-1]);

    /* --- file round-trip: Save_VTK (binary) then Read_VTK into bS, compare --- */
    Save_VTK(&ctx, ctx.Sx, ctx.Sy, ctx.Sz, 0, "/tmp/magnoom_test_roundtrip.vtk");
    Read_VTK(&ctx, ctx.bSx, ctx.bSy, ctx.bSz, "/tmp/magnoom_test_roundtrip.vtk");
    print_spin("roundtrip", 0, ctx.bSx[0], ctx.bSy[0], ctx.bSz[0]);
    print_spin("roundtrip", ctx.NOS / 2, ctx.bSx[ctx.NOS/2], ctx.bSy[ctx.NOS/2], ctx.bSz[ctx.NOS/2]);
    print_spin("roundtrip", ctx.NOS - 1, ctx.bSx[ctx.NOS-1], ctx.bSy[ctx.NOS-1], ctx.bSz[ctx.NOS-1]);
    remove("/tmp/magnoom_test_roundtrip.vtk");

    printf("TEST DONE\n");
    return 0;
}
