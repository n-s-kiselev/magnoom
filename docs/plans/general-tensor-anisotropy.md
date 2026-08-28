# General tensor anisotropy

Status: Active
Owner: Shared
Last updated: 2026-08-28

## 1. Goal

Replace the global uniaxial and cubic anisotropy parameters with symmetric
rank-2 and rank-4 anisotropy tensors defined per atom in the unit cell.

The energy convention is:

    E = -(K2[i][j] m[i] m[j] + K4[i][j][k][l] m[i] m[j] m[k] m[l])

## 2. Motivation

The existing two-axis uniaxial model and separate cubic constant cannot
represent general crystal anisotropy. Tensor components provide one common
representation for uniaxial, cubic, tetragonal, and lower-symmetry cases.

## 3. Current State

The active physics now uses per-atom tensors. The legacy easy-axis vectors and
constants remain as inputs: `anisotropy_build_from_legacy()` converts them
into per-atom tensors at startup and whenever a legacy GUI control changes.
Effective field and energy routines use the rotated tensor cache.

## 4. Requirements

- [x] Store one local and one rotated tensor per active unit-cell atom.
- [x] Store one configurable local-to-global 3x3 rotation per unit-cell atom.
- [x] Enforce rank-2 and rank-4 symmetry in component setters.
- [x] Support atom index `-1` for assignments to every active atom.
- [x] Preserve existing anisotropy behavior through the physics migration.
- [x] Parse raw tensor and rotation records from `magnoom.cfg`.
- [ ] Generate GUI controls only for components nonzero after startup parsing.

## 5. Non-Goals

- Symmetry or point-group detection.
- Named cubic or tetragonal converters.
- Tensor validation UI.
- Magnetostatics or Epstein-zeta changes.
- Dynamic allocation for anisotropy tensors.
- New committed tests or benchmarks.

## 6. Constraints

- Use C99 and fixed storage for `MAX_ATOMS_PER_BLOCK` entries in `magnoom_ctx`.
- Only entries below `AtomsPerBlock` are active.
- Tensor component assignments must use symmetry-enforcing setters.
- Rotation from local to global coordinates happens outside integration loops.
- GUI tensor edits remain unsynchronized with solver threads, matching the
  existing direct parameter-control behavior.
- Each implementation stage is committed separately and reviewed before the
  next stage begins.

## 7. Relevant Files and Components

- `magnoom.c`: context, fixed capacities, initialization, and unity includes.
- `anisotropy.c`: tensor setters, contractions, gradients, and rotations.
- `solvers.c`: effective field and energy calculations.
- `visualization.c`: configuration parser and AntTweakBar controls.
- `magnoom.cfg`: startup tensor and rotation records.
- `README.md`: user-visible configuration documentation.

## 8. Existing Patterns to Reuse

- `AtomsPerBlock` is the site-type count and the basis-atom index is the site
  type selected in solver loops.
- Startup-generated exchange controls demonstrate variable-count widgets.
- `readConfigFile()` provides strict finite numeric parsing.
- Context initialization uses fixed-capacity inline arrays where practical.

## 9. Proposed Design

`AnisotropyTensor` contains `double K2[3][3]` and
`double K4[3][3][3][3]`. The context owns fixed-capacity arrays of local
tensors, global rotated tensors, and local-to-global matrices. The matrix
convention is `m_global[i] = R[i][a] m_local[a]`.

Anisotropy mode is runtime-selectable. Global mode makes every solver lookup
use atom 0's rotated tensor; Individual mode uses the basis-atom index. Mode
changes preserve all stored tensors. A GUI action can copy atom 0's local K2
and K4 tensors to every atom while preserving each atom's rotation matrix.

Configuration records use zero-based atom indices and one-based tensor/matrix
indices. Atom `-1` applies an assignment globally. Records are processed in
order, so later records override earlier assignments.

## 10. Implementation Steps

### Step 1: Data model

Purpose:
- Establish the tensor representation without changing behavior.

Changes:
- Add context storage and identity local frames.
- Add symmetric setters, atom-aware setters, contractions, gradients, and
  local-to-global tensor rotation.

Validation:
- Existing tests, macOS build, simulated Linux compile, and diff review.

Status:
- Completed

### Step 2: Physics wiring

Purpose:
- Replace legacy anisotropy formulas while preserving their results.

Changes:
- Convert both uniaxial terms and the cubic term to tensors.
- Use rotated per-atom tensors in field and energy routines.
- Add runtime Global/Individual selection semantics, defaulting to Global.
- Temporarily retain legacy GUI controls through tensor-rebuild callbacks.

Validation:
- Compare old and new formulas with an uncommitted numerical harness.
- Existing build and platform checks.

Status:
- Completed

### Step 3: Configuration

Purpose:
- Permit arbitrary local tensor and rotation components at startup.

Changes:
- Parse `K2`, `K4`, and `R` records with atom `-1` support.
- Validate ranges, finite values, and proper rotation matrices.
- Rotate tensors once after configuration and basis selection.

Validation:
- Exercise valid, global, overriding, and malformed records manually.
- Existing build and platform checks.

Status:
- Completed

### Step 4: Sparse GUI

Purpose:
- Expose only the configured nonzero independent components.

Changes:
- Remove legacy anisotropy controls.
- Add a dedicated iconified F6 Anisotropy bar.
- Add a runtime Global/Individual selector and an atom selector shown only in
  Individual mode.
- Create six reusable K2 and fifteen reusable K4 callback controls with a
  fixed `0.000001` increment and explicit one-based index labels.
- Redirect controls to atom 0 in Global mode or the selected atom in Individual
  mode, showing only components active in that atom's sparsity mask.
- Add a button that copies atom 0's local K2/K4 tensors to all atoms, preserves
  destination rotations, refreshes all rotated caches, and expands destination
  masks with atom 0's active components.
- Update all permutations and the rotated cache from control callbacks.

Validation:
- Inspect generated controls and live component updates manually.
- Existing build and platform checks.

Status:
- Pending

## 11. Validation Strategy

- Build: `./nob`
- Existing tests: `./nob -test`
- Platform check: simulated Linux C99 compilation.
- Regression check: temporary Stage 2 old/new energy and field comparison.
- Formatting: `git diff --check`.
- Manual checks: config parsing and sparse GUI after their respective stages.

## 12. Risks and Mitigations

- Rank-4 multiplicity errors: all assignments pass through a setter that writes
  every index permutation.
- Sign regression: retain a leading minus sign in energy so positive legacy
  constants map directly to positive tensor components.
- Per-atom indexing errors: use the existing basis-atom index in nested loops
  and `i % AtomsPerBlock` in flat loops.
- Rotation cost in hot loops: maintain a pre-rotated global tensor cache.
- GUI races: accepted explicitly to match existing live parameter controls.

## 13. Open Questions

None.

## 14. Progress Log

### 2026-08-27

- Completed Stage 0 discovery.
- Confirmed the leading-minus energy convention, identity default frames,
  fixed maximum capacity, raw config syntax, and unsynchronized GUI edits.
- Added fixed-capacity local/global tensors and identity local frames.
- Added symmetry-enforcing setters, atom `-1` assignment, energy and gradient
  contractions, and local-to-global tensor rotation without solver wiring.
- Replaced the legacy uniaxial and cubic formulas in `GetEffectiveField`,
  `GetTotalEnergyFerro`, and `GetTotalEnergy` with rotated tensor
  contractions. Added runtime Global/Individual selection, defaulting to
  Global, and kept legacy GUI controls through tensor-rebuild callbacks.
- Verified old/new energy and effective-field equality with a temporary
  uncommitted harness using nonzero `Ku1`, `Ku2`, `Kc`, and non-axis-aligned
  easy axes; the harness passed and was removed.

### 2026-08-28

- Added strict ordered parsing for raw `K2`, `K4`, and local-to-global `R`
  records, including atom `-1` assignments and one-based component indices.
- Applied staged records after basis selection, validated active atom ranges
  and proper rotation matrices, and refreshed the rotated tensor cache.
- Exercised global assignment, later override, valid rotation, malformed
  component index, and invalid rotation cases with a temporary harness.

## 15. Final Result

Pending.

## 16. Remaining Limitations

Pending.
