// ////////////////////////////////////////////// GEOMETRY ////////////////////////////////////////////
// //translation vectors cubic:
// //float		abc[3][3] = {
// // 				{	1.0f, 0.0f, 0.0f }, // a
// // 				{	0.0f, 1.0f, 0.0f }, // b
// // 				{	0.0f, 0.0f, 1.0f }};// c
				
// //atom positions in basic domain: {x0,y0,z0, x1,y1,z1,...}
// //not more then 100 atoms per basic domain
// //FCC
// //float	Block[4*3] = 
// // 		{	0.0f,	0.0f,	0.0f,	
// // 			0.5f,	0.5f,	0.0f,
// // 			0.5f,	0.0f,	0.5f,
// // 			0.0f,	0.5f,	0.5f	
// // 		};
// //BCC
// /*
// float Block[2*3] = 
//     { 0.0f, 0.0f, 0.0f, 
//       0.5f, 0.5f, 0.5f, 
//     };
//  */




// //Simple Cubic 1 [001]
// // float		abc[3][3] = {
// // 				{	1.0f, 0.0f, 0.0f }, // a
// // 				{	0.0f, 1.0f, 0.0f }, // b
// // 				{	0.0f, 0.0f, 1.0f }};// c
// // float			Block[][3] = { 
// // 				{0.5f, 0.5f, 0.5f},  
// // 				};

// /*float		abc[3][3] = {
// 				{	3.0f, 0.0f, 0.0f }, // a
// 				{	0.0f, 1.0f, 0.0f }, // b
// 				{	0.0f, 0.0f, 10.0f }};// c
// float			Block[][3] = { 
// 				{0.0f, 0.5f, 0.0f}, 
// 				{1.5f, 0.5f, 0.0f}, 
// 				{3.0f, 0.5f, 0.0f}, 
// 				{0.75f, 0.0f, 0.0f},
// 				{2.25f, 0.0f, 0.0f},
// 				{0.75f, 1.0f, 0.0f},
// 				{2.25f, 1.0f, 0.0f}  
// 				};*/

// // float		abc[3][3] = {
// // 				{	1.0f, 0.0f, 0.0f }, // a
// // 				{	0.0f, 1.0f, 0.0f }, // b
// // 				{	0.0f, 0.0f, 1.0f }};// c
// // float			Block[][3] = { 
// // 				{0.5f, 0.5f, 0.5f},  
// // 				};
// // Simple Cubic 2 [011]
// // float		abc[3][3] = {
// // 				{	1.0f, 0.0f, 0.0f }, // a
// // 				{	0.0f, sqrt(2.f), 0.0f }, // b
// // 				{	0.0f, 0.0f, sqrt(2.f) }};// c				
// // float		Block[][3] = {		
// // 				{0,		      0,        0},	
// // 				{0.5,	sqrt(2.f)/2,	0}
// // 				};	

// //Simple Cubic 3 [111]
// // float		abc[3][3] = {
// // 				{	sqrt(2.f), 0.0f, 0.0f }, // a
// // 				{	0.0f, sqrt(3.f)*sqrt(2.f), 0.0f }, // b
// // 				{	0.0f, 0.0f, sqrt(3.f) }};// c				
// // float		Block[][3] = {		
// // 				{0,		      0,	        0},	
// // 				{sqrt(2.f)/2, sqrt(6.f)/2,	0},
// // 				{sqrt(2.f)/2, sqrt(6.f)/6,  sqrt(3.f)/3.f},
// // 				{0, 2*sqrt(6.f)/3,  sqrt(3.f)/3.f},
// // 				{0, 1*sqrt(6.f)/3, -sqrt(3.f)/3.f},
// // 				{sqrt(2.f)/2, 5*sqrt(6.f)/6, -sqrt(3.f)/3.f}
// // 				};

// //FCC 2
// // float		abc[3][3] = {
// // 				{	1.0f, 0.0f, 0.0f }, // a
// // 				{	0.0f, sqrt(2.f), 0.0f }, // b
// // 				{	0.0f, 0.0f, sqrt(2.f) }};// c				
// // float		Block[][3] = {		
// // 				{0,		      0,	          0},	
// // 				{0,	sqrt(2.f)/2,			  0},
// // 				{0.5,	sqrt(2.f)/4,sqrt(2.f)/4},
// // 				{0.5, 3*sqrt(2.f)/4, sqrt(2.f)/4},
// // 				{0,	sqrt(2.f)/2,	sqrt(2.f)/2}
// // 				};	
// //FCC 3
// // float		abc[3][3] = {
// // 				{sqrt(2.f), 0.0f, 0.0f }, // a
// // 				{sqrt(2.f)/2, sqrt(6.f)/2, 0.0f }, // b
// // 				{sqrt(2.f)/2, sqrt(6.f)/6, sqrt(3.f)/sqrt(2.f) }};// c				
// // float		Block[][3] = {		
// // 				{0,		      0,	          0},	
// // 				};


// //B20(1)
// // float		abc[3][3] = {
// // 				{	1.0f, 0.0f, 0.0f }, // a
// // 				{	0.0f, 1.0f, 0.0f }, // b
// // 				{	0.0f, 0.0f, 1.0f }};// c
// // #define uB20		0.138//0.138f//0.133
// // float		Block[][3] = {		
// // 				{   uB20,		   uB20,	    uB20},	
// // 				{0.5f+uB20,	0.5f-uB20,	1.0f-uB20},
// // 				{1.0f-uB20,	0.5f+uB20,	0.5f-uB20},	
// // 				{0.5f-uB20,	1.0f-uB20,	0.5f+uB20}	
// // 				};
// //B20(2)

// float		abc[3][3] = {
// 				{	1.0f, 0.0f, 0.0f }, // a
// 				{	0.0f, 1.0f, 0.0f }, // b
// 				{	0.0f, 0.0f, 1.0f }};// c

// float		Block[][3] = {		
// 				{   0.5,		   0.5,	    0.5}	
// 				};


// // #define uB20		0.138f
// // float		Block[][3] = {		
// // 				{   0.,		   0.,	    0.},	
// // 				{0.5f,	0.5f-2*uB20,	1.0f-2*uB20},
// // 				{1.0f-2*uB20,	0.5f,	0.5f-2*uB20},	
// // 				{0.5f-2*uB20,	1.0f-2*uB20,	0.5f}	
// // 				};

// //B20(Maria)
// // float		abc[3][3] = {
// // 				{	1.0f, 0.0f, 0.0f }, // a
// // 				{	0.0f, 1.0f, 0.0f }, // b
// // 				{	0.0f, 0.0f, 1.0f }};// c
// // //#define uB20		0.138f//MnSi
// // // float		Block[][3] = {		
// // // 				{       0.0f,        0.0f,        0.0f},//r1= u, u, u	
// // // 				{       0.5f, 0.5f-2*uB20,     -2*uB20},//r2= 0.5+u, 0.5-u,-u
// // // 				{    -2*uB20,        0.5f, 0.5f-2*uB20},//r3= -u, 0.5+u, 0.5 -u	
// // // 				{0.5f-2*uB20,     -2*uB20,        0.5f} //r4=0.5-u, -u, 0.5+u
// // // 				};
// // #define uB20		0.135f//FeGe
// // float		Block[][3] = {		
// // 				{       0.0f,        0.0f,        0.0f},//r1= u, u, u	
// // 				{       0.5f, 0.5f-2*uB20,     -2*uB20},//r2= 0.5+u, 0.5-u,-u
// // 				{    -2*uB20,        0.5f, 0.5f-2*uB20},//r3= -u, 0.5+u, 0.5 -u	
// // 				{0.5f-2*uB20,     -2*uB20,        0.5f} //r4=0.5-u, -u, 0.5+u
// // 				};

// //number of translations for the basic domain along a,b, and c verctors respectively 
// //int			uABC[3] = {2,147,2};//Grid dimensionality along translation vectors a, b, c; uABC[i]>0 
// int			uABC[3] = {10,10,10};//Grid dimensionality along translation vectors a, b, c; uABC[i]>0 
// //int			uABC[3] = {6,6,6};//Grid dimensionality along translation vectors a, b, c; uABC[i]>0 
// //int			uABC[3] = {80,80,20};//Grid dimensionality along translation vectors a, b, c; uABC[i]>0 


// int			ShellNumber = 1;
// int			AtomsPerBlock = sizeof(Block)/sizeof(float)/3;
// float*		RadiusOfShell = (float *)calloc(ShellNumber , sizeof(float));  
// int*		NeighborsPerAtom = (int *)calloc(AtomsPerBlock, sizeof(int));
// // total number of neighbour pairs per whole map of neighbours
// int			NOS=AtomsPerBlock*uABC[0]*uABC[1]*uABC[2]; // number of spins
// int			NOS_AL=AtomsPerBlock*uABC[1]*uABC[2]; // number of spins per A layer
// int			NOS_BL=AtomsPerBlock*uABC[0]*uABC[2]; // number of spins per B layer
// int			NOS_CL=AtomsPerBlock*uABC[0]*uABC[1]; // number of spins per C layer

// int 		NOSK = 0;

// double 		iNOS = 1.0/NOS;

// int			NOB=uABC[0]*uABC[1]*uABC[2]; // number of Blocks
// int			NOB_AL=uABC[1]*uABC[2]; // number of spins per A layer
// int			NOB_BL=uABC[0]*uABC[2]; // number of spins per B layer
// int			NOB_CL=uABC[0]*uABC[1]; // number of spins per C layer

/*****************************************************************************/
/* Coordination-shell construction: shells are distinguished by both         */
/* distance and local point symmetry, so the number of final shells is not   */
/* known before construction. These are generous compile-time caps for the   */
/* fixed/preallocated scratch buffers used while building the map (never     */
/* dynamically grown/reallocated).                                           */
/*****************************************************************************/
#define MAX_NEIGHBORS_PER_SHELL_PER_ATOM 256
#define MAX_RAW_SHELL_LEVELS 64
#define MAX_NEIGHBOR_PAIRS_TEMP (MAX_SHELL_NUM * MAX_ATOMS_PER_BLOCK * MAX_NEIGHBORS_PER_SHELL_PER_ATOM)
#define SHELL_SIGNATURE_TOL 1e-4f
/* Periodic images searched in either direction along each lattice vector. */
#define MAX_IMAGE_TRANSLATION 15

/* Set to 0 to restore the pre-migration behavior: shells grouped purely by
 * distance, with no splitting by local point symmetry (every bond at a
 * given raw distance shares one shell, regardless of angular signature or
 * which atom pair it connects). Exists to compare the two side by side --
 * e.g. for a low-symmetry structure like EuSi, where the realized shell
 * count/assignment can differ from the old distance-only grouping. Flip
 * and rebuild (`./nob`) to switch; only ClassifyCanonicalBonds() below
 * changes behavior, everything else in the neighbor-map pipeline is
 * identical either way. */
#ifndef MAGNOOM_SHELL_SYMMETRY_AWARE
#define MAGNOOM_SHELL_SYMMETRY_AWARE 0
#endif

/* One candidate bond from a central atom to a periodic image of another
 * (or the same) basis atom; dx/dy/dz is the Cartesian displacement
 * block[I0] - (block[I1] + J*a + K*b + L*c). */
typedef struct NeighborCandidate
{
	int   I1, J, K, L;
	float dx, dy, dz;
} NeighborCandidate;

/* One finished entry of the in-progress neighbor map: the bond from atom I0
 * in the home block to the image (I1, J, K, L), and the shell it belongs to. */
typedef struct NeighborPair
{
	int I0, I1, J, K, L, shell;
} NeighborPair;

/* One (angular signature, unordered atom-pair) class already assigned a
 * global shell index. Buckets live on ClassifyCanonicalBonds' stack, so they
 * start empty at every raw distance level -- angular comparisons are only
 * meaningful between bonds at the same distance. */
typedef struct ShellBucket
{
	float signature[MAX_NEIGHBORS_PER_SHELL_PER_ATOM];
	int   signatureLength;
	int   atomLo, atomHi;
	int   finalShell;
} ShellBucket;

/* Number of periodic images searched in either direction along lattice
 * vector `axis`: none at all along a degenerate (zero) vector. */
static int
ImageRange(float abc[][3], int axis)
{
	if (abc[axis][0]==0 && abc[axis][1]==0 && abc[axis][2]==0) return 0;
	return MAX_IMAGE_TRANSLATION;
}

/* Fills `bond` with the displacement from atom I0 in the home block to the
 * image (I1, J, K, L) of atom I1, and returns the length of that bond. */
static float
MeasureBond(float abc[][3], float block[][3], int I0, int I1, int J, int K, int L, NeighborCandidate *bond)
{
	float x1 = block[I1][0] + abc[0][0]*J + abc[1][0]*K + abc[2][0]*L;
	float y1 = block[I1][1] + abc[0][1]*J + abc[1][1]*K + abc[2][1]*L;
	float z1 = block[I1][2] + abc[0][2]*J + abc[1][2]*K + abc[2][2]*L;
	bond->I1 = I1;
	bond->J  = J;
	bond->K  = K;
	bond->L  = L;
	bond->dx = block[I0][0] - x1;
	bond->dy = block[I0][1] - y1;
	bond->dz = block[I0][2] - z1;
	return sqrtf(bond->dx*bond->dx + bond->dy*bond->dy + bond->dz*bond->dz);
}

static float
CosAngleBetweenBonds(const NeighborCandidate *a, const NeighborCandidate *b)
{
	float dot = a->dx*b->dx + a->dy*b->dy + a->dz*b->dz;
	float lenA = sqrtf(a->dx*a->dx + a->dy*a->dy + a->dz*a->dz);
	float lenB = sqrtf(b->dx*b->dx + b->dy*b->dy + b->dz*b->dz);
	return dot/(lenA*lenB);
}

static void
SortAscending(float *values, int count)
{
	for (int i = 1; i < count; i++) {
		float v = values[i];
		int j = i - 1;
		while (j >= 0 && values[j] > v) {
			values[j+1] = values[j];
			j--;
		}
		values[j+1] = v;
	}
}

/* Elementwise relative-tolerance comparison (tol scaled by float precision,
 * not shell.h's double_1arrs_equal's tol=1e-12 -- see migration notes). */
static bool
SignaturesEqual(const float *a, int aLen, const float *b, int bLen, float tol)
{
	if (aLen != bLen) return false;
	for (int i = 0; i < aLen; i++) {
		float diff = fabsf(a[i] - b[i]);
		float scale = fabsf(a[i]);
		if (fabsf(b[i]) > scale) scale = fabsf(b[i]);
		if (scale < 1.0f) scale = 1.0f;
		if (diff > tol*scale) return false;
	}
	return true;
}

/* Finds the smallest interatomic distance strictly greater than priorRadius,
 * searching every basic-domain atom against every periodic image within
 * ImageRange() unit-cell translations. Returns false if no larger distinct
 * distance exists. */
static bool
FindNextShellRadius(float abc[][3], float block[][3], int atomsPerBlock, float priorRadius, float *outRadius)
{
	int Jin = ImageRange(abc, 0), Kin = ImageRange(abc, 1), Lin = ImageRange(abc, 2);
	float smallestGap = 1e30f;
	bool found = false;
	for (int I0 = 0; I0 < atomsPerBlock; I0++)
	for (int J = -Jin; J <= Jin; J++)
	for (int K = -Kin; K <= Kin; K++)
	for (int L = -Lin; L <= Lin; L++)
	for (int I1 = 0; I1 < atomsPerBlock; I1++) {
		if (I0==I1 && J==0 && K==0 && L==0) continue; // not the same atom
		NeighborCandidate bond;
		float testRadius = MeasureBond(abc, block, I0, I1, J, K, L, &bond);
		float gap = testRadius - priorRadius;
		if (gap > 1e-6f && gap < smallestGap) {
			smallestGap = gap;
			*outRadius = testRadius;
			found = true;
		}
	}
	return found;
}

/* Gathers every candidate bond from atom I0 whose distance matches radius
 * (tolerance 1e-6, as elsewhere in this file). Returns false if there are
 * more than MAX_NEIGHBORS_PER_SHELL_PER_ATOM of them. */
static bool
GatherShellCandidates(float abc[][3], float block[][3], int atomsPerBlock, int I0, float radius,
	NeighborCandidate *candidates, int *outCandidateCount)
{
	int Jin = ImageRange(abc, 0), Kin = ImageRange(abc, 1), Lin = ImageRange(abc, 2);
	int candidateCount = 0;
	for (int J = -Jin; J <= Jin; J++)
	for (int K = -Kin; K <= Kin; K++)
	for (int L = -Lin; L <= Lin; L++)
	for (int I1 = 0; I1 < atomsPerBlock; I1++) {
		if (I0==I1 && J==0 && K==0 && L==0) continue; // not the same atom
		NeighborCandidate bond;
		float testRadius = MeasureBond(abc, block, I0, I1, J, K, L, &bond);
		if (fabsf(testRadius - radius) >= 1e-6f) continue; // atom I1 is not in this shell
		if (candidateCount == MAX_NEIGHBORS_PER_SHELL_PER_ATOM) return false;
		candidates[candidateCount] = bond;
		candidateCount++;
	}
	*outCandidateCount = candidateCount;
	return true;
}

/* Fixed total order over the two directions of a bond, picking exactly one
 * of (I0,I1,J,K,L) and its reciprocal (I1,I0,-J,-K,-L) as the "canonical"
 * one -- the direction that gets classified by its angular signature, the
 * other being derived from it (see pass 2 below). */
static bool
IsCanonicalDirection(int I0, int I1, int J, int K, int L)
{
	if (I1 != I0) return I1 > I0;
	if (J != 0) return J > 0;
	if (K != 0) return K > 0;
	return L > 0; // L != 0 guaranteed: the exact self-pair (0,0,0) is excluded
}

/* Appends one classified bond to the in-progress map. Returns false, loudly,
 * when the fixed-capacity scratch map is full. */
static bool
AppendPair(NeighborPair *pairs, int *pairCount, NeighborPair pair)
{
	if (*pairCount == MAX_NEIGHBOR_PAIRS_TEMP) {
		fprintf(stderr, "magnoom_ctx_build_neighbor_map: exceeded %d total neighbor pairs.\n", MAX_NEIGHBOR_PAIRS_TEMP);
		return false;
	}
	pairs[*pairCount] = pair;
	(*pairCount)++;
	return true;
}

/* Classifies the canonical direction of every bond at `radius` by its
 * angular signature and atom pair, appending each to `pairs`. Bonds whose
 * (signature, atom-pair) class was already seen AT THIS RADIUS reuse its
 * shell; every new class opens a new shell, so *shellCount grows by one per
 * class. This is where a raw distance level splits into several shells. */
static bool
ClassifyCanonicalBonds(magnoom_ctx *ctx, float radius, int *shellCount, NeighborPair *pairs, int *pairCount)
{
#if MAGNOOM_SHELL_SYMMETRY_AWARE
	ShellBucket buckets[MAX_SHELL_NUM];
	int bucketCount = 0;
#else
	// Distance-only mode: every bond at this radius shares one shell,
	// regardless of angular signature or which atom pair it connects --
	// matches the pre-migration GetShells()/CreateMapOfNeighbors() grouping.
	if (*shellCount == MAX_SHELL_NUM) {
		fprintf(stderr, "magnoom_ctx_build_neighbor_map: exceeded the maximum of %d shells.\n", MAX_SHELL_NUM);
		return false;
	}
	int levelShell = *shellCount;
	ctx->Neighbors.ShellRadius[levelShell] = radius;
	(*shellCount)++;
#endif

	for (int I0 = 0; I0 < ctx->AtomsPerBlock; I0++) {
		NeighborCandidate candidates[MAX_NEIGHBORS_PER_SHELL_PER_ATOM];
		int candidateCount;
		if (!GatherShellCandidates(ctx->abc, ctx->Block, ctx->AtomsPerBlock, I0, radius,
				candidates, &candidateCount)) {
			fprintf(stderr, "magnoom_ctx_build_neighbor_map: shell at radius %f around atom %d has more than %d neighbors.\n",
				(double)radius, I0, MAX_NEIGHBORS_PER_SHELL_PER_ATOM);
			return false;
		}

		for (int n = 0; n < candidateCount; n++) {
			const NeighborCandidate *bond = &candidates[n];
			if (!IsCanonicalDirection(I0, bond->I1, bond->J, bond->K, bond->L)) continue;

#if MAGNOOM_SHELL_SYMMETRY_AWARE
			// Rotation/reflection-invariant fingerprint of this bond: the
			// sorted cosines of its angle to every other bond in the same
			// shell around the same central atom (temp/src/shell.h's
			// symmetry_shell_signature, minus the constant self-cosine).
			float signature[MAX_NEIGHBORS_PER_SHELL_PER_ATOM];
			int sigLen = 0;
			for (int m = 0; m < candidateCount; m++) {
				if (m == n) continue;
				signature[sigLen++] = CosAngleBetweenBonds(bond, &candidates[m]);
			}
			SortAscending(signature, sigLen);

			// A canonical bond always runs from the lower to the higher
			// basis-atom index (see IsCanonicalDirection), so (I0, I1) is
			// already the ordered atom pair the bucket is keyed on.
			int matched = -1;
			for (int k = 0; k < bucketCount; k++) {
				if (buckets[k].atomLo == I0 && buckets[k].atomHi == bond->I1 &&
					SignaturesEqual(buckets[k].signature, buckets[k].signatureLength,
						signature, sigLen, SHELL_SIGNATURE_TOL)) {
					matched = k;
					break;
				}
			}
			if (matched < 0) {
				if (*shellCount == MAX_SHELL_NUM) {
					fprintf(stderr, "magnoom_ctx_build_neighbor_map: splitting the shell at radius %f exceeds the maximum of %d shells.\n",
						(double)radius, MAX_SHELL_NUM);
					return false;
				}
				memcpy(buckets[bucketCount].signature, signature, (size_t)sigLen*sizeof(float));
				buckets[bucketCount].signatureLength = sigLen;
				buckets[bucketCount].atomLo = I0;
				buckets[bucketCount].atomHi = bond->I1;
				buckets[bucketCount].finalShell = *shellCount;
				ctx->Neighbors.ShellRadius[*shellCount] = radius;
				matched = bucketCount;
				bucketCount++;
				(*shellCount)++;
			}

			NeighborPair canonical = { I0, bond->I1, bond->J, bond->K, bond->L, buckets[matched].finalShell };
#else
			NeighborPair canonical = { I0, bond->I1, bond->J, bond->K, bond->L, levelShell };
#endif
			if (!AppendPair(pairs, pairCount, canonical)) return false;
		}
	}
	return true;
}

/* Walks the distinct interatomic distances outwards, filling `pairs` with
 * both directions of every bond until targetShellCount shells exist. */
static bool
BuildShellLevels(magnoom_ctx *ctx, int targetShellCount, NeighborPair *pairs, int *outPairCount, int *outShellCount)
{
	int pairCount = 0;
	int shellCount = 0;
	float currRadius = 0.0f;

	for (int rawLevel = 0; shellCount < targetShellCount; rawLevel++) {
		if (rawLevel >= MAX_RAW_SHELL_LEVELS) {
			fprintf(stderr, "magnoom_ctx_build_neighbor_map: exceeded %d raw distance levels without reaching %d shells.\n",
				MAX_RAW_SHELL_LEVELS, targetShellCount);
			return false;
		}

		float nextRadius = 0.0f;
		if (!FindNextShellRadius(ctx->abc, ctx->Block, ctx->AtomsPerBlock, currRadius, &nextRadius)) {
			fprintf(stderr, "magnoom_ctx_build_neighbor_map: ran out of distinct interatomic distances before reaching %d shells.\n",
				targetShellCount);
			return false;
		}
		currRadius = nextRadius;

		int levelPairStart = pairCount;
		if (!ClassifyCanonicalBonds(ctx, currRadius, &shellCount, pairs, &pairCount)) return false;

		// The reciprocal of each bond just classified is the same physical
		// bond seen from the other end, so it belongs to the same shell.
		// Deriving the non-canonical directions from those records -- instead
		// of gathering and classifying them again -- is what makes that
		// guarantee exact: the two directions' per-atom angular signatures
		// need not agree when I0 and I1 sit in different local environments,
		// yet a symmetric pairwise exchange Hamiltonian needs one exchange
		// constant per bond, not one per direction.
		int levelCanonicalEnd = pairCount;
		for (int k = levelPairStart; k < levelCanonicalEnd; k++) {
			NeighborPair canonical = pairs[k];
			NeighborPair reciprocal = { canonical.I1, canonical.I0,
				-canonical.J, -canonical.K, -canonical.L, canonical.shell };
			if (!AppendPair(pairs, &pairCount, reciprocal)) return false;
		}
	}

	*outPairCount = pairCount;
	*outShellCount = shellCount;
	return true;
}

/* Builds ctx->Neighbors: the flattened template of neighbor-pair bonds in
 * the basic block, grouped into coordination shells distinguished by both
 * distance and local point symmetry. Bonds at the same distance are split
 * into separate shells whenever their rotation/reflection-invariant angular
 * signature (the sorted list of cosines between that bond and every other
 * bond in the same shell around the same central atom) or their unordered
 * (central, neighbor) atom-pair differs -- this ports the classification
 * logic of temp/src/shell.h's neighbor_map_new(), adapted to Magnoom's
 * global (not per-atom) shell index and to fixed/preallocated construction
 * instead of dynamic arrays.
 *
 * ctx->ShellNumber on entry is the desired number of final (post-split)
 * shells; on success it is overwritten with the realized count, which can
 * exceed the requested one whenever splitting occurs. */
bool
magnoom_ctx_build_neighbor_map(magnoom_ctx *ctx)
{
	if (ctx->ShellNumber <= 0 || ctx->ShellNumber > MAX_SHELL_NUM) {
		fprintf(stderr, "magnoom_ctx_build_neighbor_map: ShellNumber must be between 1 and %d.\n", MAX_SHELL_NUM);
		return false;
	}

	/* Generously sized scratch copy of the whole in-progress map, trimmed
	 * into the exact-size final arrays once the shell count is known. */
	NeighborPair *pairs = (NeighborPair *)malloc(MAX_NEIGHBOR_PAIRS_TEMP*sizeof(NeighborPair));
	if (pairs == NULL) {
		fprintf(stderr, "magnoom_ctx_build_neighbor_map: unable to allocate temporary neighbor-pair storage.\n");
		return false;
	}

	int pairCount = 0;
	int shellCount = 0;
	bool ok = BuildShellLevels(ctx, ctx->ShellNumber, pairs, &pairCount, &shellCount);

	if (ok) {
		ctx->Neighbors.AIdxBlock = (int *)malloc((size_t)pairCount*sizeof(int));
		ctx->Neighbors.NIdxBlock = (int *)malloc((size_t)pairCount*sizeof(int));
		ctx->Neighbors.NIdxGridA = (int *)malloc((size_t)pairCount*sizeof(int));
		ctx->Neighbors.NIdxGridB = (int *)malloc((size_t)pairCount*sizeof(int));
		ctx->Neighbors.NIdxGridC = (int *)malloc((size_t)pairCount*sizeof(int));
		ctx->Neighbors.ShellIdx  = (int *)malloc((size_t)pairCount*sizeof(int));
		if (ctx->Neighbors.AIdxBlock == NULL || ctx->Neighbors.NIdxBlock == NULL ||
			ctx->Neighbors.NIdxGridA == NULL || ctx->Neighbors.NIdxGridB == NULL ||
			ctx->Neighbors.NIdxGridC == NULL || ctx->Neighbors.ShellIdx == NULL) {
			fprintf(stderr, "magnoom_ctx_build_neighbor_map: unable to allocate the final neighbor-pair arrays.\n");
			free(ctx->Neighbors.AIdxBlock); free(ctx->Neighbors.NIdxBlock); free(ctx->Neighbors.NIdxGridA);
			free(ctx->Neighbors.NIdxGridB); free(ctx->Neighbors.NIdxGridC); free(ctx->Neighbors.ShellIdx);
			memset(&ctx->Neighbors, 0, sizeof(ctx->Neighbors));
			ok = false;
		} else {
			for (int i = 0; i < pairCount; i++) {
				ctx->Neighbors.AIdxBlock[i] = pairs[i].I0;
				ctx->Neighbors.NIdxBlock[i] = pairs[i].I1;
				ctx->Neighbors.NIdxGridA[i] = pairs[i].J;
				ctx->Neighbors.NIdxGridB[i] = pairs[i].K;
				ctx->Neighbors.NIdxGridC[i] = pairs[i].L;
				ctx->Neighbors.ShellIdx[i]  = pairs[i].shell;
			}
			ctx->Neighbors.PairCount = pairCount;
			ctx->Neighbors.ShellCount = shellCount;
			ctx->ShellNumber = shellCount;
		}
	}

	free(pairs);

	return ok;
}

void
GetPosition( float abc[][3], float block[][3], int I, int J, int K, int L, float * XYZ)
{
	XYZ[0]=block[I][0]+abc[0][0]*J+abc[1][0]*K+abc[2][0]*L;
	XYZ[1]=block[I][1]+abc[0][1]*J+abc[1][1]*K+abc[2][1]*L;
	XYZ[2]=block[I][2]+abc[0][2]*J+abc[1][2]*K+abc[2][2]*L;
}

void
SetExch1( float abc[][3], float block[][3], int arrSize, 
	float* Jijs, //input-> exchange coupling constant in each shell
	float* Bijs, //input-> bi-quadratic exchange coupling constant in each shell
	float* Dijs, //input-> modul of DMI vecotor in each shell
	int* aiBlock, int* niBlock, int* niGridA, int* niGridB, int* niGridC, int* sIdx, 
	float* J0, float* B0, float* D0, float* Dx, float* Dy, float* Dz)
{
	int I0,I1,J,K,L,S;
	double XYZ0[3], XYZ1[3];

    for(int i=0; i<arrSize; i++)
    {
    	I0= aiBlock[i];
    	I1= niBlock[i];
    	J = niGridA[i];
    	K = niGridB[i];
    	L = niGridC[i];
    	S = sIdx[i]; //shell index
    	//spin i position
		XYZ0[0]=block[I0][0];
		XYZ0[1]=block[I0][1];
		XYZ0[2]=block[I0][2];
		//spin j position
		XYZ1[0]=block[I1][0]+abc[0][0]*J+abc[1][0]*K+abc[2][0]*L;
		XYZ1[1]=block[I1][1]+abc[0][1]*J+abc[1][1]*K+abc[2][1]*L;
		XYZ1[2]=block[I1][2]+abc[0][2]*J+abc[1][2]*K+abc[2][2]*L;
		//r_{ij}
		XYZ1[0]-=XYZ0[0];//r1=r1-r0 /x
		XYZ1[1]-=XYZ0[1];//r1=r1-r0 /y
		XYZ1[2]-=XYZ0[2];//r1=r1-r0 /z	
		//metka
		// double norm[3]={0,0,1};
		// if (S==0) {
		// 	Cross(norm, XYZ1, XYZ0); (void)Unit( XYZ0, XYZ0);//interface induced DMI -- skyrmions of Neel type.
		// } else if (S==1){
		// 	(void)Unit( XYZ1, XYZ0);//Bloch skyrmion
		// }
		(void)Unit( XYZ1, XYZ0);//Bloch skyrmion
		Dx[i] = XYZ0[0];
		Dy[i] = XYZ0[1];
		Dz[i] = XYZ0[2];

		J0[i] = Jijs[S];
		D0[i] = Dijs[S];
		B0[i] = Bijs[S];
    }
}
