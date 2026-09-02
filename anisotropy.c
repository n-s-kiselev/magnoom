bool k2_set(double K2[3][3], int i, int j, double value)
{
	if (K2 == NULL || i < 0 || i >= 3 || j < 0 || j >= 3 || !isfinite(value)) {
		fprintf(stderr, "k2_set: invalid K2 component (%d,%d)=%g -- indices must be 0..2.\n", i, j, value);
		return false;
	}
	K2[i][j] = value;
	K2[j][i] = value;
	return true;
}

bool k4_set(double K4[3][3][3][3], int i, int j, int k, int l, double value)
{
	int index[4] = {i, j, k, l};
	if (K4 == NULL || i < 0 || i >= 3 || j < 0 || j >= 3 ||
		k < 0 || k >= 3 || l < 0 || l >= 3 || !isfinite(value)) {
		fprintf(stderr, "k4_set: invalid K4 component (%d,%d,%d,%d)=%g -- indices must be 0..2.\n",
			i, j, k, l, value);
		return false;
	}

	for (int a = 0; a < 4; ++a) {
		for (int b = 0; b < 4; ++b) {
			if (b == a) continue;
			for (int c = 0; c < 4; ++c) {
				if (c == a || c == b) continue;
				for (int d = 0; d < 4; ++d) {
					if (d == a || d == b || d == c) continue;
					K4[index[a]][index[b]][index[c]][index[d]] = value;
				}
			}
		}
	}
	return true;
}

static const int anisotropy_k2_components[ANISOTROPY_K2_COMPONENT_COUNT][2] = {
	{0, 0}, {0, 1}, {0, 2}, {1, 1}, {1, 2}, {2, 2}
};

static const int anisotropy_k4_components[ANISOTROPY_K4_COMPONENT_COUNT][4] = {
	{0, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 0, 2}, {0, 0, 1, 1}, {0, 0, 1, 2},
	{0, 0, 2, 2}, {0, 1, 1, 1}, {0, 1, 1, 2}, {0, 1, 2, 2}, {0, 2, 2, 2},
	{1, 1, 1, 1}, {1, 1, 1, 2}, {1, 1, 2, 2}, {1, 2, 2, 2}, {2, 2, 2, 2}
};

bool anisotropy_set_k2(magnoom_ctx *ctx, int atom, int i, int j, double value)
{
	if (ctx == NULL || ctx->AtomsPerBlock <= 0 || atom < -1 || atom >= ctx->AtomsPerBlock) return false;
	if (ctx->anisotropy_mode == ANISOTROPY_GLOBAL) {
		if (atom > 0) return true;
		return k2_set(ctx->anisotropy_local[0].K2, i, j, value);
	}
	if (atom >= 0) return k2_set(ctx->anisotropy_local[atom].K2, i, j, value);
	for (int site = 0; site < ctx->AtomsPerBlock; ++site) {
		if (!k2_set(ctx->anisotropy_local[site].K2, i, j, value)) return false;
	}
	return ctx->AtomsPerBlock > 0;
}

bool anisotropy_set_k4(magnoom_ctx *ctx, int atom, int i, int j, int k, int l, double value)
{
	if (ctx == NULL || ctx->AtomsPerBlock <= 0 || atom < -1 || atom >= ctx->AtomsPerBlock) return false;
	if (ctx->anisotropy_mode == ANISOTROPY_GLOBAL) {
		if (atom > 0) return true;
		return k4_set(ctx->anisotropy_local[0].K4, i, j, k, l, value);
	}
	if (atom >= 0) return k4_set(ctx->anisotropy_local[atom].K4, i, j, k, l, value);
	for (int site = 0; site < ctx->AtomsPerBlock; ++site) {
		if (!k4_set(ctx->anisotropy_local[site].K4, i, j, k, l, value)) return false;
	}
	return ctx->AtomsPerBlock > 0;
}

double anisotropy_energy(const AnisotropyTensor *tensor, const double m[3])
{
	double energy = 0.0;
	if (tensor == NULL || m == NULL) return 0.0;

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			energy -= tensor->K2[i][j]*m[i]*m[j];
			for (int k = 0; k < 3; ++k) {
				for (int l = 0; l < 3; ++l)
					energy -= tensor->K4[i][j][k][l]*m[i]*m[j]*m[k]*m[l];
			}
		}
	}
	return energy;
}

void anisotropy_energy_gradient(const AnisotropyTensor *tensor, const double m[3], double gradient[3])
{
	double spin[3];
	if (gradient == NULL) return;
	if (tensor == NULL || m == NULL) {
		gradient[0] = gradient[1] = gradient[2] = 0.0;
		return;
	}
	spin[0] = m[0]; spin[1] = m[1]; spin[2] = m[2];
	gradient[0] = gradient[1] = gradient[2] = 0.0;

	for (int p = 0; p < 3; ++p) {
		for (int j = 0; j < 3; ++j) {
			gradient[p] -= 2.0*tensor->K2[p][j]*spin[j];
			for (int k = 0; k < 3; ++k) {
				for (int l = 0; l < 3; ++l)
					gradient[p] -= 4.0*tensor->K4[p][j][k][l]*spin[j]*spin[k]*spin[l];
			}
		}
	}
}

void anisotropy_rotate_tensor(const AnisotropyTensor *local, const double rotation[3][3],
	AnisotropyTensor *global)
{
	AnisotropyTensor local_copy;
	if (local == NULL || rotation == NULL || global == NULL) return;
	if (local == global) {
		local_copy = *local;
		local = &local_copy;
	}
	memset(global, 0, sizeof(*global));

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			for (int a = 0; a < 3; ++a) {
				for (int b = 0; b < 3; ++b)
					global->K2[i][j] += rotation[i][a]*rotation[j][b]*local->K2[a][b];
			}
			for (int k = 0; k < 3; ++k) {
				for (int l = 0; l < 3; ++l) {
					for (int a = 0; a < 3; ++a) {
						for (int b = 0; b < 3; ++b) {
							for (int c = 0; c < 3; ++c) {
								for (int d = 0; d < 3; ++d) {
									global->K4[i][j][k][l] += rotation[i][a]*rotation[j][b]*
										rotation[k][c]*rotation[l][d]*local->K4[a][b][c][d];
								}
							}
						}
					}
				}
			}
		}
	}
}

bool anisotropy_normalize_quaternion(const double input[4], double output[4])
{
	if (input == NULL || output == NULL) return false;
	double scale = 0.0;
	for (int i = 0; i < 4; ++i) {
		if (!isfinite(input[i])) return false;
		double magnitude = fabs(input[i]);
		if (magnitude > scale) scale = magnitude;
	}
	if (scale == 0.0) return false;
	double norm_squared = 0.0;
	for (int i = 0; i < 4; ++i) {
		double scaled = input[i]/scale;
		norm_squared += scaled*scaled;
	}
	double inverse_norm = 1.0/sqrt(norm_squared);
	for (int i = 0; i < 4; ++i) output[i] = (input[i]/scale)*inverse_norm;
	return true;
}

void anisotropy_quaternion_to_matrix(const double quaternion[4], double rotation[3][3])
{
	double q[4];
	if (rotation == NULL) return;
	if (!anisotropy_normalize_quaternion(quaternion, q)) {
		memset(rotation, 0, 9*sizeof(double));
		rotation[0][0] = rotation[1][1] = rotation[2][2] = 1.0;
		return;
	}
	double xx = q[0]*q[0], yy = q[1]*q[1], zz = q[2]*q[2];
	double xy = q[0]*q[1], xz = q[0]*q[2], yz = q[1]*q[2];
	double wx = q[3]*q[0], wy = q[3]*q[1], wz = q[3]*q[2];
	rotation[0][0] = 1.0 - 2.0*(yy + zz);
	rotation[0][1] = 2.0*(xy - wz);
	rotation[0][2] = 2.0*(xz + wy);
	rotation[1][0] = 2.0*(xy + wz);
	rotation[1][1] = 1.0 - 2.0*(xx + zz);
	rotation[1][2] = 2.0*(yz - wx);
	rotation[2][0] = 2.0*(xz - wy);
	rotation[2][1] = 2.0*(yz + wx);
	rotation[2][2] = 1.0 - 2.0*(xx + yy);
}

void anisotropy_rotate_site(magnoom_ctx *ctx, int atom)
{
	if (ctx == NULL || atom < 0 || atom >= ctx->AtomsPerBlock) return;
	double rotation[3][3];
	anisotropy_quaternion_to_matrix(ctx->anisotropy_quaternion[atom], rotation);
	anisotropy_rotate_tensor(&ctx->anisotropy_local[atom], rotation,
		&ctx->anisotropy_global[atom]);
}

bool anisotropy_set_quaternion(magnoom_ctx *ctx, int atom, const double quaternion[4])
{
	double normalized[4];
	if (ctx == NULL || atom < 0 || atom >= ctx->AtomsPerBlock ||
		!anisotropy_normalize_quaternion(quaternion, normalized)) return false;
	memcpy(ctx->anisotropy_quaternion[atom], normalized, sizeof(normalized));
	anisotropy_rotate_site(ctx, atom);
	return true;
}

void anisotropy_quaternion_from_axis_angle(const double axis[3], double angle, double quaternion[4])
{
	if (axis == NULL || quaternion == NULL) return;
	double norm = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
	double sina2 = sin(0.5 * angle);
	if (norm <= 0.0) {
		quaternion[0] = quaternion[1] = quaternion[2] = 0.0;
		quaternion[3] = 1.0;
		return;
	}
	quaternion[0] = sina2 * axis[0] / norm;
	quaternion[1] = sina2 * axis[1] / norm;
	quaternion[2] = sina2 * axis[2] / norm;
	quaternion[3] = cos(0.5 * angle);
}

void anisotropy_quaternion_multiply(const double q1[4], const double q2[4], double out[4])
{
	if (q1 == NULL || q2 == NULL || out == NULL) return;
	double r[4];
	r[0] = q1[3]*q2[0] + q1[0]*q2[3] + q1[1]*q2[2] - q1[2]*q2[1];
	r[1] = q1[3]*q2[1] + q1[1]*q2[3] + q1[2]*q2[0] - q1[0]*q2[2];
	r[2] = q1[3]*q2[2] + q1[2]*q2[3] + q1[0]*q2[1] - q1[1]*q2[0];
	r[3] = q1[3]*q2[3] - (q1[0]*q2[0] + q1[1]*q2[1] + q1[2]*q2[2]);
	out[0] = r[0]; out[1] = r[1]; out[2] = r[2]; out[3] = r[3];
}

bool anisotropy_compose_axis_angle(magnoom_ctx *ctx, int atom, const double axis[3], double angle_degrees)
{
	if (ctx == NULL || atom < 0 || atom >= ctx->AtomsPerBlock || axis == NULL) return false;
	double angle = angle_degrees * (3.141592653589793 / 180.0);
	double q_new[4] = {0.0, 0.0, 0.0, 1.0};
	anisotropy_quaternion_from_axis_angle(axis, angle, q_new);
	double q_existing[4];
	memcpy(q_existing, ctx->anisotropy_quaternion[atom], sizeof(q_existing));
	double q_combined[4];
	anisotropy_quaternion_multiply(q_new, q_existing, q_combined);
	return anisotropy_set_quaternion(ctx, atom, q_combined);
}

bool anisotropy_reset_quaternion(magnoom_ctx *ctx, int atom)
{
	if (ctx == NULL || atom < 0 || atom >= ctx->AtomsPerBlock) return false;
	double identity[4] = {0.0, 0.0, 0.0, 1.0};
	return anisotropy_set_quaternion(ctx, atom, identity);
}

void anisotropy_rotate_all(magnoom_ctx *ctx)
{
	if (ctx == NULL) return;
	for (int atom = 0; atom < ctx->AtomsPerBlock; ++atom) anisotropy_rotate_site(ctx, atom);
}

bool anisotropy_copy_atom0_tensors(magnoom_ctx *ctx)
{
	if (ctx == NULL || ctx->AtomsPerBlock <= 0) return false;
	for (int atom = 1; atom < ctx->AtomsPerBlock; ++atom) {
		ctx->anisotropy_local[atom] = ctx->anisotropy_local[0];
	}
	anisotropy_rotate_all(ctx);
	return true;
}

bool anisotropy_apply_config_records(magnoom_ctx *ctx)
{
	if (ctx == NULL || ctx->AtomsPerBlock <= 0) return false;
	for (int record_index = 0; record_index < ctx->anisotropy_config_record_count; ++record_index) {
		const AnisotropyConfigRecord *record = &ctx->anisotropy_config_records[record_index];
		if (record->atom < -1 || record->atom >= ctx->AtomsPerBlock) {
			fprintf(stderr, "magnoom.cfg:%d targets atom %d, but the active basis has %d atoms.\n",
				record->line, record->atom, ctx->AtomsPerBlock);
			return false;
		}
		int first = record->atom < 0 ? 0 : record->atom;
		int last = record->atom < 0 ? ctx->AtomsPerBlock : record->atom + 1;
		for (int atom = first; atom < last; ++atom) {
			switch (record->kind) {
				case ANISOTROPY_RECORD_K2:
					if (!k2_set(ctx->anisotropy_local[atom].K2,
						record->index[0], record->index[1], record->value)) return false;
					break;
				case ANISOTROPY_RECORD_K4:
					if (!k4_set(ctx->anisotropy_local[atom].K4,
						record->index[0], record->index[1], record->index[2], record->index[3],
						record->value)) return false;
					break;
				case ANISOTROPY_RECORD_QUATERNION:
					ctx->anisotropy_quaternion[atom][record->index[0]] = record->value;
					break;
				default:
					return false;
			}
		}
	}

	for (int atom = 0; atom < ctx->AtomsPerBlock; ++atom) {
		double normalized[4];
		if (!anisotropy_normalize_quaternion(ctx->anisotropy_quaternion[atom], normalized)) {
			fprintf(stderr, "magnoom.cfg defines an invalid quaternion for atom %d.\n", atom);
			return false;
		}
		memcpy(ctx->anisotropy_quaternion[atom], normalized, sizeof(normalized));
	}
	anisotropy_rotate_all(ctx);
	return true;
}

int anisotropy_site_index(const magnoom_ctx *ctx, int atom)
{
	if (ctx == NULL || atom < 0 || atom >= ctx->AtomsPerBlock) return 0;
	return ctx->anisotropy_mode == ANISOTROPY_GLOBAL ? 0 : atom;
}

/*****************************************************************************/
/* Anisotropy energy map: samples anisotropy_energy() over a (theta, phi)   */
/* grid and exports it as CSV and as a legacy-VTK triangulated surface, for */
/* every atom and for the total (summed) energy. Triggered by the F6 bar's */
/* "Anisotropy map to file" button (visualization.c:CB_ExportAnisotropyMap).*/
/*****************************************************************************/

typedef struct AnisotropyEnergyGrid
{
	int    n_theta;   /* theta = 0..pi inclusive */
	int    n_phi;     /* phi = 0..2*pi, periodic (no duplicate endpoint) */
	double theta_step;
	double phi_step;
} AnisotropyEnergyGrid;

static void anisotropy_map_direction(const AnisotropyEnergyGrid *grid, int i, int j, double n[3])
{
	double theta = i*grid->theta_step;
	double phi = j*grid->phi_step;
	n[0] = sin(theta)*cos(phi);
	n[1] = sin(theta)*sin(phi);
	n[2] = cos(theta);
}

/* Fills E[grid->n_theta * grid->n_phi] (row-major, index i*n_phi+j) with
 * tensor's energy at every sampled direction, via anisotropy_energy() --
 * no duplicated energy math. */
static void anisotropy_sample_energy(const AnisotropyEnergyGrid *grid, const AnisotropyTensor *tensor,
	double *E)
{
	for (int i = 0; i < grid->n_theta; ++i) {
		for (int j = 0; j < grid->n_phi; ++j) {
			double n[3];
			anisotropy_map_direction(grid, i, j, n);
			E[i*grid->n_phi + j] = anisotropy_energy(tensor, n);
		}
	}
}

/* Min-max normalizes E[count] into En[count], so that En = 0 at the lowest
 * sampled energy and 1 at the highest (a flat map normalizes to all zeros).
 * Every consumer -- the CSV's En column and both VTK surface radii -- reads
 * this one array, so they can never disagree. */
static void anisotropy_normalize_energy(const double *E, size_t count, double *En)
{
	double e_min = HUGE_VAL, e_max = -HUGE_VAL;
	for (size_t k = 0; k < count; ++k) {
		if (E[k] < e_min) e_min = E[k];
		if (E[k] > e_max) e_max = E[k];
	}
	double denom = (e_max > e_min) ? (e_max - e_min) : 1.0;
	for (size_t k = 0; k < count; ++k) En[k] = (E[k] - e_min)/denom;
}

static bool anisotropy_write_energy_csv(const char *path, const AnisotropyEnergyGrid *grid,
	const double *E, const double *En)
{
	FILE *file = fopen(path, "w");
	if (file == NULL) {
		fprintf(stderr, "Cannot open output file '%s': %s\n", path, strerror(errno));
		return false;
	}
	fputs("theta,phi,nx,ny,nz,E,En\n", file);
	for (int i = 0; i < grid->n_theta; ++i) {
		for (int j = 0; j < grid->n_phi; ++j) {
			double n[3];
			anisotropy_map_direction(grid, i, j, n);
			int k = i*grid->n_phi + j;
			fprintf(file, "%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n",
				i*grid->theta_step, j*grid->phi_step, n[0], n[1], n[2], E[k], En[k]);
		}
	}
	fclose(file);
	return true;
}

/* Legacy-VTK binary integers are big-endian, like the floats/doubles
 * end_swap_float()/end_swap_double() already convert for Save_VTK(). */
static int32_t end_swap_int32(int32_t val)
{
	int32_t ret = 0;
	for (int i = 0; i < 4; ++i) ((char *)&ret)[3-i] = ((char *)&val)[i];
	return ret;
}

/* VTK point index of ring (1..n_theta-2), column j; j may be n_phi, one past
 * the last column, which wraps back to column 0 so that phi is periodic with
 * no duplicated seam point. Must stay in step with the order in which
 * anisotropy_write_energy_vtk() emits POINTS: north pole, then the interior
 * rings ring-major/phi-minor, then the south pole. */
static int anisotropy_vtk_ring_point(const AnisotropyEnergyGrid *grid, int ring, int j)
{
	return 1 + (ring - 1)*grid->n_phi + j%grid->n_phi;
}

static void anisotropy_write_vtk_int(FILE *file, bool binary, int value)
{
	if (binary) {
		int32_t swapped = end_swap_int32((int32_t)value);
		fwrite(&swapped, sizeof(swapped), 1, file);
	} else {
		fprintf(file, "%d ", value);
	}
}

static void anisotropy_write_vtk_triangle(FILE *file, bool binary, int a, int b, int c)
{
	anisotropy_write_vtk_int(file, binary, 3);
	anisotropy_write_vtk_int(file, binary, a);
	anisotropy_write_vtk_int(file, binary, b);
	anisotropy_write_vtk_int(file, binary, c);
	if (!binary) fputs("\n", file);
}

/* Writes grid sample (i,j) as one VTK point at radius En (or 1-En when
 * inverse) along its direction. Both output formats take their coordinates
 * from the same float triple, so ASCII and binary can never disagree in
 * value -- only in printed precision. */
static void anisotropy_write_vtk_point(FILE *file, bool binary, const AnisotropyEnergyGrid *grid,
	const double *En, bool inverse, int i, int j)
{
	double n[3];
	anisotropy_map_direction(grid, i, j, n);
	double en = En[i*grid->n_phi + j];
	double scale = inverse ? (1.0 - en) : en;
	float xyz[3] = { (float)(scale*n[0]), (float)(scale*n[1]), (float)(scale*n[2]) };
	if (binary) {
		for (int c = 0; c < 3; ++c) {
			float swapped = end_swap_float(xyz[c], sizeof(xyz[c]));
			fwrite(&swapped, sizeof(swapped), 1, file);
		}
	} else {
		fprintf(file, "%.6g %.6g %.6g\n", xyz[0], xyz[1], xyz[2]);
	}
}

/* Writes one En value as VTK point-data scalar, in the same ASCII/binary
 * shape as anisotropy_write_vtk_int() (space-separated for ASCII, so the
 * caller adds one trailing newline after the whole block). */
static void anisotropy_write_vtk_scalar(FILE *file, bool binary, double value)
{
	if (binary) {
		float swapped = end_swap_float((float)value, sizeof(float));
		fwrite(&swapped, sizeof(swapped), 1, file);
	} else {
		fprintf(file, "%.6g ", value);
	}
}

/* Writes the (theta,phi) grid as a closed, triangulated legacy-VTK
 * UNSTRUCTURED_GRID surface: point (i,j) sits at r*n(theta_i,phi_j), where
 * r = En (regular) or r = 1-En (inverse); En is the caller's already
 * normalized energy -- this function neither samples nor renormalizes
 * energy, it only lays out geometry from the array it's given. The two
 * poles (theta=0, theta=pi) are each a single vertex (every phi maps to the
 * same direction there) connected to their nearest ring by a triangle fan;
 * interior rings are connected by a 2-triangle split per phi-step, with phi
 * wrapped modulo n_phi -- a closed manifold with no duplicated or
 * degenerate cells (checked via Euler's formula V-E+F=2 during design). */
static bool anisotropy_write_energy_vtk(const char *path, VtkFormatMode format,
	const AnisotropyEnergyGrid *grid, const double *En, bool inverse, const char *title)
{
	bool binary = (format == LEGACY_BINARY_VTK);
	int n_theta = grid->n_theta, n_phi = grid->n_phi;
	int point_count = (n_theta - 2)*n_phi + 2;
	int triangle_count = n_phi*(2*n_theta - 4);

	FILE *file = fopen(path, "wb");
	if (file == NULL) {
		fprintf(stderr, "Cannot open output file '%s': %s\n", path, strerror(errno));
		return false;
	}

	fputs("# vtk DataFile Version 2.0\n", file);
	fprintf(file, "%s\n", title);
	fputs(binary ? "BINARY\n" : "ASCII\n", file);
	fputs("DATASET UNSTRUCTURED_GRID\n", file);

	/* Point order defines the indices anisotropy_vtk_ring_point() returns:
	 * point 0 is the north pole, then one point per interior-ring sample
	 * (ring-major, phi-minor), then the south pole last. */
	int north_pole = 0;
	int south_pole = point_count - 1;
	fprintf(file, "POINTS %d float\n", point_count);
	anisotropy_write_vtk_point(file, binary, grid, En, inverse, 0, 0);
	for (int ring = 1; ring <= n_theta - 2; ++ring) {
		for (int j = 0; j < n_phi; ++j) {
			anisotropy_write_vtk_point(file, binary, grid, En, inverse, ring, j);
		}
	}
	anisotropy_write_vtk_point(file, binary, grid, En, inverse, n_theta - 1, 0);
	if (binary) fputs("\n", file);

	fprintf(file, "CELLS %d %d\n", triangle_count, triangle_count*4);
	for (int j = 0; j < n_phi; ++j) {
		anisotropy_write_vtk_triangle(file, binary, north_pole,
			anisotropy_vtk_ring_point(grid, 1, j), anisotropy_vtk_ring_point(grid, 1, j+1));
	}
	for (int ring = 1; ring <= n_theta - 3; ++ring) {
		for (int j = 0; j < n_phi; ++j) {
			int a = anisotropy_vtk_ring_point(grid, ring, j);
			int b = anisotropy_vtk_ring_point(grid, ring, j+1);
			int c = anisotropy_vtk_ring_point(grid, ring+1, j);
			int d = anisotropy_vtk_ring_point(grid, ring+1, j+1);
			anisotropy_write_vtk_triangle(file, binary, a, b, c);
			anisotropy_write_vtk_triangle(file, binary, b, d, c);
		}
	}
	for (int j = 0; j < n_phi; ++j) {
		anisotropy_write_vtk_triangle(file, binary, south_pole,
			anisotropy_vtk_ring_point(grid, n_theta - 2, j+1), anisotropy_vtk_ring_point(grid, n_theta - 2, j));
	}
	if (binary) fputs("\n", file);

	fprintf(file, "CELL_TYPES %d\n", triangle_count);
	for (int t = 0; t < triangle_count; ++t) anisotropy_write_vtk_int(file, binary, 5 /* VTK_TRIANGLE */);
	if (!binary) fputs("\n", file);

	/* Per-vertex normalized energy, for ParaView to color the surface by
	 * (interpolated across triangle faces) -- matches the radius each point
	 * was placed at (En for the regular surface, 1-En for the inverse one),
	 * so the scalar and the geometry always agree. Same point order as
	 * POINTS above, so scalar k always describes point k. */
	fprintf(file, "POINT_DATA %d\n", point_count);
	fprintf(file, "SCALARS %s float 1\n", inverse ? "1-En" : "En");
	fputs("LOOKUP_TABLE default\n", file);
	anisotropy_write_vtk_scalar(file, binary, inverse ? 1.0 - En[0] : En[0]);
	for (int ring = 1; ring <= n_theta - 2; ++ring) {
		for (int j = 0; j < n_phi; ++j) {
			double en = En[ring*n_phi + j];
			anisotropy_write_vtk_scalar(file, binary, inverse ? 1.0 - en : en);
		}
	}
	{
		double en = En[(n_theta - 1)*n_phi];
		anisotropy_write_vtk_scalar(file, binary, inverse ? 1.0 - en : en);
	}
	if (!binary) fputs("\n", file);

	fclose(file);
	return true;
}

/* Writes one energy map's three files: "<base_name>.csv", the regular surface
 * "<base_name>.vtk", and the inverse surface "<base_name>_inverse.vtk". E and
 * En are the same sampled and normalized arrays for all three, so the CSV and
 * both surfaces always describe one identical energy map. */
static bool anisotropy_write_map_files(magnoom_ctx *ctx, const AnisotropyEnergyGrid *grid,
	const char *base_name, const char *title, const char *inverse_title,
	const double *E, const double *En)
{
	VtkFormatMode format = ctx->anisotropy_map_vtk_format;
	char name[64], path[MAGNOOM_PATH_CAPACITY];

	snprintf(name, sizeof(name), "%s.csv", base_name);
	if (!magnoom_resolve_output_path(ctx, name, path, sizeof(path))) return false;
	if (!anisotropy_write_energy_csv(path, grid, E, En)) return false;

	snprintf(name, sizeof(name), "%s.vtk", base_name);
	if (!magnoom_resolve_output_path(ctx, name, path, sizeof(path))) return false;
	if (!anisotropy_write_energy_vtk(path, format, grid, En, false, title)) return false;

	snprintf(name, sizeof(name), "%s_inverse.vtk", base_name);
	if (!magnoom_resolve_output_path(ctx, name, path, sizeof(path))) return false;
	if (!anisotropy_write_energy_vtk(path, format, grid, En, true, inverse_title)) return false;

	return true;
}

/* Recomputes and overwrites every anisotropy-map output file (CSV + regular
 * and inverse-energy VTK surfaces, per atom and for the total), using the
 * live ctx->anisotropy_global tensors (already rotated -- see
 * anisotropy_rotate_site()) and the atom each one actually uses in the
 * running simulation (anisotropy_site_index(), matching solvers.c). */
bool anisotropy_export_energy_maps(magnoom_ctx *ctx)
{
	if (ctx == NULL) return false;
	if (!isfinite(ctx->anisotropy_map_theta_step) || ctx->anisotropy_map_theta_step <= 0.0 ||
		!isfinite(ctx->anisotropy_map_phi_step) || ctx->anisotropy_map_phi_step <= 0.0) {
		fprintf(stderr, "Anisotropy map: theta and phi steps must be positive.\n");
		return false;
	}

	AnisotropyEnergyGrid grid;
	grid.theta_step = ctx->anisotropy_map_theta_step;
	grid.phi_step = ctx->anisotropy_map_phi_step;
	grid.n_theta = (int)llround(PI/grid.theta_step) + 1;
	grid.n_phi = (int)llround(2.0*PI/grid.phi_step);
	if (grid.n_theta < 3 || grid.n_phi < 3) {
		fprintf(stderr, "Anisotropy map: the theta/phi steps produce too coarse a grid (n_theta=%d, n_phi=%d).\n",
			grid.n_theta, grid.n_phi);
		return false;
	}

	size_t sample_count = (size_t)grid.n_theta*(size_t)grid.n_phi;
	double *E_total = (double *)calloc(sample_count, sizeof(double));
	double *E_atom = (double *)malloc(sample_count*sizeof(double));
	double *En = (double *)malloc(sample_count*sizeof(double));
	if (E_total == NULL || E_atom == NULL || En == NULL) {
		fprintf(stderr, "Anisotropy map: unable to allocate the sampling grid.\n");
		free(E_total); free(E_atom); free(En);
		return false;
	}

	bool ok = true;
	for (int atom = 0; ok && atom < ctx->AtomsPerBlock; ++atom) {
		const AnisotropyTensor *tensor = &ctx->anisotropy_global[anisotropy_site_index(ctx, atom)];
		anisotropy_sample_energy(&grid, tensor, E_atom);
		anisotropy_normalize_energy(E_atom, sample_count, En);
		for (size_t k = 0; k < sample_count; ++k) E_total[k] += E_atom[k];

		char base_name[64];
		snprintf(base_name, sizeof(base_name), "anisotropy_atom_%d", atom);
		ok = anisotropy_write_map_files(ctx, &grid, base_name, "Anisotropy energy map",
			"Inverse anisotropy energy map", E_atom, En);
	}

	if (ok) {
		anisotropy_normalize_energy(E_total, sample_count, En);
		ok = anisotropy_write_map_files(ctx, &grid, "anisotropy_total", "Total anisotropy energy map",
			"Inverse total anisotropy energy map", E_total, En);
	}

	free(E_total);
	free(E_atom);
	free(En);
	return ok;
}
