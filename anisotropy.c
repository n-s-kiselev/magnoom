bool k2_set(double K2[3][3], int i, int j, double value)
{
	if (K2 == NULL || i < 0 || i >= 3 || j < 0 || j >= 3 || !isfinite(value)) return false;
	K2[i][j] = value;
	K2[j][i] = value;
	return true;
}

bool k4_set(double K4[3][3][3][3], int i, int j, int k, int l, double value)
{
	int index[4] = {i, j, k, l};
	if (K4 == NULL || i < 0 || i >= 3 || j < 0 || j >= 3 ||
		k < 0 || k >= 3 || l < 0 || l >= 3 || !isfinite(value)) return false;

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
	if (ctx == NULL || atom < 0 || atom >= ctx->AtomsPerBlock) return false;
	double angle = angle_degrees * (3.141592653589793 / 180.0);
	double q_new[4];
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
