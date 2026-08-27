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

void anisotropy_rotate_site(magnoom_ctx *ctx, int atom)
{
	if (ctx == NULL || atom < 0 || atom >= ctx->AtomsPerBlock) return;
	anisotropy_rotate_tensor(&ctx->anisotropy_local[atom], ctx->anisotropy_rotation[atom],
		&ctx->anisotropy_global[atom]);
}

void anisotropy_rotate_all(magnoom_ctx *ctx)
{
	if (ctx == NULL) return;
	for (int atom = 0; atom < ctx->AtomsPerBlock; ++atom) anisotropy_rotate_site(ctx, atom);
}

int anisotropy_site_index(const magnoom_ctx *ctx, int atom)
{
	if (ctx == NULL || atom < 0 || atom >= ctx->AtomsPerBlock) return 0;
	return ctx->anisotropy_mode == ANISOTROPY_GLOBAL ? 0 : atom;
}

bool anisotropy_build_from_legacy(magnoom_ctx *ctx)
{
	if (ctx == NULL || ctx->AtomsPerBlock <= 0 ||
		!isfinite(ctx->Ku1) || !isfinite(ctx->Ku2) || !isfinite(ctx->Kc)) return false;
	for (int i = 0; i < 3; ++i) {
		if (!isfinite(ctx->VKu1[i]) || !isfinite(ctx->VKu2[i])) return false;
	}

	for (int atom = 0; atom < ctx->AtomsPerBlock; ++atom) {
		AnisotropyTensor *tensor = &ctx->anisotropy_local[atom];
		memset(tensor, 0, sizeof(*tensor));
		for (int i = 0; i < 3; ++i) {
			for (int j = i; j < 3; ++j) {
				double value = (double)ctx->Ku1*ctx->VKu1[i]*ctx->VKu1[j] +
					(double)ctx->Ku2*ctx->VKu2[i]*ctx->VKu2[j];
				if (!k2_set(tensor->K2, i, j, value)) return false;
			}
			if (!k4_set(tensor->K4, i, i, i, i, ctx->Kc)) return false;
		}
	}
	anisotropy_rotate_all(ctx);
	return true;
}
