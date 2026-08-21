static void magnoom_worker_barrier(magnoom_ctx *ctx, int thread)
{
	if (thread == THREADS_NUMBER - 1) {
		sem_post(ctx->sem_in[0]);
		sem_wait(ctx->sem_in[thread]);
		sem_post(ctx->sem_out[0]);
		sem_wait(ctx->sem_out[thread]);
	} else {
		sem_wait(ctx->sem_in[thread]);
		sem_post(ctx->sem_in[thread + 1]);
		sem_wait(ctx->sem_out[thread]);
		sem_post(ctx->sem_out[thread + 1]);
	}
}

static int magnoom_first_nonfinite_spin(const magnoom_ctx *ctx)
{
	for (int i = 0; i < ctx->NOS; ++i) {
		if (ctx->Kind[i] != 0 &&
			(!isfinite(VEC_X(ctx->S, i)) || !isfinite(VEC_Y(ctx->S, i)) ||
			 !isfinite(VEC_Z(ctx->S, i)))) return i;
	}
	return -1;
}

void GetEffectiveField(	magnoom_ctx *ctx, const double* s,
					int naini, 	int nafin,
					int nbini, 	int nbfin,
					int ncini, 	int ncfin)
{
	const int numNeighbors = ctx->NeighborPairs;
	const int *aidxBlock = ctx->AIdxBlock;
	const int *nidxBlock = ctx->NIdxBlock;
	const int *nidxGridA = ctx->NIdxGridA;
	const int *nidxGridB = ctx->NIdxGridB;
	const int *nidxGridC = ctx->NIdxGridC;
	const int *shellIdx = ctx->SIdx;
	const float *Jij = ctx->Jij;
	const float *Bij = ctx->Bij;
	const float *Dij = ctx->Dij;
	const float *VDMx = ctx->VDMx;
	const float *VDMy = ctx->VDMy;
	const float *VDMz = ctx->VDMz;
	const float *vku1 = ctx->VKu1;
	const float ku1 = ctx->Ku1;
	const float *vku2 = ctx->VKu2;
	const float ku2 = ctx->Ku2;
	const float kc = ctx->Kc;
	const float *BextDCDirection = ctx->BextDCDirection;
	const float BextDCMagnitude = ctx->BextDCMagnitude;
	double *heffx = ctx->Heffx;
	double *heffy = ctx->Heffy;
	double *heffz = ctx->Heffz;
	double tmp0;
	int Ip, I, J, K, L;
	int S;
	float Je, Bq, dx, dy, dz, DM;
	int i,j;
	int na1, Na = ctx->uABC[0];
	int nb1, Nb = ctx->uABC[1];
	int nc1, Nc = ctx->uABC[2];
	//single spin interactions (or potentila terms): Zeeman and Anizotropy:
	for (int Ip=0; Ip<ctx->AtomsPerBlock; Ip++)
	{
		//for (int nc=0; nc<Nc; nc++)//nc(a,b)=neghbor in the direction of c(a,b)-vector
		for (int nc=ncini; nc<ncfin; nc++)//nc(a,b)=neghbor in the direction of c(a,b)-vector
		{
			nc1 = Na * Nb * nc;
			for (int nb=nbini; nb<nbfin; nb++)
			{
				nb1 = Na * nb;
				for (int na=naini; na<nafin; na++)
				{
					i = Ip + ctx->AtomsPerBlock * ( na + nb1 + nc1 );// index of spin "i"
					// External-field Zeeman contribution:
					heffx[i] = BextDCMagnitude*BextDCDirection[0]+ctx->BextACEnabled*ctx->BextACScalar*ctx->BextACDirection[0];
					heffy[i] = BextDCMagnitude*BextDCDirection[1]+ctx->BextACEnabled*ctx->BextACScalar*ctx->BextACDirection[1];
					heffz[i] = BextDCMagnitude*BextDCDirection[2]+ctx->BextACEnabled*ctx->BextACScalar*ctx->BextACDirection[2];
					if(nc==0||nc==Nc-1){
					heffx[i] /= 2;
					heffy[i] /= 2;
					heffz[i] /= 2;
					}

					//uniaxial anisotropy1:
					tmp0 = VEC_X(s,i)*vku1[0] + VEC_Y(s,i)*vku1[1] + VEC_Z(s,i)*vku1[2];
					heffx[i]+= 2 * ku1 * vku1[0] * tmp0;
					heffy[i]+= 2 * ku1 * vku1[1] * tmp0;
					heffz[i]+= 2 * ku1 * vku1[2] * tmp0;
					//uniaxial anisotropy2:
					tmp0 = VEC_X(s,i)*vku2[0] + VEC_Y(s,i)*vku2[1] + VEC_Z(s,i)*vku2[2];
					heffx[i]+= 2 * ku2 * vku2[0] * tmp0;
					heffy[i]+= 2 * ku2 * vku2[1] * tmp0;
					heffz[i]+= 2 * ku2 * vku2[2] * tmp0;
					//cubic anisotropy:
					heffx[i]+= 4 * kc * VEC_X(s,i)*VEC_X(s,i)*VEC_X(s,i);
					heffy[i]+= 4 * kc * VEC_Y(s,i)*VEC_Y(s,i)*VEC_Y(s,i);
					heffz[i]+= 4 * kc * VEC_Z(s,i)*VEC_Z(s,i)*VEC_Z(s,i);			
				}
			}
		}
	}
	// pairwise spin interactions
	int bc_a; // boundary condition along "a"
	int bc_b; // boundary condition along "b"
	int bc_c; // boundary condition along "c"
	int bc_f=1; // boundary condition factor 
	for (int ni=0; ni<numNeighbors; ni++)//over the whole pairs 
	{
		Ip= aidxBlock[ni];//index I^prime of spin in the block
		I = nidxBlock[ni];//index I of neghbor spin in the block
		J = nidxGridA[ni];//relative index J of the block along a-vector for neghbor I
		K = nidxGridB[ni];//relative index K of the block along b-vector for neghbor I
		L = nidxGridC[ni];//relative index L of the block along c-vector for neghbor I
		S = shellIdx[ni];
		Je= Jij[S];//Comapre to the expression for the total energy
		//Je= Jij[S]*Jexc[ni];
		Bq= Bij[S];//where there is a factor 1/2 which commes from the double summation
		DM= Dij[S];//Here one has not the factor 1/2 because Ei=Heff*Si
		dx= VDMx[ni];
		dy= VDMy[ni];
		dz= VDMz[ni];
		//for (int nc=0; nc<Nc; nc++)//nc(a,b)=neghbor in the direction of c(a,b)-vector
		for (int nc=ncini; nc<ncfin; nc++)//nc(a,b)=neghbor in the direction of c(a,b)-vector
		{
			bc_c = 1 - (1-ctx->Boundary[2])*(( (2*Nc) + (nc+L) )/Nc )%2; // boundary condition along "c"
			for (int nb=nbini; nb<nbfin; nb++)
			{
				bc_b = 1 - (1-ctx->Boundary[1])*(( (2*Nb) + (nb+K) )/Nb )%2; // boundary condition along "b"
				for (int na=naini; na<nafin; na++)
				{
					bc_a = 1 - (1-ctx->Boundary[0])*(( (2*Na) + (na+J) )/Na )%2; // boundary condition along "a"
					bc_f = bc_a*bc_b*bc_c;
					na1 = na;
					nb1 = Na * nb;
					nc1 = Na * Nb * nc;
					i = Ip + ctx->AtomsPerBlock * ( na1 + nb1 + nc1 );// index of spin "i"
					na1 = (Na + na + J)%Na;
					nb1 = Na * ((Nb + nb + K)%Nb);
					nc1 = Na * Nb * ((Nc + nc + L)%Nc);					
					j = I  + ctx->AtomsPerBlock * ( na1 + nb1 + nc1 );// index of neighbouring spin "j"
					// Symmetric Heisenberg exchange:
					// heffx[i]+= Je * VEC_X(s,j) * bc_f;
					// heffy[i]+= Je * VEC_Y(s,j) * bc_f;
					// heffz[i]+= Je * VEC_Z(s,j) * bc_f;
					// Bi-quadratic exchange:
					// tmp0 = (VEC_X(s,i)*VEC_X(s,j) + VEC_Y(s,i)*VEC_Y(s,j) + VEC_Z(s,i)*VEC_Z(s,j)) * bc_f;
					// heffx[i]+= 2 * Bq * VEC_X(s,j) * tmp0;
					// heffy[i]+= 2 * Bq * VEC_Y(s,j) * tmp0;
					// heffz[i]+= 2 * Bq * VEC_Z(s,j) * tmp0;
					// Dzyaloshinskii-Moriya interaction (antisymmetric exchange):
					// heffx[i]+= DM*(VEC_Y(s,j)*dz - VEC_Z(s,j)*dy) * bc_f;
					// heffy[i]+= DM*(VEC_Z(s,j)*dx - VEC_X(s,j)*dz) * bc_f;
					// heffz[i]+= DM*(VEC_X(s,j)*dy - VEC_Y(s,j)*dx) * bc_f;	

					// a bit shorter version:
					
					tmp0 = Bq * 2 * (VEC_X(s,i)*VEC_X(s,j) + VEC_Y(s,i)*VEC_Y(s,j) + VEC_Z(s,i)*VEC_Z(s,j));				
					heffx[i] = heffx[i] + (Je * VEC_X(s,j) +  VEC_X(s,j) * tmp0 + DM*(VEC_Y(s,j)*dz - VEC_Z(s,j)*dy) )* bc_f;
					heffy[i] = heffy[i] + (Je * VEC_Y(s,j) +  VEC_Y(s,j) * tmp0 + DM*(VEC_Z(s,j)*dx - VEC_X(s,j)*dz) )* bc_f;
					heffz[i] = heffz[i] + (Je * VEC_Z(s,j) +  VEC_Z(s,j) * tmp0 + DM*(VEC_X(s,j)*dy - VEC_Y(s,j)*dx) )* bc_f;
				}
			}
		}
	}
}

double GetTotalEnergyMoment(magnoom_ctx *ctx)
{
	const double *s = ctx->bS;
	const double *Hx = ctx->Heffx;
	const double *Hy = ctx->Heffy;
	const double *Hz = ctx->Heffz;
	double *Etot = ctx->Etot0;
	double *Mtot = ctx->outputMtotal;
	const int N = ctx->NOS;
	double tmp0 = 0;
	Mtot[0] = 0;
	Mtot[1] = 0;
	Mtot[2] = 0;
	for (int i=0; i<N; i++)
	{

		Etot[i] = 1.0-Hx[i]*VEC_X(s,i) - Hy[i]*VEC_Y(s,i) - Hz[i]*VEC_Z(s,i);
		// metka test stochastic LLG
		double vtemp[3];
		// opposite to the rotation of the vector of external field
		RotateVector(VEC_X(s,i),VEC_Y(s,i),VEC_Z(s,i),0,0,1,-ctx->BextDCPhi,vtemp); //rotate about y by theta of the external field
		RotateVector(vtemp[0],vtemp[1],vtemp[2],0,1,0,-ctx->BextDCTheta,vtemp); //rotate about z by phi of the external field
		Mtot[0] = Mtot[0] + vtemp[0];
		Mtot[1] = Mtot[1] + vtemp[1];
		Mtot[2] = Mtot[2] + vtemp[2];
		// Mtot[0] = Mtot[0] + sx[i];
		// Mtot[1] = Mtot[1] + sy[i];
		// Mtot[2] = Mtot[2] + sz[i];
		tmp0 = tmp0 + Etot[i];
	}
	return tmp0/N;
}

double GetTotalEnergyFerro(magnoom_ctx *ctx)
{
	double sx = ctx->BextDCDirection[0];
	double sy = ctx->BextDCDirection[1];
	double sz = ctx->BextDCDirection[2];
	const int numNeighbors = ctx->NeighborPairs;
	const int *aidxBlock = ctx->AIdxBlock;
	const int *nidxGridA = ctx->NIdxGridA;
	const int *nidxGridB = ctx->NIdxGridB;
	const int *nidxGridC = ctx->NIdxGridC;
	const int *shellIdx = ctx->SIdx;
	const float *Jij = ctx->Jij;
	const float *Bij = ctx->Bij;
	const float *Dij = ctx->Dij;
	const float *VDMx = ctx->VDMx;
	const float *VDMy = ctx->VDMy;
	const float *VDMz = ctx->VDMz;
	const float *vku1 = ctx->VKu1;
	const float ku1 = ctx->Ku1;
	const float *vku2 = ctx->VKu2;
	const float ku2 = ctx->Ku2;
	const float kc = ctx->Kc;
	const float *BextDCDirection = ctx->BextDCDirection;
	const float BextDCMagnitude = ctx->BextDCMagnitude;
	double *Etot = ctx->Etot;
	const int N = ctx->NOS;
	double tmp0=1.0/sqrt(sx*sx+sy*sy+sz*sz);
	sx = sx/tmp0;
	sy = sy/tmp0;
	sz = sz/tmp0;
	tmp0=0;
	//single spin interactions (or potentila terms): Zeeman and Anizotropy:
	for (int i=0; i<N; i++)
	{
	// External-field Zeeman contribution:
	//Etot[i] =-BextDCMagnitude*(BextDCDirection[0]*sx+BextDCDirection[1]*sy+BextDCDirection[2]*sz);
	Etot[i] =-BextDCMagnitude*((BextDCDirection[0]+ctx->BextACEnabled*ctx->BextACScalar*ctx->BextACDirection[0])*sx+(BextDCDirection[1]+ctx->BextACEnabled*ctx->BextACScalar*ctx->BextACDirection[1])*sy+(BextDCDirection[2]+ctx->BextACEnabled*ctx->BextACScalar*ctx->BextACDirection[2])*sz);
	//uniaxial anisotropy1:
	tmp0 = sx*vku1[0] + sy*vku1[1] + sz*vku1[2];
	Etot[i]-= ku1 * tmp0 * tmp0;
	//uniaxial anisotropy2:
	tmp0 = sx*vku2[0] + sy*vku2[1] + sz*vku2[2];
	Etot[i]-= ku2 * tmp0 * tmp0;
	//cubic anisotropy:
	Etot[i]-= kc * (sx*sx*sx*sx + sy*sy*sy*sy + sz*sz*sz*sz);	
	}
	// pairwise spin interactions
	int Ip, J, K, L;
	int S;
	double Je, Bq, dx, dy, dz, DM;
	int na1, Na = ctx->uABC[0];
	int nb1, Nb = ctx->uABC[1];
	int nc1, Nc = ctx->uABC[2];
	int bc_a; // boundary condition along "a"
	int bc_b; // boundary condition along "b"
	int bc_c; // boundary condition along "c"
	int bc_f=1; // boundary condition factor 
	for (int ni=0; ni<numNeighbors; ni++)
	{
		Ip= aidxBlock[ni];//index I^prime of spin in the block
		// I = nidxBlock[ni];//index I of neghbor spin in the block
		J = nidxGridA[ni];//relative index J of block along a-vector for neghbor I
		K = nidxGridB[ni];//relative index K of block along b-vector for neghbor I
		L = nidxGridC[ni];//relative index L of block along c-vector for neghbor I
		S = shellIdx[ni];
		Je= 0.5*Jij[S];//Note, here we have factor 1/2 which comes from bouble summation over all spins
		//Je= 0.5*Jij[S]*Jexc[ni];
		Bq= 0.5*Bij[S];//because Jij, Bij, and Dij are coupling constants per PAIR of spins
		DM= 0.5*Dij[S];//Note, factor 1/2 is mising in effective field expression!!!
		dx= VDMx[ni];// in general vector (dx, dy, dz) are assumed to be unit vector
		dy= VDMy[ni];// but not necessarily,
		dz= VDMz[ni];// and if so, then DM (= Dij) should be set to 1.0
		for (int nc=0; nc<Nc; nc++)
		{
			bc_c = 1 - (1-ctx->Boundary[2])*(( (2*Nc) + (nc+L) )/Nc )%2; // boundary condition along "c"
			for (int nb=0; nb<Nb; nb++)
			{
				bc_b = 1 - (1-ctx->Boundary[1])*(( (2*Nb) + (nb+K) )/Nb )%2; // boundary condition along "b"
				for (int na=0; na<Na; na++)
				{
					bc_a = 1 - (1-ctx->Boundary[0])*(( (2*Na) + (na+J) )/Na )%2; // boundary condition along "a"
					bc_f = bc_a*bc_b*bc_c;
					na1 = na;
					nb1 = Na * nb;
					nc1 = Na * Nb * nc;
					int i = Ip + ctx->AtomsPerBlock * ( na1 + nb1 + nc1 );// index of spin "i"
					// na1 = (Na + na + J)%Na;
					// nb1 = Na * ((Nb + nb + K)%Nb);
					// nc1 = Na * Nb * ((Nc + nc + L)%Nc);
					// j = I  + AtomsPerBlock * ( na1 + nb1 + nc1 );// index of neighbouring spin "j"
					tmp0 = sx*sx + sy*sy + sz*sz;
					Etot[i] = Etot[i] - bc_f*( Je * tmp0 + Bq * tmp0 * tmp0 + DM*(dx*(sy*sz - sz*sy) + dy*(sz*sx - sx*sz) + dz*(sx*sy - sy*sx)) );
				}
			}
		}
	}
	tmp0=0;
	for (int i=0;i<ctx->NOS;i++)
	{
		tmp0 = tmp0 + Etot[i];
	}
	return tmp0;
}

double
GetTotalEnergy(magnoom_ctx *ctx)
{
	const double *s = ctx->bS;
	const int numNeighbors = ctx->NeighborPairs;
	const int *aidxBlock = ctx->AIdxBlock;
	const int *nidxBlock = ctx->NIdxBlock;
	const int *nidxGridA = ctx->NIdxGridA;
	const int *nidxGridB = ctx->NIdxGridB;
	const int *nidxGridC = ctx->NIdxGridC;
	const int *shellIdx = ctx->SIdx;
	const float *Jij = ctx->Jij;
	const float *Bij = ctx->Bij;
	const float *Dij = ctx->Dij;
	const float *VDMx = ctx->VDMx;
	const float *VDMy = ctx->VDMy;
	const float *VDMz = ctx->VDMz;
	const float *vku1 = ctx->VKu1;
	const float ku1 = ctx->Ku1;
	const float *vku2 = ctx->VKu2;
	const float ku2 = ctx->Ku2;
	const float kc = ctx->Kc;
	const float *BextDCDirection = ctx->BextDCDirection;
	const float BextDCMagnitude = ctx->BextDCMagnitude;
	double *Etot = ctx->Etot;
	double *Mtot = ctx->Mtot;
	const int N = ctx->NOS;
	double tmp0=0;
	Mtot[0] = 0;
	Mtot[1] = 0;
	Mtot[2] = 0;
	//single spin interactions (or potentila terms): Zeeman and Anizotropy:
	for (int i=0; i<N; i++)
	{
	// External-field Zeeman contribution:
	Etot[i] =-BextDCMagnitude*((BextDCDirection[0]+ctx->BextACEnabled*ctx->BextACScalar*ctx->BextACDirection[0])*VEC_X(s,i)+(BextDCDirection[1]+ctx->BextACEnabled*ctx->BextACScalar*ctx->BextACDirection[1])*VEC_Y(s,i)+(BextDCDirection[2]+ctx->BextACEnabled*ctx->BextACScalar*ctx->BextACDirection[2])*VEC_Z(s,i));
	//uniaxial anisotropy:
	tmp0 = VEC_X(s,i)*vku1[0] + VEC_Y(s,i)*vku1[1] + VEC_Z(s,i)*vku1[2];
	Etot[i]-= ku1 * tmp0 * tmp0;
	//uniaxial anisotropy:
	tmp0 = VEC_X(s,i)*vku2[0] + VEC_Y(s,i)*vku2[1] + VEC_Z(s,i)*vku2[2];
	Etot[i]-= ku2 * tmp0 * tmp0;
	//cubic anisotropy:
	Etot[i]-= kc * (VEC_X(s,i)*VEC_X(s,i)*VEC_X(s,i)*VEC_X(s,i) + VEC_Y(s,i)*VEC_Y(s,i)*VEC_Y(s,i)*VEC_Y(s,i) + VEC_Z(s,i)*VEC_Z(s,i)*VEC_Z(s,i)*VEC_Z(s,i));	
	Mtot[0] = Mtot[0] + VEC_X(s,i);
	Mtot[1] = Mtot[1] + VEC_Y(s,i);
	Mtot[2] = Mtot[2] + VEC_Z(s,i);
	}

	// pairwise spin interactions
	int Ip, I, J, K, L;
	int S;
	float Je, Bq, dx, dy, dz, DM;
	int i,j;
	int na1, Na = ctx->uABC[0];
	int nb1, Nb = ctx->uABC[1];
	int nc1, Nc = ctx->uABC[2];
	int bc_a; // boundary condition along "a"
	int bc_b; // boundary condition along "b"
	int bc_c; // boundary condition along "c"
	int bc_f=1; // boundary condition factor 
	for (int ni=0; ni<numNeighbors; ni++)
	{
		Ip= aidxBlock[ni];//index I^prime of spin in the block
		I = nidxBlock[ni];//index I of neghbor spin in the block
		J = nidxGridA[ni];//relative index J of block along a-vector for neghbor I
		K = nidxGridB[ni];//relative index K of block along b-vector for neghbor I
		L = nidxGridC[ni];//relative index L of block along c-vector for neghbor I
		S = shellIdx[ni];
		Je= 0.5*Jij[S];//Note, here we have factor 1/2 which comes from bouble summation over all spins
		//Je= 0.5*Jij[S]*Jexc[ni];
		Bq= 0.5*Bij[S];//because Jij, Bij, and Dij are coupling constants per PAIR of spins
		DM= 0.5*Dij[S];//Note, factor 1/2 is mising in effective field expression!!!
		dx= VDMx[ni];// in general vector (dx, dy, dz) are assumed to be unit vector
		dy= VDMy[ni];// but not necessarily,
		dz= VDMz[ni];// and if so, then DM (= Dij) should be set to 1.0
		for (int nc=0; nc<Nc; nc++)
		{
			bc_c = 1 - (1-ctx->Boundary[2])*(( (2*Nc) + (nc+L) )/Nc )%2; // boundary condition along "c"
			for (int nb=0; nb<Nb; nb++)
			{
				bc_b = 1 - (1-ctx->Boundary[1])*(( (2*Nb) + (nb+K) )/Nb )%2; // boundary condition along "b"
				for (int na=0; na<Na; na++)
				{
					bc_a = 1 - (1-ctx->Boundary[0])*(( (2*Na) + (na+J) )/Na )%2; // boundary condition along "a"
					bc_f = bc_a*bc_b*bc_c;
					na1 = na;
					nb1 = Na * nb;
					nc1 = Na * Nb * nc;
					i = Ip + ctx->AtomsPerBlock * ( na1 + nb1 + nc1 );// index of spin "i"
					na1 = (Na + na + J)%Na;
					nb1 = Na * ((Nb + nb + K)%Nb);
					nc1 = Na * Nb * ((Nc + nc + L)%Nc);;
					j = I  + ctx->AtomsPerBlock * ( na1 + nb1 + nc1 );// index of neighbouring spin "j"
					//Symmetric Heisenberg exchange:
					// Etot[i]-= Je * (VEC_X(s,i)*VEC_X(s,j)+VEC_Y(s,i)*VEC_Y(s,j)+VEC_Z(s,i)*VEC_Z(s,j));
					// //bi-quadratic exchange:
					// tmp0 = VEC_X(s,i)*VEC_X(s,j) + VEC_Y(s,i)*VEC_Y(s,j) + VEC_Z(s,i)*VEC_Z(s,j);
					// Etot[i]-= Bq * tmp0 * tmp0;
					// //Dzyaloshinskii-Moriya interaction (antisymmetric exchange):
					// Etot[i]-= DM*(dx*(VEC_Y(s,i)*VEC_Z(s,j) - VEC_Z(s,i)*VEC_Y(s,j)) + dy*(VEC_Z(s,i)*VEC_X(s,j) - VEC_X(s,i)*VEC_Z(s,j)) + dz*(VEC_X(s,i)*VEC_Y(s,j) - VEC_Y(s,i)*VEC_X(s,j)));
					tmp0 = VEC_X(s,i)*VEC_X(s,j) + VEC_Y(s,i)*VEC_Y(s,j) + VEC_Z(s,i)*VEC_Z(s,j);
					Etot[i] = Etot[i] - bc_f*( Je * tmp0 + Bq * tmp0 * tmp0 + DM*(dx*(VEC_Y(s,i)*VEC_Z(s,j) - VEC_Z(s,i)*VEC_Y(s,j)) + dy*(VEC_Z(s,i)*VEC_X(s,j) - VEC_X(s,i)*VEC_Z(s,j)) + dz*(VEC_X(s,i)*VEC_Y(s,j) - VEC_Y(s,i)*VEC_X(s,j))) );
				}
			}
		}
	}
	tmp0=0;
	for (i=0;i<ctx->NOS;i++)
	{
		tmp0 = tmp0 + Etot[i];
	}
	return tmp0;
}

// void
// GetFluctuations( float* rx, float* ry, float* rz, int N ){
// 	for (int i=0; i<N; i++){
// 		float s = 2;
// 		float x,y,z;
// 		while(s>1){
// 			x = 2.0 * (0.5 - rand() / (float)RAND_MAX);
// 			y = 2.0 * (0.5 - rand() / (float)RAND_MAX);
// 			z = 2.0 * (0.5 - rand() / (float)RAND_MAX);
// 			s = x*x+y*y+z*z;
// 		}
// 		rx[i] = x;
// 		ry[i] = y;
// 		rz[i] = z;		
// 	}
// }

void
GetFluctuations(magnoom_ctx *ctx){
	float *rx = ctx->RNx;
	float *ry = ctx->RNy;
	float *rz = ctx->RNz;
	const int N = ctx->NOS;
	for (int i=0; i<N; i++){
		float r[3];
		float A=sqrt(2*fabs(log(ctx->t_step)));
		for(int j=0; j<3; j++){
			float U = rand() / (float)RAND_MAX;
			float V = rand() / (float)RAND_MAX;
			r[j] = sqrt(-2*log(U))*cos(2*PI*V);	
			if(r[j]>A){
				r[j]=A;
			}else if(r[j]<-A){
				r[j]=-A;
			}	
		}	
		rx[i] = r[0];
 		ry[i] = r[1];
		rz[i] = r[2];
	}
}

void
StochasticLLG(	magnoom_ctx *ctx, int thread,
				int naini, 	int nafin,
				int nbini, 	int nbfin,
				int ncini, 	int ncfin)
// The semi-implicit midpoint (SIB) solver for the LLG-equations:
// J.H. Mentink et al, J. Phys.: Condens. Matter, 22, 176001 (2010).
// The method consists of two steps: prediction step and final step, see Eq.(18).
// At each step one has to find component of vector n1 solwing the equation: 
// n_1x - (n_1y*Az - n_1z*A_y)   n_0x + (n_0y*Az - n_0z*A_y)    pay attention  
// n_1x - (n_1y*Az - n_1z*A_y) = n_0x + (n_0y*Az - n_0z*A_y) <--to the sign in
// n_1x - (n_1y*Az - n_1z*A_y)   n_0x + (n_0y*Az - n_0z*A_y)    front of brackets
// where n0 is the unit vectors on the previouce step and
// A = h/2 * [-Heff(n) - alpha * (n0 cross Heff(n))]+sqrt(h)/2 * [-R - alpha * (n0 cross R)].
// In matrix form: 
//               M*n_1 = Mt*n_0,  (***)
// where Mt=transpose(M)
//        | 1  -Az  Ay |        | 1   Az -Ay |
//    M = | Az  1  -Ax |,  Mt = |-Az  1   Ax |
//        |-Ay  Ax  1  |        | Ay -Ax  1  |
// Let's introduce vector a = Mt*n_0, then 
//               M*n_1 = a. 
// The solution for n_1x, n_1y, n_1z in (***) 
// can be found with Cramers' rule: A*x=B, x_i = A_i/detA
// where  A_i  is the matrix formed by replacing the i-th column of A by the column vector B.
// Note: determinant M, detM = 1 + Ax^2 + Ay^2 + Az^2.
{
	double *in = ctx->S;
	double *tn = ctx->tS;
	const float *rx = ctx->RNx;
	const float *ry = ctx->RNy;
	const float *rz = ctx->RNz;
	const float alpha = ctx->damping;
	const float h = ctx->t_step;
	const float temperature = ctx->Temperature;
	double rh=sqrt(h);// h = \delta t
	double nx, ny, nz;// components of the unit vector
	double Hx, Hy, Hz;// components of the effective field
	float Cx, Cy, Cz;// spin-torque term
	float Rx, Ry, Rz;// random variables - thermal fluctuations
	double ax, ay, az;// temp
	double Ax, Ay, Az;// total matrix
	double detMi;	 // detMi = 1/detM
	double D = sqrt(2.0 * alpha / (1.0 + alpha * alpha) * temperature);
	double Alpha_d = alpha;// / (1.0 + alpha * alpha * Precession);
	double Alpha_p = 1.0f;// / (1.0 + alpha * alpha);

	//electric DC current vector (VCu) and density (Cu)
	Cx = ctx->VCu[0] * ctx->Cu;
	Cy = ctx->VCu[1] * ctx->Cu;
	Cz = ctx->VCu[2] * ctx->Cu;
	GetEffectiveField(ctx, in, naini, nafin, nbini, nbfin, ncini, ncfin);
	//prediction step of midpoint solver:
	int Na = ctx->uABC[0];
	int Nb = ctx->uABC[1];
	int nb1, nc1;
	int i;
	for (int Ip=0; Ip<ctx->AtomsPerBlock; Ip++)
	{
		//for (int nc=0; nc<Nc; nc++)//nc(a,b)=neghbor in the direction of c(a,b)-vector
		for (int nc=ncini; nc<ncfin; nc++)//nc(a,b)=neghbor in the direction of c(a,b)-vector
		{
			nc1 = Na * Nb * nc;
			for (int nb=nbini; nb<nbfin; nb++)
			{
				nb1 = Na * nb;
				for (int na=naini; na<nafin; na++)
				{
					i = Ip + ctx->AtomsPerBlock * ( na + nb1 + nc1 );// index of spin "i"
					nx = VEC_X(in,i);	ny = VEC_Y(in,i);	nz = VEC_Z(in,i);
					Hx = ctx->Heffx[i];	Hy = ctx->Heffy[i];	Hz = ctx->Heffz[i];
					Rx = rx[i];		Ry = ry[i];		Rz = rz[i];

					// deterministic terms of Landau–Lifshitz equation:
					Ax = 0.5f * h * ( - Alpha_p * ctx->Precession * Hx - Alpha_d * (ny * Hz - nz * Hy) );
					Ay = 0.5f * h * ( - Alpha_p * ctx->Precession * Hy - Alpha_d * (nz * Hx - nx * Hz) );
					Az = 0.5f * h * ( - Alpha_p * ctx->Precession * Hz - Alpha_d * (nx * Hy - ny * Hx) );

					// Spin-torque term
					if (ctx->Cu!=0) {
						Ax = Ax + 0.5f * h * ( - Alpha_d * Cx + (ny * Cz - nz * Cy) ); //pay attention to the signe and factors
						Ay = Ay + 0.5f * h * ( - Alpha_d * Cy + (nz * Cx - nx * Cz) );
						Az = Az + 0.5f * h * ( - Alpha_d * Cz + (nx * Cy - ny * Cx) );
					}

					// stochastic terms of Landau–Lifshitz equation:
					if (temperature>0) {
						Ax = Ax + 0.5f * rh * D * ( - Alpha_p * Rx - Alpha_d * (ny * Rz - nz * Ry) );
						Ay = Ay + 0.5f * rh * D * ( - Alpha_p * Ry - Alpha_d * (nz * Rx - nx * Rz) );
						Az = Az + 0.5f * rh * D * ( - Alpha_p * Rz - Alpha_d * (nx * Ry - ny * Rx) );
					}

					ax = nx + ny * Az - nz * Ay;
					ay = ny + nz * Ax - nx * Az;
					az = nz + nx * Ay - ny * Ax;

					// detMi  = 1.0f/(1.0f + Ax * Ax + Ay * Ay + Az * Az);
					// VEC_X(tn,i) = (ax*(1+Ax*Ax)+ay*(Ax*Ay+Az)+az*(Ax*Az-Ay))*detMi
					// VEC_Y(tn,i) = (ax*(Ay*Ax-Az)+ay*(1+Ay*Ay)+az*(Ay*Az+Ax))*detMi
					// VEC_Z(tn,i) = (ax*(Az*Ax+Ay)+ay*(Az*Ay-Ax)+az*(1+Az*Az))*detMi
					// VEC_X(tn,i) = ( VEC_X(tn,i) + nx ) / 2;
					// VEC_Y(tn,i) = ( VEC_X(tn,i) + nx ) / 2;
					// VEC_Z(tn,i) = ( VEC_X(tn,i) + nx ) / 2;

					// let's do it in a a bit more efficient way by using Hx, Hy, Hz, Rx, Ry and Rz as temporal variables:
					Hx = Ax * Ax;
					Hy = Ay * Ay;
					Hz = Az * Az;
					Rx = Ay * Az;
					Ry = Ax * Az;
					Rz = Ax * Ay;
					
					detMi = 1.0f/(1.0f + Hx + Hy + Hz);
					
					nx = nx + ( ax * (1. + Hx) + ay * (Rz + Az) + az * (Ry - Ay) ) * detMi;
					ny = ny + ( ax * (Rz - Az) + ay * (1. + Hy) + az * (Rx + Ax) ) * detMi;
					nz = nz + ( ax * (Ry + Ay) + ay * (Rx - Ax) + az * (1. + Hz) ) * detMi;

					VEC_X(tn,i) = nx * 0.5f;
					VEC_Y(tn,i) = ny * 0.5f;
					VEC_Z(tn,i) = nz * 0.5f;
				}
			}
		}
	}

	magnoom_worker_barrier(ctx, thread);
	GetEffectiveField(ctx, tn, naini, nafin, nbini, nbfin, ncini, ncfin);

	//final step of midpoint solver:
	for (int Ip=0; Ip<ctx->AtomsPerBlock; Ip++)
	{
		for (int nc=ncini; nc<ncfin; nc++)//nc(a,b)=neghbor in the direction of c(a,b)-vector
		{
			nc1 = Na * Nb * nc;
			for (int nb=nbini; nb<nbfin; nb++)
			{
				nb1 = Na * nb;
				for (int na=naini; na<nafin; na++)
				{
					i = Ip + ctx->AtomsPerBlock * ( na + nb1 + nc1 );// index of spin "i"
					nx = VEC_X(tn,i);	ny = VEC_Y(tn,i);	nz = VEC_Z(tn,i);	// <-- compare this line to corresponding one in prediction step
					Hx = ctx->Heffx[i];	Hy = ctx->Heffy[i];	Hz = ctx->Heffz[i];	// <-- they are new values for heff
					Rx = rx[i];		Ry = ry[i];		Rz = rz[i];		// <-- they are the same values as in prediction step
					// deterministic terms of Landau–Lifshitz equation:
					Ax = 0.5f * h * ( - Alpha_p * ctx->Precession * Hx - Alpha_d * (ny * Hz - nz * Hy) );
					Ay = 0.5f * h * ( - Alpha_p * ctx->Precession * Hy - Alpha_d * (nz * Hx - nx * Hz) );
					Az = 0.5f * h * ( - Alpha_p * ctx->Precession * Hz - Alpha_d * (nx * Hy - ny * Hx) );
					// Spin-torque term
					Ax = Ax + 0.5f * h * ( - Alpha_d * Cx + (ny * Cz - nz * Cy) ); //pay attention to the signe and factors
					Ay = Ay + 0.5f * h * ( - Alpha_d * Cy + (nz * Cx - nx * Cz) );
					Az = Az + 0.5f * h * ( - Alpha_d * Cz + (nx * Cy - ny * Cx) );
					// stochastic terms of Landau–Lifshitz equation:
					Ax = Ax + 0.5f * rh * D * ( - Alpha_p * Rx - Alpha_d * (ny * Rz - nz * Ry) );
					Ay = Ay + 0.5f * rh * D * ( - Alpha_p * Ry - Alpha_d * (nz * Rx - nx * Rz) );
					Az = Az + 0.5f * rh * D * ( - Alpha_p * Rz - Alpha_d * (nx * Ry - ny * Rx) );

					nx = VEC_X(in,i);	ny = VEC_Y(in,i);	nz = VEC_Z(in,i);	//<-- pay attention to this line!

					ax = nx + ny * Az - nz * Ay;
					ay = ny + nz * Ax - nx * Az;
					az = nz + nx * Ay - ny * Ax;

					Hx = Ax * Ax;
					Hy = Ay * Ay;
					Hz = Az * Az;
					Rx = Ay * Az;
					Ry = Ax * Az;
					Rz = Ax * Ay;
					
					detMi = 1.0f/(1.0f + Hx + Hy + Hz);
					
					VEC_X(in,i) = ( ax * (1. + Hx) + ay * (Rz + Az) + az * (Ry - Ay) ) * detMi;// <-- back to the array of spins new values
					VEC_Y(in,i) = ( ax * (Rz - Az) + ay * (1. + Hy) + az * (Rx + Ax) ) * detMi;
					VEC_Z(in,i) = ( ax * (Ry + Ay) + ay * (Rx - Ax) + az * (1. + Hz) ) * detMi;

					// find max torque:
					detMi = ctx->Heffx[i]*VEC_X(in,i)+ctx->Heffy[i]*VEC_Y(in,i)+ctx->Heffz[i]*VEC_Z(in,i);
					Hx = ctx->Heffx[i]-detMi*VEC_X(in,i);
					Hy = ctx->Heffy[i]-detMi*VEC_Y(in,i);	
					Hz = ctx->Heffz[i]-detMi*VEC_Z(in,i);
					detMi = sqrt(Hx*Hx + Hy*Hy + Hz*Hz);
					if (detMi > ctx->Max_torque[thread]) ctx->Max_torque[thread] = detMi;
					VEC_X(ctx->bS,i)=VEC_X(in,i);
					VEC_Y(ctx->bS,i)=VEC_Y(in,i);
					VEC_Z(ctx->bS,i)=VEC_Z(in,i);
				}
			}
		}
	}	
}


void
StochasticLLG_Heun(	magnoom_ctx *ctx, int thread,
				int naini, 	int nafin,
				int nbini, 	int nbfin,
				int ncini, 	int ncfin)

{
	double *in = ctx->S;
	double *tn = ctx->tS;
	const float *rx = ctx->RNx;
	const float *ry = ctx->RNy;
	const float *rz = ctx->RNz;
	const float alpha = ctx->damping;
	const float h = ctx->t_step;
	const float temperature = ctx->Temperature;
	double rh=sqrt(h);// h = \delta t
	double nx, ny, nz;// components of the unit vector
	double Hx, Hy, Hz;// components of the effective field
	float Cx, Cy, Cz;// spin-torque term
	float Rx, Ry, Rz;// random variables - thermal fluctuations
	double Fx, Fy, Fz;// total matrix
	double detMi;	 // detMi = 1/detM
	double D = sqrt(2.0 * alpha / (1.0 + alpha * alpha) * temperature);
	double Alpha_d = alpha / (1.0 + alpha * alpha * ctx->Precession);
	double Alpha_p = 1.0f / (1.0 + alpha * alpha);

	//electric DC current vector (VCu) and density (Cu)
	Cx = ctx->VCu[0] * ctx->Cu;
	Cy = ctx->VCu[1] * ctx->Cu;
	Cz = ctx->VCu[2] * ctx->Cu;
	GetEffectiveField(ctx, in, naini, nafin, nbini, nbfin, ncini, ncfin);
	//Predictor step:
	int Na = ctx->uABC[0];
	int Nb = ctx->uABC[1];
	int nb1, nc1;
	int i;
	for (int Ip=0; Ip<ctx->AtomsPerBlock; Ip++)
	{
		for (int nc=ncini; nc<ncfin; nc++)//nc(or na or nb) is a neghbor in the direction of c (or a or b)-vector
		{
			nc1 = Na * Nb * nc;
			for (int nb=nbini; nb<nbfin; nb++)
			{
				nb1 = Na * nb;
				for (int na=naini; na<nafin; na++)
				{
					i = Ip + ctx->AtomsPerBlock * ( na + nb1 + nc1 );// index of spin "i"
					nx = VEC_X(in,i);	ny = VEC_Y(in,i);	nz = VEC_Z(in,i);
					Hx = ctx->Heffx[i];	Hy = ctx->Heffy[i];	Hz = ctx->Heffz[i];
					Rx = rx[i];		Ry = ry[i];		Rz = rz[i];

					// deterministic terms of Landau–Lifshitz equation:
					Fx = Alpha_p * Hx + Alpha_d * (ny * Hz - nz * Hy);
					Fy = Alpha_p * Hy + Alpha_d * (nz * Hx - nx * Hz);
					Fz = Alpha_p * Hz + Alpha_d * (nx * Hy - ny * Hx);

					// Spin-torque term
					if (ctx->Cu!=0) {
						Fx = Fx + h * ( Alpha_d * Cx - (ny * Cz - nz * Cy) ); //pay attention to the signe and factors
						Fy = Fy + h * ( Alpha_d * Cy - (nz * Cx - nx * Cz) );
						Fz = Fz + h * ( Alpha_d * Cz - (nx * Cy - ny * Cx) );
					}

					// stochastic terms of Landau–Lifshitz equation:
					if (temperature>0) {
						Fx = Fx + rh * D * ( Alpha_p * Rx + Alpha_d * (ny * Rz - nz * Ry) );
						Fy = Fy + rh * D * ( Alpha_p * Ry + Alpha_d * (nz * Rx - nx * Rz) );
						Fz = Fz + rh * D * ( Alpha_p * Rz + Alpha_d * (nx * Ry - ny * Rx) );
					}
					
					VEC_X(tn,i) = nx - h * (ny * Fz - nz * Fy);
					VEC_Y(tn,i) = ny - h * (nz * Fx - nx * Fz);
					VEC_Z(tn,i) = nz - h * (nx * Fy - ny * Fx);
				}
			}
		}
	}

	magnoom_worker_barrier(ctx, thread);
	GetEffectiveField(ctx, tn, naini, nafin, nbini, nbfin, ncini, ncfin);

	//corrector step:
	for (int Ip=0; Ip<ctx->AtomsPerBlock; Ip++)
	{
		for (int nc=ncini; nc<ncfin; nc++)//nc(a,b)=neghbor in the direction of c(a,b)-vector
		{
			nc1 = Na * Nb * nc;
			for (int nb=nbini; nb<nbfin; nb++)
			{
				nb1 = Na * nb;
				for (int na=naini; na<nafin; na++)
				{
					i = Ip + ctx->AtomsPerBlock * ( na + nb1 + nc1 );	// index of spin "i"
					nx = VEC_X(tn,i);	ny = VEC_Y(tn,i);	nz = VEC_Z(tn,i);	// <-- compare this line to corresponding one in prediction step
					Hx = ctx->Heffx[i];	Hy = ctx->Heffy[i];	Hz = ctx->Heffz[i];	// <-- they are new values for heff
					Rx = rx[i];		Ry = ry[i];		Rz = rz[i];		// <-- they are the same values as in prediction step

					// deterministic terms of Landau–Lifshitz equation:
					Fx = Alpha_p * Hx + Alpha_d * (ny * Hz - nz * Hy);
					Fy = Alpha_p * Hy + Alpha_d * (nz * Hx - nx * Hz);
					Fz = Alpha_p * Hz + Alpha_d * (nx * Hy - ny * Hx);

					// Spin-torque term
					if (ctx->Cu!=0) {
						Fx = Fx + h * ( Alpha_d * Cx - (ny * Cz - nz * Cy) ); //pay attention to the signe and factors
						Fy = Fy + h * ( Alpha_d * Cy - (nz * Cx - nx * Cz) );
						Fz = Fz + h * ( Alpha_d * Cz - (nx * Cy - ny * Cx) );
					}

					// stochastic terms of Landau–Lifshitz equation:
					if (temperature>0) {
						Fx = Fx + rh * D * ( Alpha_p * Rx + Alpha_d * (ny * Rz - nz * Ry) );
						Fy = Fy + rh * D * ( Alpha_p * Ry + Alpha_d * (nz * Rx - nx * Rz) );
						Fz = Fz + rh * D * ( Alpha_p * Rz + Alpha_d * (nx * Ry - ny * Rx) );
					}
					
					nx = 0.5 * ( nx + VEC_X(in,i) - h * (ny * Fz - nz * Fy) );
					ny = 0.5 * ( ny + VEC_Y(in,i) - h * (nz * Fx - nx * Fz) );
					nz = 0.5 * ( nz + VEC_Z(in,i) - h * (nx * Fy - ny * Fx) );

					//normalize spin
					detMi = 1.0 / sqrt(nx*nx + ny*ny + nz*nz);
					VEC_X(in,i) = nx * detMi;
					VEC_Y(in,i) = ny * detMi;
					VEC_Z(in,i) = nz * detMi;

					//find max torque
					detMi = ctx->Heffx[i]*VEC_X(in,i)+ctx->Heffy[i]*VEC_Y(in,i)+ctx->Heffz[i]*VEC_Z(in,i);
					Hx = ctx->Heffx[i]-detMi*VEC_X(in,i);
					Hy = ctx->Heffy[i]-detMi*VEC_Y(in,i);	
					Hz = ctx->Heffz[i]-detMi*VEC_Z(in,i);
					detMi = sqrt(Hx*Hx + Hy*Hy + Hz*Hz);
					if (detMi > ctx->Max_torque[thread]) ctx->Max_torque[thread] = detMi;
					VEC_X(ctx->bS,i)=VEC_X(in,i);
					VEC_Y(ctx->bS,i)=VEC_Y(in,i);
					VEC_Z(ctx->bS,i)=VEC_Z(in,i);
				}
			}
		}
	}	
}



void
StochasticLLG_RK23(	magnoom_ctx *ctx, int thread,
				int naini, 	int nafin,
				int nbini, 	int nbfin,
				int ncini, 	int ncfin)

{
	double *in = ctx->S;
	double *tn = ctx->tS;
	const float *rx = ctx->RNx;
	const float *ry = ctx->RNy;
	const float *rz = ctx->RNz;
	const float alpha = ctx->damping;
	const float h = ctx->t_step;
	const float temperature = ctx->Temperature;
	const float Xi = ctx->Xi;
	const float Curr_u = ctx->Curr_u;
	double rh=sqrt(h);// h = \delta t
	double nx, ny, nz;// components of the unit vector
	double Hx, Hy, Hz;// components of the effective field
	double Hx1, Hy1, Hz1;// components of the effective field
	double Hx2, Hy2, Hz2;// components of the effective field
	float Cx, Cy, Cz;// spin-torque term
	float Rx, Ry, Rz;// random variables - thermal fluctuations
	double Fx, Fy, Fz;// total matrix
	double detMi;	 // detMi = 1/detM
	double D = sqrt(2.0 * alpha / (1.0 + alpha * alpha) * temperature);
	double Alpha_d = alpha / (1.0 + alpha * alpha * ctx->Precession);
	double Alpha_p = 1.0f / (1.0 + alpha * alpha);

	double xi = Xi; 
	double u = Curr_u;
	double c1 = (xi-alpha)*u;
	double c2 = u == 0.0 || alpha == 0.0f ? 0.0 : (1+xi*alpha)*u/alpha;
	// c2 = (alpha>1e-20)? (1+xi*alpha)*u/alpha : 0;

	//electric DC current vector (VCu) and density (Cu)
	Cx = ctx->VCu[0] * ctx->Cu;
	Cy = ctx->VCu[1] * ctx->Cu;
	Cz = ctx->VCu[2] * ctx->Cu;
	GetEffectiveField(ctx, in, naini, nafin, nbini, nbfin, ncini, ncfin);
	//k1:
	int Na = ctx->uABC[0];
	int Nb = ctx->uABC[1];
	int nb1, nc1;
	int i;
	for (int Ip=0; Ip<ctx->AtomsPerBlock; Ip++)
	{
		for (int nc=ncini; nc<ncfin; nc++)//nc(or na or nb) is a neghbor in the direction of c (or a or b)-vector
		{
			nc1 = Na * Nb * nc;
			for (int nb=nbini; nb<nbfin; nb++)
			{
				nb1 = Na * nb;
				for (int na=naini; na<nafin; na++)
				{
					i = Ip + ctx->AtomsPerBlock * ( na + nb1 + nc1 );// index of spin "i"
					nx = VEC_X(in,i);	ny = VEC_Y(in,i);	nz = VEC_Z(in,i);
					Hx = ctx->Heffx[i];	Hy = ctx->Heffy[i];	Hz = ctx->Heffz[i];

					//Li-Zhang first order
					int ip = Ip + ctx->AtomsPerBlock * ( (na +1)%Na + nb1 + nc1 );
					int im = Ip + ctx->AtomsPerBlock * ( (na-1+Na)%Na + nb1 + nc1 );
				
					Hx1 = Hx + c1*(VEC_X(in,ip)-VEC_X(in,im))/2;
					Hy1 = Hy + c1*(VEC_Y(in,ip)-VEC_Y(in,im))/2;
					Hz1 = Hz + c1*(VEC_Z(in,ip)-VEC_Z(in,im))/2;

					Hx2 = Hx + c2*(VEC_X(in,ip)-VEC_X(in,im))/2;
					Hy2 = Hy + c2*(VEC_Y(in,ip)-VEC_Y(in,im))/2;
					Hz2 = Hz + c2*(VEC_Z(in,ip)-VEC_Z(in,im))/2;

					//Li-Zhang second order
					// int ip = Ip + ctx->AtomsPerBlock * ( (na +1)%Na + nb1 + nc1 );
					// int im = Ip + ctx->AtomsPerBlock * ( (na-1+Na)%Na + nb1 + nc1 );
					// int ipp = Ip + ctx->AtomsPerBlock * ( (na +2)%Na + nb1 + nc1 );
					// int imm = Ip + ctx->AtomsPerBlock * ( (na-2+Na)%Na + nb1 + nc1 );

					// Hx1 = Hx + c1*(VEC_X(in,imm)/12 - 2*VEC_X(in,im)/3 + 2*VEC_X(in,ip)/3 -VEC_X(in,ipp)/12);
					// Hy1 = Hy + c1*(VEC_Y(in,imm)/12 - 2*VEC_Y(in,im)/3 + 2*VEC_Y(in,ip)/3 -VEC_Y(in,ipp)/12);
					// Hz1 = Hz + c1*(VEC_Z(in,imm)/12 - 2*VEC_Z(in,im)/3 + 2*VEC_Z(in,ip)/3 -VEC_Z(in,ipp)/12);

					// Hx2 = Hx + c2*(VEC_X(in,imm)/12 - 2*VEC_X(in,im)/3 + 2*VEC_X(in,ip)/3 -VEC_X(in,ipp)/12);
					// Hy2 = Hy + c2*(VEC_Y(in,imm)/12 - 2*VEC_Y(in,im)/3 + 2*VEC_Y(in,ip)/3 -VEC_Y(in,ipp)/12);
					// Hz2 = Hz + c2*(VEC_Z(in,imm)/12 - 2*VEC_Z(in,im)/3 + 2*VEC_Z(in,ip)/3 -VEC_Z(in,ipp)/12);

					Rx = rx[i];		Ry = ry[i];		Rz = rz[i];

					// deterministic terms of Landau–Lifshitz equation:
					// Fx = Alpha_p * Hx + Alpha_d * (ny * Hz - nz * Hy);
					// Fy = Alpha_p * Hy + Alpha_d * (nz * Hx - nx * Hz);
					// Fz = Alpha_p * Hz + Alpha_d * (nx * Hy - ny * Hx);
					Fx = Alpha_p * Hx1 + Alpha_d * (ny * Hz2 - nz * Hy2);
					Fy = Alpha_p * Hy1 + Alpha_d * (nz * Hx2 - nx * Hz2);
					Fz = Alpha_p * Hz1 + Alpha_d * (nx * Hy2 - ny * Hx2);

					// Spin-torque term
					if (ctx->Cu!=0) {
						Fx = Fx + h * ( Alpha_d * Cx - (ny * Cz - nz * Cy) ); //pay attention to the signe and factors
						Fy = Fy + h * ( Alpha_d * Cy - (nz * Cx - nx * Cz) );
						Fz = Fz + h * ( Alpha_d * Cz - (nx * Cy - ny * Cx) );
					}

					// stochastic terms of Landau–Lifshitz equation:
					if (temperature>0) {
						Fx = Fx + rh * D * ( Alpha_p * Rx + Alpha_d * (ny * Rz - nz * Ry) );
						Fy = Fy + rh * D * ( Alpha_p * Ry + Alpha_d * (nz * Rx - nx * Rz) );
						Fz = Fz + rh * D * ( Alpha_p * Rz + Alpha_d * (nx * Ry - ny * Rx) );
					}
					
					VEC_X(tn,i) = nx - 0.5*h * (ny * Fz - nz * Fy);
					VEC_Y(tn,i) = ny - 0.5*h * (nz * Fx - nx * Fz);
					VEC_Z(tn,i) = nz - 0.5*h * (nx * Fy - ny * Fx);
				}
			}
		}
	}

	magnoom_worker_barrier(ctx, thread);
	GetEffectiveField(ctx, tn, naini, nafin, nbini, nbfin, ncini, ncfin);

	//corrector step:
	for (int Ip=0; Ip<ctx->AtomsPerBlock; Ip++)
	{
		for (int nc=ncini; nc<ncfin; nc++)//nc(a,b)=neghbor in the direction of c(a,b)-vector
		{
			nc1 = Na * Nb * nc;
			for (int nb=nbini; nb<nbfin; nb++)
			{
				nb1 = Na * nb;
				for (int na=naini; na<nafin; na++)
				{
					i = Ip + ctx->AtomsPerBlock * ( na + nb1 + nc1 );	// index of spin "i"
					int ip = Ip + ctx->AtomsPerBlock * ( (na +1)%Na + nb1 + nc1 );
					int im = Ip + ctx->AtomsPerBlock * ( (na-1+Na)%Na + nb1 + nc1 );
					// int ipp = Ip + ctx->AtomsPerBlock * ( (na +2)%Na + nb1 + nc1 );
					// int imm = Ip + ctx->AtomsPerBlock * ( (na-2+Na)%Na + nb1 + nc1 );
					nx = VEC_X(tn,i);	ny = VEC_Y(tn,i);	nz = VEC_Z(tn,i);	// <-- compare this line to corresponding one in prediction step
					Hx = ctx->Heffx[i];	Hy = ctx->Heffy[i];	Hz = ctx->Heffz[i];	// <-- they are new values for heff
					Hx1 = Hx + c1*(VEC_X(tn,ip)-VEC_X(tn,im))/2;
					Hy1 = Hy + c1*(VEC_Y(tn,ip)-VEC_Y(tn,im))/2;
					Hz1 = Hz + c1*(VEC_Z(tn,ip)-VEC_Z(tn,im))/2;

					Hx2 = Hx + c2*(VEC_X(tn,ip)-VEC_X(tn,im))/2;
					Hy2 = Hy + c2*(VEC_Y(tn,ip)-VEC_Y(tn,im))/2;
					Hz2 = Hz + c2*(VEC_Z(tn,ip)-VEC_Z(tn,im))/2;
					Rx = rx[i];		Ry = ry[i];		Rz = rz[i];		// <-- they are the same values as in prediction step

					// deterministic terms of Landau–Lifshitz equation:
					// Fx = Alpha_p * Hx + Alpha_d * (ny * Hz - nz * Hy);
					// Fy = Alpha_p * Hy + Alpha_d * (nz * Hx - nx * Hz);
					// Fz = Alpha_p * Hz + Alpha_d * (nx * Hy - ny * Hx);
					Fx = Alpha_p * Hx1 + Alpha_d * (ny * Hz2 - nz * Hy2);
					Fy = Alpha_p * Hy1 + Alpha_d * (nz * Hx2 - nx * Hz2);
					Fz = Alpha_p * Hz1 + Alpha_d * (nx * Hy2 - ny * Hx2);

					// Spin-torque term
					if (ctx->Cu!=0) {
						Fx = Fx + h * ( Alpha_d * Cx - (ny * Cz - nz * Cy) ); //pay attention to the signe and factors
						Fy = Fy + h * ( Alpha_d * Cy - (nz * Cx - nx * Cz) );
						Fz = Fz + h * ( Alpha_d * Cz - (nx * Cy - ny * Cx) );
					}

					// stochastic terms of Landau–Lifshitz equation:
					if (temperature>0) {
						Fx = Fx + rh * D * ( Alpha_p * Rx + Alpha_d * (ny * Rz - nz * Ry) );
						Fy = Fy + rh * D * ( Alpha_p * Ry + Alpha_d * (nz * Rx - nx * Rz) );
						Fz = Fz + rh * D * ( Alpha_p * Rz + Alpha_d * (nx * Ry - ny * Rx) );
					}
					//(Cx,Cy,Cz here are used as temp variables)
					Cx = VEC_X(in,i) - h * (ny * Fz - nz * Fy);
					Cy = VEC_Y(in,i) - h * (nz * Fx - nx * Fz);
					Cz = VEC_Z(in,i) - h * (nx * Fy - ny * Fx);

					//normalize spin
					detMi = 1.0 / sqrt(Cx*Cx + Cy*Cy + Cz*Cz);
					VEC_X(in,i) = Cx * detMi;
					VEC_Y(in,i) = Cy * detMi;
					VEC_Z(in,i) = Cz * detMi;

					//find max torque
					detMi = ctx->Heffx[i]*VEC_X(in,i)+ctx->Heffy[i]*VEC_Y(in,i)+ctx->Heffz[i]*VEC_Z(in,i);
					Hx = ctx->Heffx[i]-detMi*VEC_X(in,i);
					Hy = ctx->Heffy[i]-detMi*VEC_Y(in,i);	
					Hz = ctx->Heffz[i]-detMi*VEC_Z(in,i);
					detMi = sqrt(Hx*Hx + Hy*Hy + Hz*Hz);
					if (detMi > ctx->Max_torque[thread]) ctx->Max_torque[thread] = detMi;
					VEC_X(ctx->bS,i)=VEC_X(in,i);
					VEC_Y(ctx->bS,i)=VEC_Y(in,i);
					VEC_Z(ctx->bS,i)=VEC_Z(in,i);
				}
			}
		}
	}	
}


void
StochasticLLG_RK4(	magnoom_ctx *ctx, int thread,
				int naini, 	int nafin,
				int nbini, 	int nbfin,
				int ncini, 	int ncfin)
{
	double *in = ctx->S;
	double *tn = ctx->tS;
	double *next = ctx->rkS;
	const float *rx = ctx->RNx;
	const float *ry = ctx->RNy;
	const float *rz = ctx->RNz;
	const float alpha = ctx->damping;
	const float h = ctx->t_step;
	const float temperature = ctx->Temperature;
	const float Xi = ctx->Xi;
	const float Curr_u = ctx->Curr_u;
	double rh=sqrt(h);// h = \delta t
	double nx, ny, nz;// components of the unit vector
	double Hx, Hy, Hz;// components of the effective field
	double Hx1, Hy1, Hz1;// components of the effective field
	double Hx2, Hy2, Hz2;// components of the effective field
	float Cx, Cy, Cz;// spin-torque term
	float Rx, Ry, Rz;// random variables - thermal fluctuations
	double Fx, Fy, Fz;// total matrix
	double detMi;	 // detMi = 1/detM
	double D = sqrt(2.0 * alpha / (1.0 + alpha * alpha) * temperature);
	double Alpha_d = alpha / (1.0 + alpha * alpha * ctx->Precession);
	double Alpha_p = 1.0f / (1.0 + alpha * alpha);

	//electric DC current vector (VCu) and density (Cu)
	Cx = ctx->VCu[0] * ctx->Cu;
	Cy = ctx->VCu[1] * ctx->Cu;
	Cz = ctx->VCu[2] * ctx->Cu;
	GetEffectiveField(ctx, in, naini, nafin, nbini, nbfin, ncini, ncfin);
	//Predictor step:
	int Na = ctx->uABC[0];
	int Nb = ctx->uABC[1];
	int nb1, nc1;
	int i;

	double xi = Xi; 
	double u = Curr_u;
	double c1 = (xi-alpha)*u;
	double c2 = u == 0.0 || alpha == 0.0f ? 0.0 : (1+xi*alpha)*u/alpha;

	//k1:
	for (int Ip=0; Ip<ctx->AtomsPerBlock; Ip++)
	{
		for (int nc=ncini; nc<ncfin; nc++)//nc(or na or nb) is a neghbor in the direction of c (or a or b)-vector
		{
			nc1 = Na * Nb * nc;
			for (int nb=nbini; nb<nbfin; nb++)
			{
				nb1 = Na * nb;
				for (int na=naini; na<nafin; na++)
				{
					i = Ip + ctx->AtomsPerBlock * ( na + nb1 + nc1 );// index of spin "i"
					nx = VEC_X(in,i);	ny = VEC_Y(in,i);	nz = VEC_Z(in,i);
					Hx = ctx->Heffx[i];	Hy = ctx->Heffy[i];	Hz = ctx->Heffz[i];
					Rx = rx[i];		Ry = ry[i];		Rz = rz[i];

					int ip = Ip + ctx->AtomsPerBlock * ( (na +1)%Na + nb1 + nc1 );
					int im = Ip + ctx->AtomsPerBlock * ( (na-1+Na)%Na + nb1 + nc1 );
					// int ipp = Ip + ctx->AtomsPerBlock * ( (na +2)%Na + nb1 + nc1 );
					// int imm = Ip + ctx->AtomsPerBlock * ( (na-2+Na)%Na + nb1 + nc1 );
					// Hx1 = Hx + c1*(VEC_X(in,imm)/12 - 2*VEC_X(in,im)/3 + 2*VEC_X(in,ip)/3 -VEC_X(in,ipp)/12);
					// Hy1 = Hy + c1*(VEC_Y(in,imm)/12 - 2*VEC_Y(in,im)/3 + 2*VEC_Y(in,ip)/3 -VEC_Y(in,ipp)/12);
					// Hz1 = Hz + c1*(VEC_Z(in,imm)/12 - 2*VEC_Z(in,im)/3 + 2*VEC_Z(in,ip)/3 -VEC_Z(in,ipp)/12);

					// Hx2 = Hx + c2*(VEC_X(in,imm)/12 - 2*VEC_X(in,im)/3 + 2*VEC_X(in,ip)/3 -VEC_X(in,ipp)/12);
					// Hy2 = Hy + c2*(VEC_Y(in,imm)/12 - 2*VEC_Y(in,im)/3 + 2*VEC_Y(in,ip)/3 -VEC_Y(in,ipp)/12);
					// Hz2 = Hz + c2*(VEC_Z(in,imm)/12 - 2*VEC_Z(in,im)/3 + 2*VEC_Z(in,ip)/3 -VEC_Z(in,ipp)/12);
					Hx1 = Hx + c1*(VEC_X(in,ip)-VEC_X(in,im))/2;
					Hy1 = Hy + c1*(VEC_Y(in,ip)-VEC_Y(in,im))/2;
					Hz1 = Hz + c1*(VEC_Z(in,ip)-VEC_Z(in,im))/2;

					Hx2 = Hx + c2*(VEC_X(in,ip)-VEC_X(in,im))/2;
					Hy2 = Hy + c2*(VEC_Y(in,ip)-VEC_Y(in,im))/2;
					Hz2 = Hz + c2*(VEC_Z(in,ip)-VEC_Z(in,im))/2;

					// deterministic terms of Landau–Lifshitz equation:
					// Fx = Alpha_p * Hx + Alpha_d * (ny * Hz - nz * Hy);
					// Fy = Alpha_p * Hy + Alpha_d * (nz * Hx - nx * Hz);
					// Fz = Alpha_p * Hz + Alpha_d * (nx * Hy - ny * Hx);
					Fx = Alpha_p * Hx1 + Alpha_d * (ny * Hz2 - nz * Hy2);
					Fy = Alpha_p * Hy1 + Alpha_d * (nz * Hx2 - nx * Hz2);
					Fz = Alpha_p * Hz1 + Alpha_d * (nx * Hy2 - ny * Hx2);

					// Spin-torque term
					if (ctx->Cu!=0) {
						Fx = Fx + h * ( Alpha_d * Cx - (ny * Cz - nz * Cy) ); //pay attention to the signe and factors
						Fy = Fy + h * ( Alpha_d * Cy - (nz * Cx - nx * Cz) );
						Fz = Fz + h * ( Alpha_d * Cz - (nx * Cy - ny * Cx) );
					}

					// stochastic terms of Landau–Lifshitz equation:
					if (temperature>0) {
						Fx = Fx + rh * D * ( Alpha_p * Rx + Alpha_d * (ny * Rz - nz * Ry) );
						Fy = Fy + rh * D * ( Alpha_p * Ry + Alpha_d * (nz * Rx - nx * Rz) );
						Fz = Fz + rh * D * ( Alpha_p * Rz + Alpha_d * (nx * Ry - ny * Rx) );
					}
					//k1 (Cx,Cy,Cz here are used as temp variables)
					Cx = - h * (ny * Fz - nz * Fy);
					Cy = - h * (nz * Fx - nx * Fz);
					Cz = - h * (nx * Fy - ny * Fx);
					//save k1/6 in global temp array
					VEC_X(ctx->t2S,i) = Cx/6.0;
					VEC_Y(ctx->t2S,i) = Cy/6.0;
					VEC_Z(ctx->t2S,i) = Cz/6.0;
					//y_n+k1/2 will be used on the next step
					VEC_X(tn,i) = nx + Cx*0.5;
					VEC_Y(tn,i) = ny + Cy*0.5;
					VEC_Z(tn,i) = nz + Cz*0.5;				
				}
			}
		}
	}
	//Heff(y_n+k1/2):
	magnoom_worker_barrier(ctx, thread);
	GetEffectiveField(ctx, tn, naini, nafin, nbini, nbfin, ncini, ncfin);

	//k2:
	for (int Ip=0; Ip<ctx->AtomsPerBlock; Ip++)
	{
		for (int nc=ncini; nc<ncfin; nc++)//nc(or na or nb) is a neghbor in the direction of c (or a or b)-vector
		{
			nc1 = Na * Nb * nc;
			for (int nb=nbini; nb<nbfin; nb++)
			{
				nb1 = Na * nb;
				for (int na=naini; na<nafin; na++)
				{
					i = Ip + ctx->AtomsPerBlock * ( na + nb1 + nc1 );// index of spin "i"
					nx = VEC_X(tn,i);	ny = VEC_Y(tn,i);	nz = VEC_Z(tn,i);
					Hx = ctx->Heffx[i];	Hy = ctx->Heffy[i];	Hz = ctx->Heffz[i];
					Rx = rx[i];		Ry = ry[i];		Rz = rz[i];

					int ip = Ip + ctx->AtomsPerBlock * ( (na +1)%Na + nb1 + nc1 );
					int im = Ip + ctx->AtomsPerBlock * ( (na-1+Na)%Na + nb1 + nc1 );
					// int ipp = Ip + ctx->AtomsPerBlock * ( (na +2)%Na + nb1 + nc1 );
					// int imm = Ip + ctx->AtomsPerBlock * ( (na-2+Na)%Na + nb1 + nc1 );
					// Hx1 = Hx + c1*(VEC_X(tn,imm)/12 - 2*VEC_X(tn,im)/3 + 2*VEC_X(tn,ip)/3 -VEC_X(tn,ipp)/12);
					// Hy1 = Hy + c1*(VEC_Y(tn,imm)/12 - 2*VEC_Y(tn,im)/3 + 2*VEC_Y(tn,ip)/3 -VEC_Y(tn,ipp)/12);
					// Hz1 = Hz + c1*(VEC_Z(tn,imm)/12 - 2*VEC_Z(tn,im)/3 + 2*VEC_Z(tn,ip)/3 -VEC_Z(tn,ipp)/12);

					// Hx2 = Hx + c2*(VEC_X(tn,imm)/12 - 2*VEC_X(tn,im)/3 + 2*VEC_X(tn,ip)/3 -VEC_X(tn,ipp)/12);
					// Hy2 = Hy + c2*(VEC_Y(tn,imm)/12 - 2*VEC_Y(tn,im)/3 + 2*VEC_Y(tn,ip)/3 -VEC_Y(tn,ipp)/12);
					// Hz2 = Hz + c2*(VEC_Z(tn,imm)/12 - 2*VEC_Z(tn,im)/3 + 2*VEC_Z(tn,ip)/3 -VEC_Z(tn,ipp)/12);

					Hx1 = Hx + c1*(VEC_X(tn,ip)-VEC_X(tn,im))/2;
					Hy1 = Hy + c1*(VEC_Y(tn,ip)-VEC_Y(tn,im))/2;
					Hz1 = Hz + c1*(VEC_Z(tn,ip)-VEC_Z(tn,im))/2;

					Hx2 = Hx + c2*(VEC_X(tn,ip)-VEC_X(tn,im))/2;
					Hy2 = Hy + c2*(VEC_Y(tn,ip)-VEC_Y(tn,im))/2;
					Hz2 = Hz + c2*(VEC_Z(tn,ip)-VEC_Z(tn,im))/2;

					// deterministic terms of Landau–Lifshitz equation:
					// Fx = Alpha_p * Hx + Alpha_d * (ny * Hz - nz * Hy);
					// Fy = Alpha_p * Hy + Alpha_d * (nz * Hx - nx * Hz);
					// Fz = Alpha_p * Hz + Alpha_d * (nx * Hy - ny * Hx);
					Fx = Alpha_p * Hx1 + Alpha_d * (ny * Hz2 - nz * Hy2);
					Fy = Alpha_p * Hy1 + Alpha_d * (nz * Hx2 - nx * Hz2);
					Fz = Alpha_p * Hz1 + Alpha_d * (nx * Hy2 - ny * Hx2);


					// Spin-torque term
					if (ctx->Cu!=0) {
						Fx = Fx + h * ( Alpha_d * Cx - (ny * Cz - nz * Cy) ); //pay attention to the signe and factors
						Fy = Fy + h * ( Alpha_d * Cy - (nz * Cx - nx * Cz) );
						Fz = Fz + h * ( Alpha_d * Cz - (nx * Cy - ny * Cx) );
					}

					// stochastic terms of Landau–Lifshitz equation:
					if (temperature>0) {
						Fx = Fx + rh * D * ( Alpha_p * Rx + Alpha_d * (ny * Rz - nz * Ry) );
						Fy = Fy + rh * D * ( Alpha_p * Ry + Alpha_d * (nz * Rx - nx * Rz) );
						Fz = Fz + rh * D * ( Alpha_p * Rz + Alpha_d * (nx * Ry - ny * Rx) );
					}
					//k2 (Cx,Cy,Cz here are used as temp variables)
					Cx = - h * (ny * Fz - nz * Fy);
					Cy = - h * (nz * Fx - nx * Fz);
					Cz = - h * (nx * Fy - ny * Fx);
					//save k2/3 in global temp array
					VEC_X(ctx->t2S,i)+= Cx/3.0;
					VEC_Y(ctx->t2S,i)+= Cy/3.0;
					VEC_Z(ctx->t2S,i)+= Cz/3.0;
					//y_n+k2/2 will be used on the next step
					VEC_X(next,i) = VEC_X(in,i) + Cx*0.5;
					VEC_Y(next,i) = VEC_Y(in,i) + Cy*0.5;
					VEC_Z(next,i) = VEC_Z(in,i) + Cz*0.5;
				}
			}
		}
	}
	//Heff(y_n+k2/2):
	magnoom_worker_barrier(ctx, thread);
	GetEffectiveField(ctx, next, naini, nafin, nbini, nbfin, ncini, ncfin);
	//k3:
	for (int Ip=0; Ip<ctx->AtomsPerBlock; Ip++)
	{
		for (int nc=ncini; nc<ncfin; nc++)//nc(or na or nb) is a neghbor in the direction of c (or a or b)-vector
		{
			nc1 = Na * Nb * nc;
			for (int nb=nbini; nb<nbfin; nb++)
			{
				nb1 = Na * nb;
				for (int na=naini; na<nafin; na++)
				{
					i = Ip + ctx->AtomsPerBlock * ( na + nb1 + nc1 );// index of spin "i"
					nx = VEC_X(next,i);	ny = VEC_Y(next,i);	nz = VEC_Z(next,i);
					Hx = ctx->Heffx[i];	Hy = ctx->Heffy[i];	Hz = ctx->Heffz[i];
					Rx = rx[i];		Ry = ry[i];		Rz = rz[i];

					int ip = Ip + ctx->AtomsPerBlock * ( (na +1)%Na + nb1 + nc1 );
					int im = Ip + ctx->AtomsPerBlock * ( (na-1+Na)%Na + nb1 + nc1 );
					// int ipp = Ip + ctx->AtomsPerBlock * ( (na +2)%Na + nb1 + nc1 );
					// int imm = Ip + ctx->AtomsPerBlock * ( (na-2+Na)%Na + nb1 + nc1 );
					// Hx1 = Hx + c1*(VEC_X(tn,imm)/12 - 2*VEC_X(tn,im)/3 + 2*VEC_X(tn,ip)/3 -VEC_X(tn,ipp)/12);
					// Hy1 = Hy + c1*(VEC_Y(tn,imm)/12 - 2*VEC_Y(tn,im)/3 + 2*VEC_Y(tn,ip)/3 -VEC_Y(tn,ipp)/12);
					// Hz1 = Hz + c1*(VEC_Z(tn,imm)/12 - 2*VEC_Z(tn,im)/3 + 2*VEC_Z(tn,ip)/3 -VEC_Z(tn,ipp)/12);

					// Hx2 = Hx + c2*(VEC_X(tn,imm)/12 - 2*VEC_X(tn,im)/3 + 2*VEC_X(tn,ip)/3 -VEC_X(tn,ipp)/12);
					// Hy2 = Hy + c2*(VEC_Y(tn,imm)/12 - 2*VEC_Y(tn,im)/3 + 2*VEC_Y(tn,ip)/3 -VEC_Y(tn,ipp)/12);
					// Hz2 = Hz + c2*(VEC_Z(tn,imm)/12 - 2*VEC_Z(tn,im)/3 + 2*VEC_Z(tn,ip)/3 -VEC_Z(tn,ipp)/12);
					Hx1 = Hx + c1*(VEC_X(next,ip)-VEC_X(next,im))/2;
					Hy1 = Hy + c1*(VEC_Y(next,ip)-VEC_Y(next,im))/2;
					Hz1 = Hz + c1*(VEC_Z(next,ip)-VEC_Z(next,im))/2;

					Hx2 = Hx + c2*(VEC_X(next,ip)-VEC_X(next,im))/2;
					Hy2 = Hy + c2*(VEC_Y(next,ip)-VEC_Y(next,im))/2;
					Hz2 = Hz + c2*(VEC_Z(next,ip)-VEC_Z(next,im))/2;

					// deterministic terms of Landau–Lifshitz equation:
					// Fx = Alpha_p * Hx + Alpha_d * (ny * Hz - nz * Hy);
					// Fy = Alpha_p * Hy + Alpha_d * (nz * Hx - nx * Hz);
					// Fz = Alpha_p * Hz + Alpha_d * (nx * Hy - ny * Hx);
					Fx = Alpha_p * Hx1 + Alpha_d * (ny * Hz2 - nz * Hy2);
					Fy = Alpha_p * Hy1 + Alpha_d * (nz * Hx2 - nx * Hz2);
					Fz = Alpha_p * Hz1 + Alpha_d * (nx * Hy2 - ny * Hx2);

					// Spin-torque term
					if (ctx->Cu!=0) {
						Fx = Fx + h * ( Alpha_d * Cx - (ny * Cz - nz * Cy) ); //pay attention to the signe and factors
						Fy = Fy + h * ( Alpha_d * Cy - (nz * Cx - nx * Cz) );
						Fz = Fz + h * ( Alpha_d * Cz - (nx * Cy - ny * Cx) );
					}

					// stochastic terms of Landau–Lifshitz equation:
					if (temperature>0) {
						Fx = Fx + rh * D * ( Alpha_p * Rx + Alpha_d * (ny * Rz - nz * Ry) );
						Fy = Fy + rh * D * ( Alpha_p * Ry + Alpha_d * (nz * Rx - nx * Rz) );
						Fz = Fz + rh * D * ( Alpha_p * Rz + Alpha_d * (nx * Ry - ny * Rx) );
					}
					//k3 (Cx,Cy,Cz here are used as temp variables)
					Cx = - h * (ny * Fz - nz * Fy);
					Cy = - h * (nz * Fx - nx * Fz);
					Cz = - h * (nx * Fy - ny * Fx);
					//save k2/3 in global temp array
					VEC_X(ctx->t2S,i)+= Cx/3.0;
					VEC_Y(ctx->t2S,i)+= Cy/3.0;
					VEC_Z(ctx->t2S,i)+= Cz/3.0;
					//y_n+k3 will be used on the next step
					VEC_X(tn,i) = VEC_X(in,i) + Cx;
					VEC_Y(tn,i) = VEC_Y(in,i) + Cy;
					VEC_Z(tn,i) = VEC_Z(in,i) + Cz;							
				}
			}
		}
	}
	//Heff(y_n+k3):
	magnoom_worker_barrier(ctx, thread);
	GetEffectiveField(ctx, tn, naini, nafin, nbini, nbfin, ncini, ncfin);
	//k4:
	for (int Ip=0; Ip<ctx->AtomsPerBlock; Ip++)
	{
		for (int nc=ncini; nc<ncfin; nc++)//nc(or na or nb) is a neghbor in the direction of c (or a or b)-vector
		{
			nc1 = Na * Nb * nc;
			for (int nb=nbini; nb<nbfin; nb++)
			{
				nb1 = Na * nb;
				for (int na=naini; na<nafin; na++)
				{
					i = Ip + ctx->AtomsPerBlock * ( na + nb1 + nc1 );// index of spin "i"
					nx = VEC_X(tn,i);	ny = VEC_Y(tn,i);	nz = VEC_Z(tn,i);
					Hx = ctx->Heffx[i];	Hy = ctx->Heffy[i];	Hz = ctx->Heffz[i];
					Rx = rx[i];		Ry = ry[i];		Rz = rz[i];

					int ip = Ip + ctx->AtomsPerBlock * ( (na +1)%Na + nb1 + nc1 );
					int im = Ip + ctx->AtomsPerBlock * ( (na-1+Na)%Na + nb1 + nc1 );
					// int ipp = Ip + ctx->AtomsPerBlock * ( (na +2)%Na + nb1 + nc1 );
					// int imm = Ip + ctx->AtomsPerBlock * ( (na-2+Na)%Na + nb1 + nc1 );
					// Hx1 = Hx + c1*(VEC_X(tn,imm)/12 - 2*VEC_X(tn,im)/3 + 2*VEC_X(tn,ip)/3 -VEC_X(tn,ipp)/12);
					// Hy1 = Hy + c1*(VEC_Y(tn,imm)/12 - 2*VEC_Y(tn,im)/3 + 2*VEC_Y(tn,ip)/3 -VEC_Y(tn,ipp)/12);
					// Hz1 = Hz + c1*(VEC_Z(tn,imm)/12 - 2*VEC_Z(tn,im)/3 + 2*VEC_Z(tn,ip)/3 -VEC_Z(tn,ipp)/12);

					// Hx2 = Hx + c2*(VEC_X(tn,imm)/12 - 2*VEC_X(tn,im)/3 + 2*VEC_X(tn,ip)/3 -VEC_X(tn,ipp)/12);
					// Hy2 = Hy + c2*(VEC_Y(tn,imm)/12 - 2*VEC_Y(tn,im)/3 + 2*VEC_Y(tn,ip)/3 -VEC_Y(tn,ipp)/12);
					// Hz2 = Hz + c2*(VEC_Z(tn,imm)/12 - 2*VEC_Z(tn,im)/3 + 2*VEC_Z(tn,ip)/3 -VEC_Z(tn,ipp)/12);
					Hx1 = Hx + c1*(VEC_X(tn,ip)-VEC_X(tn,im))/2;
					Hy1 = Hy + c1*(VEC_Y(tn,ip)-VEC_Y(tn,im))/2;
					Hz1 = Hz + c1*(VEC_Z(tn,ip)-VEC_Z(tn,im))/2;

					Hx2 = Hx + c2*(VEC_X(tn,ip)-VEC_X(tn,im))/2;
					Hy2 = Hy + c2*(VEC_Y(tn,ip)-VEC_Y(tn,im))/2;
					Hz2 = Hz + c2*(VEC_Z(tn,ip)-VEC_Z(tn,im))/2;
					// deterministic terms of Landau–Lifshitz equation:
					// Fx = Alpha_p * Hx + Alpha_d * (ny * Hz - nz * Hy);
					// Fy = Alpha_p * Hy + Alpha_d * (nz * Hx - nx * Hz);
					// Fz = Alpha_p * Hz + Alpha_d * (nx * Hy - ny * Hx);
					Fx = Alpha_p * Hx1 + Alpha_d * (ny * Hz2 - nz * Hy2);
					Fy = Alpha_p * Hy1 + Alpha_d * (nz * Hx2 - nx * Hz2);
					Fz = Alpha_p * Hz1 + Alpha_d * (nx * Hy2 - ny * Hx2);

					// Spin-torque term
					if (ctx->Cu!=0) {
						Fx = Fx + h * ( Alpha_d * Cx - (ny * Cz - nz * Cy) ); //pay attention to the signe and factors
						Fy = Fy + h * ( Alpha_d * Cy - (nz * Cx - nx * Cz) );
						Fz = Fz + h * ( Alpha_d * Cz - (nx * Cy - ny * Cx) );
					}

					// stochastic terms of Landau–Lifshitz equation:
					if (temperature>0) {
						Fx = Fx + rh * D * ( Alpha_p * Rx + Alpha_d * (ny * Rz - nz * Ry) );
						Fy = Fy + rh * D * ( Alpha_p * Ry + Alpha_d * (nz * Rx - nx * Rz) );
						Fz = Fz + rh * D * ( Alpha_p * Rz + Alpha_d * (nx * Ry - ny * Rx) );
					}
					//k4 (Cx,Cy,Cz here are used as temp variables)
					Cx = - h * (ny * Fz - nz * Fy);
					Cy = - h * (nz * Fx - nx * Fz);
					Cz = - h * (nx * Fy - ny * Fx);
					//save k4/6 in global temp array
					VEC_X(ctx->t2S,i)+= Cx/6.0;
					VEC_Y(ctx->t2S,i)+= Cy/6.0;
					VEC_Z(ctx->t2S,i)+= Cz/6.0;
					//y_{n+1}=y_n+k1/6+k2/3+k3/3+k4/6 - final step:
					VEC_X(in,i)+= VEC_X(ctx->t2S,i);
					VEC_Y(in,i)+= VEC_Y(ctx->t2S,i);
					VEC_Z(in,i)+= VEC_Z(ctx->t2S,i);
					nx = VEC_X(in,i);	ny = VEC_Y(in,i);	nz = VEC_Z(in,i);
					//normalize spin
					detMi = 1.0 / sqrt(nx*nx + ny*ny + nz*nz);
					VEC_X(in,i) = nx * detMi;
					VEC_Y(in,i) = ny * detMi;
					VEC_Z(in,i) = nz * detMi;

					//find max torque
					detMi = ctx->Heffx[i]*VEC_X(in,i)+ctx->Heffy[i]*VEC_Y(in,i)+ctx->Heffz[i]*VEC_Z(in,i);
					Hx = ctx->Heffx[i]-detMi*VEC_X(in,i);
					Hy = ctx->Heffy[i]-detMi*VEC_Y(in,i);	
					Hz = ctx->Heffz[i]-detMi*VEC_Z(in,i);
					detMi = sqrt(Hx*Hx + Hy*Hy + Hz*Hz);
					if (detMi > ctx->Max_torque[thread]) ctx->Max_torque[thread] = detMi;	
					VEC_X(ctx->bS,i)=VEC_X(in,i);
					VEC_Y(ctx->bS,i)=VEC_Y(in,i);
					VEC_Z(ctx->bS,i)=VEC_Z(in,i);					
				}
			}
		}
	}
}

void Relax(	magnoom_ctx *ctx, int thread,
				int naini, 	int nafin,
				int nbini, 	int nbfin,
				int ncini, 	int ncfin)
{
	double *in = ctx->S;
	bool *proj = ctx->Proj;
	
	GetEffectiveField(ctx, in, naini, nafin, nbini, nbfin, ncini, ncfin);

	double Hx, Hy, Hz, temp;// components of the effective field
	int Na = ctx->uABC[0];
	int Nb = ctx->uABC[1];
	int nb1, nc1;
	int i,s;

	double ALPHA,g1,g2, d1,d2;
	ALPHA = ctx->damping;

	for (int Ip=0; Ip<ctx->AtomsPerBlock; Ip++)
	{
		for (int nc=ncini; nc<ncfin; nc++)//nc(or na or nb) is a neghbor in the direction of c (or a or b)-vector
		{
			nc1 = Na * Nb * nc;
			for (int nb=nbini; nb<nbfin; nb++)
			{
				nb1 = Na * nb;
				for (int na=naini; na<nafin; na++)
				{
					i = Ip + ctx->AtomsPerBlock * ( na + nb1 + nc1 );// index of spin "i"
					Hx = ctx->Heffx[i];	Hy = ctx->Heffy[i];	Hz = ctx->Heffz[i];
					//find max torque
					temp = ctx->Heffx[i]*VEC_X(in,i)+ctx->Heffy[i]*VEC_Y(in,i)+ctx->Heffz[i]*VEC_Z(in,i);
					Hx = ctx->Heffx[i]-temp*VEC_X(in,i);
					Hy = ctx->Heffy[i]-temp*VEC_Y(in,i);
					Hz = ctx->Heffz[i]-temp*VEC_Z(in,i);
					temp = sqrt(Hx*Hx + Hy*Hy + Hz*Hz);
					if (temp > ctx->Max_torque[thread]) ctx->Max_torque[thread] = temp;
					// constant step descent
					proj[i] = (VEC_Z(in,i)>0)? true : false;
					s = (proj[i])? 1 : -1;
					g1 = VEC_X(in,i)/(1+s*VEC_Z(in,i));
					g2 = VEC_Y(in,i)/(1+s*VEC_Z(in,i));

					Hx = ctx->Heffx[i];	Hy = ctx->Heffy[i];	Hz = ctx->Heffz[i];
					d1 = (Hx*(VEC_Y(in,i)*VEC_Y(in,i)+s*VEC_Z(in,i)*(1+s*VEC_Z(in,i))) - Hy*VEC_X(in,i)*VEC_Y(in,i) - Hz*VEC_X(in,i)*(s+VEC_Z(in,i)));
					d2 = (-Hx*VEC_X(in,i)*VEC_Y(in,i) + Hy*(VEC_X(in,i)*VEC_X(in,i)+s*VEC_Z(in,i)*(1+s*VEC_Z(in,i))) - Hz*VEC_Y(in,i)*(s+VEC_Z(in,i)));

					g1 += ALPHA*d1;
					g2 += ALPHA*d2;

					VEC_X(in,i) = 2*g1/(1+g1*g1+g2*g2);
					VEC_Y(in,i) = 2*g2/(1+g1*g1+g2*g2);
					VEC_Z(in,i) = s*(1-g1*g1-g2*g2)/(1+g1*g1+g2*g2);

					VEC_X(ctx->bS,i)=VEC_X(in,i);
					VEC_Y(ctx->bS,i)=VEC_Y(in,i);
					VEC_Z(ctx->bS,i)=VEC_Z(in,i);
				}
			}
		}
	}
}



double EvaluateBextAC(magnoom_ctx *ctx)
{
	double R=0;
	double temp;
	switch (ctx->BextACWaveform){
		case BEXT_AC_SIN:
			R = ctx->BextACEnabled*ctx->BextACAmplitude*sin(ctx->BextACOmega*ctx->t_step*ctx->ITERATION);
		break;

		case BEXT_AC_GAUSSIAN:
			temp = (ctx->t_step*ctx->ITERATION-ctx->BextACTimeOffset)/ctx->BextACPulseWidth;
			R = ctx->BextACEnabled*ctx->BextACAmplitude*exp(-0.5*temp*temp );
		break;

		case BEXT_AC_SINC:
			if (ctx->t_step*ctx->ITERATION<=ctx->BextACPulseWidth){
				if (ctx->t_step*ctx->ITERATION==ctx->BextACTimeOffset){
					R = ctx->BextACEnabled*ctx->BextACAmplitude*1;
				}else{
					R = ctx->BextACEnabled*ctx->BextACAmplitude*sin((ctx->BextACOmega*(ctx->t_step*ctx->ITERATION-ctx->BextACTimeOffset)))/(ctx->BextACOmega*(ctx->t_step*ctx->ITERATION-ctx->BextACTimeOffset));
				}
			}else{ R = 0.0; }
		break;

		case BEXT_AC_CIRCULAR:
			if(ctx->BextACEnabled!=0){
				ctx->BextACDirection[0] = cos(ctx->BextACOmega*ctx->t_step*ctx->ITERATION);
				ctx->BextACDirection[1] = sin(ctx->BextACOmega*ctx->t_step*ctx->ITERATION);
				ctx->BextACDirection[2] = 0.0f;
				R = ctx->BextACEnabled*ctx->BextACAmplitude;
			}
		break;
	}
	ctx->BextAC[0] = R * ctx->BextACDirection[0];
	ctx->BextAC[1] = R * ctx->BextACDirection[1];
	ctx->BextAC[2] = R * ctx->BextACDirection[2];
	return R;
}

static void RecordBextACMode(magnoom_ctx *ctx)
{
	// Save mode snapshot.
	if (ctx->BextACModeRecording*ctx->BextACEnabled!=0){
		float phase = ctx->BextACOmega*ctx->t_step*ctx->ITERATION*iTPI;
		phase = phase - floor(phase);
		int Im=-1;
		float tolerance=0.5*ctx->t_step/ctx->BextACPeriod;

		for (int i=0; i<=ctx->num_images; i++){
			if ( ABS(phase - i/(float)ctx->num_images)<tolerance){
				Im = i%ctx->num_images;
			}
		}
		if (Im!=-1){
			ctx->current_rec_num_mode++;
			for (int i=0; i<ctx->NOS; i++){
				// for delta m:
				IMAGE_COMPONENT(ctx->dImage,Im,ctx->NOS,i,0)=IMAGE_COMPONENT(ctx->dImage,Im,ctx->NOS,i,0)+(VEC_X(ctx->S,i)-VEC_X(ctx->t3S,i));
				IMAGE_COMPONENT(ctx->dImage,Im,ctx->NOS,i,1)=IMAGE_COMPONENT(ctx->dImage,Im,ctx->NOS,i,1)+(VEC_Y(ctx->S,i)-VEC_Y(ctx->t3S,i));
				IMAGE_COMPONENT(ctx->dImage,Im,ctx->NOS,i,2)=IMAGE_COMPONENT(ctx->dImage,Im,ctx->NOS,i,2)+(VEC_Z(ctx->S,i)-VEC_Z(ctx->t3S,i));
				// for m:
				IMAGE_COMPONENT(ctx->Image,Im,ctx->NOS,i,0) = IMAGE_COMPONENT(ctx->Image,Im,ctx->NOS,i,0) + VEC_X(ctx->S,i);
				IMAGE_COMPONENT(ctx->Image,Im,ctx->NOS,i,1) = IMAGE_COMPONENT(ctx->Image,Im,ctx->NOS,i,1) + VEC_Y(ctx->S,i);
				IMAGE_COMPONENT(ctx->Image,Im,ctx->NOS,i,2) = IMAGE_COMPONENT(ctx->Image,Im,ctx->NOS,i,2) + VEC_Z(ctx->S,i);
			}
			printf("tolerance %1.8f \n", tolerance);
			// printf("ABS(phase - i/(float)num_images)= %1.8f \n", ABS(phase - Im/(float)num_images));
			printf("Phase %1.8f (%d) (%1.8f)\n", phase,Im,Im/(float)ctx->num_images);
			printf("Image %d\n", ctx->current_rec_num_mode);
			if (ctx->current_rec_num_mode==ctx->num_images*ctx->rec_num_mode){
				ctx->Play=0;
				pthread_mutex_lock(&ctx->culc_mutex);
				ctx->ENGINE_MUTEX=STOP_REQUESTED;
				pthread_mutex_unlock(&ctx->culc_mutex);
				ctx->BextACModeRecording=0;
				ctx->BextACEnabled=0;
				ctx->current_rec_num_mode=0;
				char phase_output_directory[MAGNOOM_PATH_CAPACITY];
				pthread_mutex_lock(&ctx->record_mutex);
				magnoom_copy_path(phase_output_directory, sizeof(phase_output_directory),
					ctx->output_directory);
				pthread_mutex_unlock(&ctx->record_mutex);
				for (int j=0; j<ctx->num_images; j++){
					for (int i=0; i<ctx->NOS; i++){
					// for dm:
						IMAGE_COMPONENT(ctx->dImage,j,ctx->NOS,i,0)=IMAGE_COMPONENT(ctx->dImage,j,ctx->NOS,i,0)/ctx->rec_num_mode;
						IMAGE_COMPONENT(ctx->dImage,j,ctx->NOS,i,1)=IMAGE_COMPONENT(ctx->dImage,j,ctx->NOS,i,1)/ctx->rec_num_mode;
						IMAGE_COMPONENT(ctx->dImage,j,ctx->NOS,i,2)=IMAGE_COMPONENT(ctx->dImage,j,ctx->NOS,i,2)/ctx->rec_num_mode;
					// for m:
						IMAGE_COMPONENT(ctx->Image,j,ctx->NOS,i,0)=IMAGE_COMPONENT(ctx->Image,j,ctx->NOS,i,0)/ctx->rec_num_mode;
						IMAGE_COMPONENT(ctx->Image,j,ctx->NOS,i,1)=IMAGE_COMPONENT(ctx->Image,j,ctx->NOS,i,1)/ctx->rec_num_mode;
						IMAGE_COMPONENT(ctx->Image,j,ctx->NOS,i,2)=IMAGE_COMPONENT(ctx->Image,j,ctx->NOS,i,2)/ctx->rec_num_mode;
						float Norm=sqrt(IMAGE_COMPONENT(ctx->Image,j,ctx->NOS,i,0)*IMAGE_COMPONENT(ctx->Image,j,ctx->NOS,i,0)+IMAGE_COMPONENT(ctx->Image,j,ctx->NOS,i,1)*IMAGE_COMPONENT(ctx->Image,j,ctx->NOS,i,1)+IMAGE_COMPONENT(ctx->Image,j,ctx->NOS,i,2)*IMAGE_COMPONENT(ctx->Image,j,ctx->NOS,i,2));
						IMAGE_COMPONENT(ctx->Image,j,ctx->NOS,i,0)=IMAGE_COMPONENT(ctx->Image,j,ctx->NOS,i,0)/Norm;
						IMAGE_COMPONENT(ctx->Image,j,ctx->NOS,i,1)=IMAGE_COMPONENT(ctx->Image,j,ctx->NOS,i,1)/Norm;
						IMAGE_COMPONENT(ctx->Image,j,ctx->NOS,i,2)=IMAGE_COMPONENT(ctx->Image,j,ctx->NOS,i,2)/Norm;

					}
					/* only dm */
						// char ovf_filename[64] = "";
						// snprintf(ovf_filename,64,"dm%d.ovf",j);
						// Save_OVF_b8(Image_x[j], Image_y[j], Image_z[j], ovf_filename);
					/* only m */
						// char vtk_filename[64] = "";
						// snprintf(vtk_filename,64,"phase%d.vtk",j);
						// Save_VTK(Image_x[j], Image_y[j], Image_z[j],0, vtk_filename);
					/* m and dm */
						char vtk_filename[64] = "";
						char output_path[MAGNOOM_PATH_CAPACITY];
						snprintf(vtk_filename,64,"phase%d.vtk",j);
						if (magnoom_resolve_path(output_path, sizeof(output_path),
							phase_output_directory, vtk_filename)) {
							Save_VTK_6(ctx, &ctx->Image[(size_t)j*ctx->NOS*3], &ctx->dImage[(size_t)j*ctx->NOS*3], 0, output_path);
						}
					/*dTheta dPhi*/
						for (int i=0; i<ctx->NOS; i++){
						// get theta and phi for the equilibrium state:
							float T=acos(VEC_Z(ctx->t3S,i));
							float F=atan2(VEC_Y(ctx->t3S,i),VEC_X(ctx->t3S,i));
						// get spin i
							vec4 s, r;
							mat4x4 My, Mz;
							s[0]=(float)IMAGE_COMPONENT(ctx->Image,j,ctx->NOS,i,0);
							s[1]=(float)IMAGE_COMPONENT(ctx->Image,j,ctx->NOS,i,1);
							s[2]=(float)IMAGE_COMPONENT(ctx->Image,j,ctx->NOS,i,2);
							s[3]=0.f;
							mat4x4_identity(My);
							mat4x4_rotate_Y(My, My, T);
							mat4x4_identity(Mz);
							mat4x4_rotate_Z(Mz, Mz, F);
							mat4x4_mul_vec4(r, My, s);
							mat4x4_mul_vec4(s, Mz, r);
							IMAGE_COMPONENT(ctx->dImage,j,ctx->NOS,i,0)=s[0];
							IMAGE_COMPONENT(ctx->dImage,j,ctx->NOS,i,1)=s[1];
							IMAGE_COMPONENT(ctx->dImage,j,ctx->NOS,i,2)=s[2];
						}
							snprintf(vtk_filename,64,"dTdF%d.vtk",j);
							if (magnoom_resolve_path(output_path, sizeof(output_path),
								phase_output_directory, vtk_filename)) {
								Save_VTK(ctx, &ctx->dImage[(size_t)j*ctx->NOS*3],0, output_path);
							}
					}
			}
		}
	}
}
	

/* this function is run by the distinct thread */
void *CALC_THREAD(void *void_ptr)
{
    calc_thread_arg *arg = (calc_thread_arg *) void_ptr;
    int threadindex = arg->index;
    magnoom_ctx *ctx = arg->ctx;
    // printf("threadindex =%d\n", threadindex );
    int dNa=0;
    int dNb=0;
    int dNc=0;

    if (ctx->uABC[0]%THREADS_NUMBER==0){
    		dNa = ctx->uABC[0]/THREADS_NUMBER;
    }else{	dNa = (int)ctx->uABC[0]/THREADS_NUMBER+1;}

    if (ctx->uABC[1]%THREADS_NUMBER==0){
    		dNb = ctx->uABC[1]/THREADS_NUMBER;
    }else{	dNb = (int)ctx->uABC[1]/THREADS_NUMBER+1;}

    if (ctx->uABC[2]%THREADS_NUMBER==0){
    		dNc = ctx->uABC[2]/THREADS_NUMBER;
    }else{	dNc = (int)ctx->uABC[2]/THREADS_NUMBER+1;}

    int naini=0;
    int nafin=0;
    int nbini=0; 
    int nbfin=0;
    int ncini=0;
    int ncfin=0;

    if (ctx->uABC[0]>=ctx->uABC[1]&&ctx->uABC[0]>=ctx->uABC[2]){      //a-axis is the longest side of the box
    	naini = dNa*threadindex; 
    	if (dNa*(threadindex+1)<ctx->uABC[0]){
    		nafin = dNa*(threadindex+1);
    	}else{
    		nafin = ctx->uABC[0];
    	}
    	nbini = 0; nbfin = ctx->uABC[1];
    	// nbini = 1; nbfin = ctx->uABC[1];//metka
    	ncini = 0; ncfin = ctx->uABC[2];
    }else if (ctx->uABC[2]>=ctx->uABC[0]&&ctx->uABC[2]>=ctx->uABC[1]){//c-axis is the longest side of the box
    	ncini = dNc*threadindex; 
    	if ( (dNc*(threadindex+1))<ctx->uABC[2]){
    		ncfin = dNc*(threadindex+1);
    	}else{
    		ncfin = ctx->uABC[2];
    	}
    	naini = 0; nafin = ctx->uABC[0];	
    	nbini = 0; nbfin = ctx->uABC[1];
    }else if (ctx->uABC[1]>=ctx->uABC[0]&&ctx->uABC[1]>=ctx->uABC[2]){//b-axis is the longest side of the box
    	nbini = dNb*threadindex; 
    	if (dNb*(threadindex+1)<ctx->uABC[1]){
    		nbfin = dNb*(threadindex+1);
    	}else{
    		nbfin = ctx->uABC[1];
    	}
    	naini = 0; nafin = ctx->uABC[0];	
    	ncini = 0; ncfin = ctx->uABC[2];    	
    }

	while(!ctx->EngineShutdown)
	{
		int engine_state;
		do {
			pthread_mutex_lock(&ctx->culc_mutex);
			engine_state = ctx->ENGINE_MUTEX;
			pthread_mutex_unlock(&ctx->culc_mutex);
			if (engine_state == WAIT && !ctx->EngineShutdown)
				glfwSleep((double)ctx->SleepTime / 1000000.0);
		} while (engine_state == WAIT && !ctx->EngineShutdown);
		if (ctx->EngineShutdown) break;
		ctx->BextACScalar = EvaluateBextAC(ctx);
		
		switch (ctx->WhichIntegrationScheme){
			case HEUN: 
				// SimpleMinimizer(ctx->Sx,ctx->Sy,ctx->Sz,tSx,tSy,tSz,Heffx,Heffy,Heffz,RNx,RNy,RNz,ctx->NOS,damping,t_step,Temperature, threadindex,naini,nafin,nbini,nbfin,ncini,ncfin);
				StochasticLLG_Heun(ctx, threadindex, naini, nafin, nbini, nbfin, ncini, ncfin);
			break;

			case SIB: 
				StochasticLLG(ctx, threadindex, naini, nafin, nbini, nbfin, ncini, ncfin);
			break;

			case RK23: 
				StochasticLLG_RK23(ctx, threadindex, naini, nafin, nbini, nbfin, ncini, ncfin);
			break;

			case RK4:
				StochasticLLG_RK4(ctx, threadindex, naini, nafin, nbini, nbfin, ncini, ncfin);
			break;

			case RELAX: 	
				Relax(ctx, threadindex, naini, nafin, nbini, nbfin, ncini, ncfin);
			break;
		}
		if (threadindex==THREADS_NUMBER-1){ 
			
			//first thread opens the first (in) door in the next (second) thread
			sem_post(ctx->sem_in[(threadindex+1)%THREADS_NUMBER]);
			// first (in)door will be open from the last thread (first sem_post)
			sem_wait(ctx->sem_in[threadindex]);

			int nonfinite_spin = magnoom_first_nonfinite_spin(ctx);
			bool finite_step = nonfinite_spin < 0;
			if (!finite_step) {
				fprintf(stderr,
					"Solver stopped after iteration %d: spin %d became nonfinite (dt=%g).\n",
					ctx->ITERATION, nonfinite_spin, (double)ctx->t_step);
				ctx->Play = 0;
			}
			if (finite_step && ctx->Temperature > 0) GetFluctuations(ctx);
			
			ctx->MAX_TORQUE=0;
			for (int i=0;i<THREADS_NUMBER;i++){
				if (ctx->Max_torque[i] > ctx->MAX_TORQUE) ctx->MAX_TORQUE = ctx->Max_torque[i];
				ctx->Max_torque[i] = 0;
			}

			ctx->ITERATION++;
			if (ctx->ITERATION==ctx->Max_Numb_Iteration){
			    ctx->Play=0;
			    ctx->currentIteration=ctx->ITERATION;
			    for (int i=0;i<ctx->NOS;i++){
					VEC_X(ctx->bS,i)=VEC_X(ctx->S,i);
					VEC_Y(ctx->bS,i)=VEC_Y(ctx->S,i);
					VEC_Z(ctx->bS,i)=VEC_Z(ctx->S,i);
				}
				pthread_mutex_lock(&ctx->culc_mutex);
		            ctx->ENGINE_MUTEX=WAIT;
		            ctx->SleepTime=3000;
				pthread_mutex_unlock(&ctx->culc_mutex);
			}

			//normalize all spins every 100 iterations
			if (finite_step && ctx->WhichIntegrationScheme != RELAX && ctx->ITERATION%100==0)
			{
				// printf("ich bin hier!\n");
				for (int i=0;i<ctx->NOS;i++)
				{
					if (ctx->Kind[i]!=0){
					double absS = 1.0f/sqrt(VEC_X(ctx->S,i)*VEC_X(ctx->S,i)+VEC_Y(ctx->S,i)*VEC_Y(ctx->S,i)+VEC_Z(ctx->S,i)*VEC_Z(ctx->S,i));
					VEC_X(ctx->S,i) = VEC_X(ctx->S,i) * absS;
					VEC_Y(ctx->S,i) = VEC_Y(ctx->S,i) * absS;
					VEC_Z(ctx->S,i) = VEC_Z(ctx->S,i) * absS;
					}		
				}
			}

			//save to file if recording is on
			if (finite_step && ctx->Record!=0 && ctx->ITERATION%ctx->rec_iteration == 0){

				ctx->outputEtotal = GetTotalEnergyMoment(ctx);
				pthread_mutex_lock(&ctx->record_mutex);
				ctx->BigDataBank[0][ctx->recordsCounter] = (float)ctx->ITERATION;
				ctx->BigDataBank[1][ctx->recordsCounter] = ctx->outputMtotal[0]*ctx->iNOS;
				ctx->BigDataBank[2][ctx->recordsCounter] = ctx->outputMtotal[1]*ctx->iNOS;
				ctx->BigDataBank[3][ctx->recordsCounter] = ctx->outputMtotal[2]*ctx->iNOS;
				ctx->BigDataBank[4][ctx->recordsCounter] = ctx->outputEtotal;
				//metka test LLG
				// BigDataBank[0][recordsCounter] = ITERATION*t_step;
				// BigDataBank[1][recordsCounter] = bSx[0];
				// BigDataBank[2][recordsCounter] = bSy[0];
				// BigDataBank[3][recordsCounter] = bSz[1];
				// BigDataBank[4][recordsCounter] = bSx[1];
				// BigDataBank[5][recordsCounter] = bSy[1];
				// BigDataBank[6][recordsCounter] = bSz[1];				
				ctx->recordsCounter++;
				if (ctx->recordsCounter==100){
					magnoom_flush_records_locked(ctx);
				}
				pthread_mutex_unlock(&ctx->record_mutex);
			}

			if (finite_step) RecordBextACMode(ctx);

			if (ctx->DATA_TRANSFER_MUTEX==WAIT_DATA){
				for (int i=0;i<ctx->NOS;i++){
					VEC_X(ctx->bS,i)=VEC_X(ctx->S,i);
					VEC_Y(ctx->bS,i)=VEC_Y(ctx->S,i);
					VEC_Z(ctx->bS,i)=VEC_Z(ctx->S,i);
				}
				pthread_mutex_lock(&ctx->show_mutex);
					ctx->DATA_TRANSFER_MUTEX=TAKE_DATA;
					ctx->currentIteration=ctx->ITERATION;
				pthread_mutex_unlock(&ctx->show_mutex);	
			}
			// All workers leave together after completing the current iteration.
			if (ctx->EngineShutdownRequested) ctx->EngineShutdown = true;
			pthread_mutex_lock(&ctx->culc_mutex);
			if (!finite_step || ctx->ENGINE_MUTEX == STOP_REQUESTED)
				ctx->ENGINE_MUTEX = WAIT;
			pthread_mutex_unlock(&ctx->culc_mutex);

			// now it opens the second (out) door in the next (second) thread
			sem_post(ctx->sem_out[(threadindex+1)%THREADS_NUMBER]);
			// second (out)door will be open from the last thread (second sem_post)
			sem_wait(ctx->sem_out[threadindex]);
			pthread_mutex_lock(&ctx->culc_mutex);
			if (ctx->ENGINE_MUTEX == WAIT) ctx->EngineIdle = true;
			pthread_mutex_unlock(&ctx->culc_mutex);
		}else{
			//all other calculation threads
			sem_wait(ctx->sem_in[threadindex]);
			// first button which open the first door in the next (second) thread
			sem_post(ctx->sem_in[(threadindex+1)%THREADS_NUMBER]);
			// second door will be open from the last thread (second button)
			sem_wait(ctx->sem_out[threadindex]);
			// second button which open the second door in the next (second) thread
			sem_post(ctx->sem_out[(threadindex+1)%THREADS_NUMBER]);
		}
	}
/* the function must return something - NULL will do */
return NULL;
}
