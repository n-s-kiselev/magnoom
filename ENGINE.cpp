void
GetEffectiveField(	double* sx, double* sy, double* sz, 
					int numNeighbors, int* aidxBlock, int* nidxBlock, int* nidxGridA, int* nidxGridB, int* nidxGridC, int* shellIdx,
					float* Jij, float* Bij, float* Dij, float* VDMx, float* VDMy, float* VDMz, float* vku1, float ku1, float* vku2, float ku2, float kc, float* VHfield, float Hfield,
					double* heffx, double* heffy, double* heffz, int N, 
					int naini, 	int nafin,
					int nbini, 	int nbfin,
					int ncini, 	int ncfin)
{
	double tmp0;
	int Ip, I, J, K, L;
	int S;
	float Je, Bq, dx, dy, dz, DM;
	int i,j;
	int na1, Na = uABC[0];
	int nb1, Nb = uABC[1];
	int nc1, Nc = uABC[2];
	//single spin interactions (or potentila terms): Zeeman and Anizotropy:
	for (int Ip=0; Ip<AtomsPerBlock; Ip++)
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
					i = Ip + AtomsPerBlock * ( na + nb1 + nc1 );// index of spin "i"
					//H-field (Zeeman energy):
					heffx[i] = Hfield*VHfield[0]+AC_FIELD_ON*HacTime*VHac[0];
					heffy[i] = Hfield*VHfield[1]+AC_FIELD_ON*HacTime*VHac[1];
					heffz[i] = Hfield*VHfield[2]+AC_FIELD_ON*HacTime*VHac[2];
					//uniaxial anisotropy1:
					tmp0 = sx[i]*vku1[0] + sy[i]*vku1[1] + sz[i]*vku1[2];
					heffx[i]+= 2 * ku1 * vku1[0] * tmp0;
					heffy[i]+= 2 * ku1 * vku1[1] * tmp0;
					heffz[i]+= 2 * ku1 * vku1[2] * tmp0;
					//uniaxial anisotropy2:
					tmp0 = sx[i]*vku2[0] + sy[i]*vku2[1] + sz[i]*vku2[2];
					heffx[i]+= 2 * ku2 * vku2[0] * tmp0;
					heffy[i]+= 2 * ku2 * vku2[1] * tmp0;
					heffz[i]+= 2 * ku2 * vku2[2] * tmp0;
					//cubic anisotropy:
					heffx[i]+= 4 * kc * sx[i]*sx[i]*sx[i];
					heffy[i]+= 4 * kc * sy[i]*sy[i]*sy[i];
					heffz[i]+= 4 * kc * sz[i]*sz[i]*sz[i];			
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
			bc_c = 1 - (1-Boundary[2])*(( (2*Nc) + (nc+L) )/Nc )%2; // boundary condition along "c"
			for (int nb=nbini; nb<nbfin; nb++)
			{
				bc_b = 1 - (1-Boundary[1])*(( (2*Nb) + (nb+K) )/Nb )%2; // boundary condition along "b"
				for (int na=naini; na<nafin; na++)
				{
					bc_a = 1 - (1-Boundary[0])*(( (2*Na) + (na+J) )/Na )%2; // boundary condition along "a"
					bc_f = bc_a*bc_b*bc_c;
					na1 = na;
					nb1 = Na * nb;
					nc1 = Na * Nb * nc;
					i = Ip + AtomsPerBlock * ( na1 + nb1 + nc1 );// index of spin "i"
					na1 = (Na + na + J)%Na;
					nb1 = Na * ((Nb + nb + K)%Nb);
					nc1 = Na * Nb * ((Nc + nc + L)%Nc);					
					j = I  + AtomsPerBlock * ( na1 + nb1 + nc1 );// index of neighbouring spin "j"
					// Symmetric Heisenberg exchange:
					// heffx[i]+= Je * sx[j] * bc_f;
					// heffy[i]+= Je * sy[j] * bc_f;
					// heffz[i]+= Je * sz[j] * bc_f;
					// Bi-quadratic exchange:
					// tmp0 = (sx[i]*sx[j] + sy[i]*sy[j] + sz[i]*sz[j]) * bc_f;
					// heffx[i]+= 2 * Bq * sx[j] * tmp0;
					// heffy[i]+= 2 * Bq * sy[j] * tmp0;
					// heffz[i]+= 2 * Bq * sz[j] * tmp0;
					// Dzyaloshinskii-Moriya interaction (antisymmetric exchange):
					// heffx[i]+= DM*(sy[j]*dz - sz[j]*dy) * bc_f;
					// heffy[i]+= DM*(sz[j]*dx - sx[j]*dz) * bc_f;
					// heffz[i]+= DM*(sx[j]*dy - sy[j]*dx) * bc_f;	

					// a bit shorter version:
					tmp0 = Bq * 2 * (sx[i]*sx[j] + sy[i]*sy[j] + sz[i]*sz[j]);				
					heffx[i] = heffx[i] + (Je * sx[j] +  sx[j] * tmp0 + DM*(sy[j]*dz - sz[j]*dy) )* bc_f;
					heffy[i] = heffy[i] + (Je * sy[j] +  sy[j] * tmp0 + DM*(sz[j]*dx - sx[j]*dz) )* bc_f;
					heffz[i] = heffz[i] + (Je * sz[j] +  sz[j] * tmp0 + DM*(sx[j]*dy - sy[j]*dx) )* bc_f;
				}
			}
		}
	}
}

double
GetTotalEnergyMoment(	double* sx, double* sy, double* sz, double* Hx, double* Hy, double* Hz, double* Etot, double* Mtot, int N)
{
	double tmp0 = 0;
	Mtot[0] = 0;
	Mtot[1] = 0;
	Mtot[2] = 0;
	for (int i=0; i<N; i++)
	{
		Etot[i] = -Hx[i]*sx[i] - Hy[i]*sy[i] - Hz[i]*sz[i];
		Mtot[0] = Mtot[0] + sx[i];
		Mtot[1] = Mtot[1] + sy[i];
		Mtot[2] = Mtot[2] + sz[i];
		tmp0 = tmp0 + 0.5*Etot[i];
	}
	return tmp0;
}

double
GetTotalEnergyFerro(double sx, double sy, double sz, 
					int numNeighbors, int* aidxBlock, int* nidxBlock, int* nidxGridA, int* nidxGridB, int* nidxGridC, int* shellIdx,
					float* Jij, float* Bij, float* Dij, float* VDMx, float* VDMy, float* VDMz, float* vku1, float ku1, float* vku2, float ku2, float kc, float* VHfield, float Hfield,
					double* Etot, int N)
{
	double tmp0=1.0/sqrt(sx*sx+sy*sy+sz*sz);
	sx = sx/tmp0;
	sy = sy/tmp0;
	sz = sz/tmp0;
	tmp0=0;
	//single spin interactions (or potentila terms): Zeeman and Anizotropy:
	for (int i=0; i<N; i++)
	{
	//H-field (Zeeman energy):
	//Etot[i] =-Hfield*(VHfield[0]*sx+VHfield[1]*sy+VHfield[2]*sz);
	Etot[i] =-Hfield*((VHfield[0]+AC_FIELD_ON*HacTime*VHac[0])*sx+(VHfield[1]+AC_FIELD_ON*HacTime*VHac[1])*sy+(VHfield[2]+AC_FIELD_ON*HacTime*VHac[2])*sz);
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
	int Ip, I, J, K, L;
	int S;
	double Je, Bq, dx, dy, dz, DM;
	int i,j;
	int na1, Na = uABC[0];
	int nb1, Nb = uABC[1];
	int nc1, Nc = uABC[2];
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
			bc_c = 1 - (1-Boundary[2])*(( (2*Nc) + (nc+L) )/Nc )%2; // boundary condition along "c"
			for (int nb=0; nb<Nb; nb++)
			{
				bc_b = 1 - (1-Boundary[1])*(( (2*Nb) + (nb+K) )/Nb )%2; // boundary condition along "b"
				for (int na=0; na<Na; na++)
				{
					bc_a = 1 - (1-Boundary[0])*(( (2*Na) + (na+J) )/Na )%2; // boundary condition along "a"
					bc_f = bc_a*bc_b*bc_c;
					na1 = na;
					nb1 = Na * nb;
					nc1 = Na * Nb * nc;
					i = Ip + AtomsPerBlock * ( na1 + nb1 + nc1 );// index of spin "i"
					na1 = (Na + na + J)%Na;
					nb1 = Na * ((Nb + nb + K)%Nb);
					nc1 = Na * Nb * ((Nc + nc + L)%Nc);;
					j = I  + AtomsPerBlock * ( na1 + nb1 + nc1 );// index of neighbouring spin "j"
					tmp0 = sx*sx + sy*sy + sz*sz;
					Etot[i] = Etot[i] - bc_f*( Je * tmp0 + Bq * tmp0 * tmp0 + DM*(dx*(sy*sz - sz*sy) + dy*(sz*sx - sx*sz) + dz*(sx*sy - sy*sx)) );
				}
			}
		}
	}
	tmp0=0;
	for (i=0;i<NOS;i++)
	{
		tmp0 = tmp0 + Etot[i];
	}
	return tmp0;
}

double
GetTotalEnergy(	double* sx, double* sy, double* sz, 
					int numNeighbors, int* aidxBlock, int* nidxBlock, int* nidxGridA, int* nidxGridB, int* nidxGridC, int* shellIdx,
					float* Jij, float* Bij, float* Dij, float* VDMx, float* VDMy, float* VDMz, float* vku1, float ku1, float* vku2, float ku2, float kc, float* VHfield, float Hfield,
					double* Etot,double* Mtot, int N)
{
	double tmp0=0;
	Mtot[0] = 0;
	Mtot[1] = 0;
	Mtot[2] = 0;
	//single spin interactions (or potentila terms): Zeeman and Anizotropy:
	for (int i=0; i<N; i++)
	{
	//H-field (Zeeman energy):
	Etot[i] =-Hfield*((VHfield[0]+AC_FIELD_ON*HacTime*VHac[0])*sx[i]+(VHfield[1]+AC_FIELD_ON*HacTime*VHac[1])*sy[i]+(VHfield[2]+AC_FIELD_ON*HacTime*VHac[2])*sz[i]);
	//uniaxial anisotropy:
	tmp0 = sx[i]*vku1[0] + sy[i]*vku1[1] + sz[i]*vku1[2];
	Etot[i]-= ku1 * tmp0 * tmp0;
	//uniaxial anisotropy:
	tmp0 = sx[i]*vku2[0] + sy[i]*vku2[1] + sz[i]*vku2[2];
	Etot[i]-= ku2 * tmp0 * tmp0;
	//cubic anisotropy:
	Etot[i]-= kc * (sx[i]*sx[i]*sx[i]*sx[i] + sy[i]*sy[i]*sy[i]*sy[i] + sz[i]*sz[i]*sz[i]*sz[i]);	
	Mtot[0] = Mtot[0] + sx[i];
	Mtot[1] = Mtot[1] + sy[i];
	Mtot[2] = Mtot[2] + sz[i];
	}

	// pairwise spin interactions
	int Ip, I, J, K, L;
	int S;
	float Je, Bq, dx, dy, dz, DM;
	int i,j;
	int na1, Na = uABC[0];
	int nb1, Nb = uABC[1];
	int nc1, Nc = uABC[2];
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
			bc_c = 1 - (1-Boundary[2])*(( (2*Nc) + (nc+L) )/Nc )%2; // boundary condition along "c"
			for (int nb=0; nb<Nb; nb++)
			{
				bc_b = 1 - (1-Boundary[1])*(( (2*Nb) + (nb+K) )/Nb )%2; // boundary condition along "b"
				for (int na=0; na<Na; na++)
				{
					bc_a = 1 - (1-Boundary[0])*(( (2*Na) + (na+J) )/Na )%2; // boundary condition along "a"
					bc_f = bc_a*bc_b*bc_c;
					na1 = na;
					nb1 = Na * nb;
					nc1 = Na * Nb * nc;
					i = Ip + AtomsPerBlock * ( na1 + nb1 + nc1 );// index of spin "i"
					na1 = (Na + na + J)%Na;
					nb1 = Na * ((Nb + nb + K)%Nb);
					nc1 = Na * Nb * ((Nc + nc + L)%Nc);;
					j = I  + AtomsPerBlock * ( na1 + nb1 + nc1 );// index of neighbouring spin "j"
					//Symmetric Heisenberg exchange:
					// Etot[i]-= Je * (sx[i]*sx[j]+sy[i]*sy[j]+sz[i]*sz[j]);
					// //bi-quadratic exchange:
					// tmp0 = sx[i]*sx[j] + sy[i]*sy[j] + sz[i]*sz[j];
					// Etot[i]-= Bq * tmp0 * tmp0;
					// //Dzyaloshinskii-Moriya interaction (antisymmetric exchange):
					// Etot[i]-= DM*(dx*(sy[i]*sz[j] - sz[i]*sy[j]) + dy*(sz[i]*sx[j] - sx[i]*sz[j]) + dz*(sx[i]*sy[j] - sy[i]*sx[j]));
					tmp0 = sx[i]*sx[j] + sy[i]*sy[j] + sz[i]*sz[j];
					Etot[i] = Etot[i] - bc_f*( Je * tmp0 + Bq * tmp0 * tmp0 + DM*(dx*(sy[i]*sz[j] - sz[i]*sy[j]) + dy*(sz[i]*sx[j] - sx[i]*sz[j]) + dz*(sx[i]*sy[j] - sy[i]*sx[j])) );
				}
			}
		}
	}
	tmp0=0;
	for (i=0;i<NOS;i++)
	{
		tmp0 = tmp0 + Etot[i];
	}
	return tmp0;
}


void
GetFluctuations( float* rx, float* ry, float* rz, int N )
{
	for (int i=0; i<N; i++)
	{
	rx[i] = 2.0 * (0.5 - rand() / (float)RAND_MAX);
	ry[i] = 2.0 * (0.5 - rand() / (float)RAND_MAX);
	rz[i] = 2.0 * (0.5 - rand() / (float)RAND_MAX);
	}
}

void
StochasticLLG(	double* inx,		double* iny,		double* inz,		// input vector field
				double* tnx,		double* tny,		double* tnz,		// temporal storage
				double* heffx,	double* heffy,	double* heffz,	// effective field
				float* rx,		float* ry,		float* rz,		// random numbers 
				int nos,		float alpha, 	float h,		// number of spins, damping, time step
				float temperature, 	int thread,					
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
	double rh=sqrt(h);// h = \delta t
	double nx, ny, nz;// components of the unit vector
	double Hx, Hy, Hz;// components of the effective field
	float Cx, Cy, Cz;// spin-torque term
	float Rx, Ry, Rz;// random variables - thermal fluctuations
	double ax, ay, az;// temp
	double Ax, Ay, Az;// total matrix
	double detMi;	 // detMi = 1/detM
	double D = sqrt(2.0 * alpha / (1.0 + alpha * alpha) * temperature);
	double Alpha_d = alpha / (1.0 + alpha * alpha * Precession);
	double Alpha_p = 1.0f / (1.0 + alpha * alpha);

	//if (temperature>0) GetFluctuations( rx, ry, rz, nos );
	//electric DC current vector (VCu) and density (Cu)
	Cx = VCu[0] * Cu;
	Cy = VCu[1] * Cu;
	Cz = VCu[2] * Cu;
	GetEffectiveField( 	inx, iny, inz, 
						NeighborPairs, AIdxBlock, NIdxBlock, NIdxGridA, NIdxGridB, NIdxGridC, SIdx,
						Jij, Bij, Dij, VDMx, VDMy, VDMz, VKu1, Ku1, VKu2, Ku2, Kc, VHf, Hf, Heffx, Heffy, Heffz, NOS,
						naini, nafin, nbini, nbfin, ncini, ncfin);
	//prediction step of midpoint solver:
	int Na = uABC[0];
	int Nb = uABC[1];
	int nb1, nc1;
	int i;
	for (int Ip=0; Ip<AtomsPerBlock; Ip++)
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
					i = Ip + AtomsPerBlock * ( na + nb1 + nc1 );// index of spin "i"
					nx = inx[i];	ny = iny[i];	nz = inz[i];
					Hx = Heffx[i];	Hy = Heffy[i];	Hz = Heffz[i];
					Rx = rx[i];		Ry = ry[i];		Rz = rz[i];

					// deterministic terms of Landau–Lifshitz equation:
					Ax = 0.5f * h * ( - Alpha_p * Precession * Hx - Alpha_d * (ny * Hz - nz * Hy) );
					Ay = 0.5f * h * ( - Alpha_p * Precession * Hy - Alpha_d * (nz * Hx - nx * Hz) );
					Az = 0.5f * h * ( - Alpha_p * Precession * Hz - Alpha_d * (nx * Hy - ny * Hx) );

					// Spin-torque term
					if (Cu!=0) {
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
					// tnx[i] = (ax*(1+Ax*Ax)+ay*(Ax*Ay+Az)+az*(Ax*Az-Ay))*detMi
					// tny[i] = (ax*(Ay*Ax-Az)+ay*(1+Ay*Ay)+az*(Ay*Az+Ax))*detMi
					// tnz[i] = (ax*(Az*Ax+Ay)+ay*(Az*Ay-Ax)+az*(1+Az*Az))*detMi
					// tnx[i] = ( tnx[i] + nx ) / 2;
					// tny[i] = ( tnx[i] + nx ) / 2;
					// tnz[i] = ( tnx[i] + nx ) / 2;

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

					tnx[i] = nx * 0.5f;
					tny[i] = ny * 0.5f;
					tnz[i] = nz * 0.5f;
				}
			}
		}
	}

	GetEffectiveField( 	tnx, tny, tnz, 
						NeighborPairs, AIdxBlock, NIdxBlock, NIdxGridA, NIdxGridB, NIdxGridC, SIdx,
						Jij, Bij, Dij, VDMx, VDMy, VDMz, VKu1, Ku1, VKu2, Ku2, Kc, VHf, Hf, Heffx, Heffy, Heffz, NOS,
						naini, nafin, nbini, nbfin, ncini, ncfin);

	//final step of midpoint solver:
	for (int Ip=0; Ip<AtomsPerBlock; Ip++)
	{
		for (int nc=ncini; nc<ncfin; nc++)//nc(a,b)=neghbor in the direction of c(a,b)-vector
		{
			nc1 = Na * Nb * nc;
			for (int nb=nbini; nb<nbfin; nb++)
			{
				nb1 = Na * nb;
				for (int na=naini; na<nafin; na++)
				{
					i = Ip + AtomsPerBlock * ( na + nb1 + nc1 );// index of spin "i"
					nx = tnx[i];	ny = tny[i];	nz = tnz[i];	// <-- compare this line to corresponding one in prediction step
					Hx = Heffx[i];	Hy = Heffy[i];	Hz = Heffz[i];	// <-- they are new values for heff
					Rx = rx[i];		Ry = ry[i];		Rz = rz[i];		// <-- they are the same values as in prediction step
					// deterministic terms of Landau–Lifshitz equation:
					Ax = 0.5f * h * ( - Alpha_p * Precession * Hx - Alpha_d * (ny * Hz - nz * Hy) );
					Ay = 0.5f * h * ( - Alpha_p * Precession * Hy - Alpha_d * (nz * Hx - nx * Hz) );
					Az = 0.5f * h * ( - Alpha_p * Precession * Hz - Alpha_d * (nx * Hy - ny * Hx) );
					// Spin-torque term
					Ax = Ax + 0.5f * h * ( - Alpha_d * Cx + (ny * Cz - nz * Cy) ); //pay attention to the signe and factors
					Ay = Ay + 0.5f * h * ( - Alpha_d * Cy + (nz * Cx - nx * Cz) );
					Az = Az + 0.5f * h * ( - Alpha_d * Cz + (nx * Cy - ny * Cx) );
					// stochastic terms of Landau–Lifshitz equation:
					Ax = Ax + 0.5f * rh * D * ( - Alpha_p * Rx - Alpha_d * (ny * Rz - nz * Ry) );
					Ay = Ay + 0.5f * rh * D * ( - Alpha_p * Ry - Alpha_d * (nz * Rx - nx * Rz) );
					Az = Az + 0.5f * rh * D * ( - Alpha_p * Rz - Alpha_d * (nx * Ry - ny * Rx) );

					nx = inx[i];	ny = iny[i];	nz = inz[i];	//<-- pay attention to this line!

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
					
					inx[i] = ( ax * (1. + Hx) + ay * (Rz + Az) + az * (Ry - Ay) ) * detMi;// <-- back to the array of spins new values
					iny[i] = ( ax * (Rz - Az) + ay * (1. + Hy) + az * (Rx + Ax) ) * detMi;
					inz[i] = ( ax * (Ry + Ay) + ay * (Rx - Ax) + az * (1. + Hz) ) * detMi;

					// find max torque:
					detMi = Heffx[i]*inx[i]+Heffy[i]*iny[i]+Heffz[i]*inz[i];
					Hx = Heffx[i]-detMi*inx[i];
					Hy = Heffy[i]-detMi*iny[i];	
					Hz = Heffz[i]-detMi*inz[i];
					detMi = sqrt(Hx*Hx + Hy*Hy + Hz*Hz);
					if (detMi > Max_torque[thread]) Max_torque[thread] = detMi;
					bSx[i]=inx[i];
					bSy[i]=iny[i];
					bSz[i]=inz[i];
				}
			}
		}
	}	
}


void
StochasticLLG_Heun(	double* inx,		double* iny,		double* inz,		// input vector field
				double* tnx,		double* tny,		double* tnz,		// temporal storage
				double* heffx,	double* heffy,	double* heffz,	// effective field
				float* rx,		float* ry,		float* rz,		// random numbers 
				int nos,		float alpha, 	float h,		// number of spins, damping, time step
				float temperature, 	int thread,					
				int naini, 	int nafin,
				int nbini, 	int nbfin,
				int ncini, 	int ncfin)

{
	double rh=sqrt(h);// h = \delta t
	double nx, ny, nz;// components of the unit vector
	double Hx, Hy, Hz;// components of the effective field
	float Cx, Cy, Cz;// spin-torque term
	float Rx, Ry, Rz;// random variables - thermal fluctuations
	double Fx, Fy, Fz;// total matrix
	double detMi;	 // detMi = 1/detM
	double D = sqrt(2.0 * alpha / (1.0 + alpha * alpha) * temperature);
	double Alpha_d = alpha / (1.0 + alpha * alpha * Precession);
	double Alpha_p = 1.0f / (1.0 + alpha * alpha);

	//if (temperature>0) GetFluctuations( rx, ry, rz, nos );
	//electric DC current vector (VCu) and density (Cu)
	Cx = VCu[0] * Cu;
	Cy = VCu[1] * Cu;
	Cz = VCu[2] * Cu;
	GetEffectiveField( 	inx, iny, inz, 
						NeighborPairs, AIdxBlock, NIdxBlock, NIdxGridA, NIdxGridB, NIdxGridC, SIdx,
						Jij, Bij, Dij, VDMx, VDMy, VDMz, VKu1, Ku1, VKu2, Ku2, Kc, VHf, Hf, Heffx, Heffy, Heffz, NOS,
						naini, nafin, nbini, nbfin, ncini, ncfin);
	//Predictor step:
	int Na = uABC[0];
	int Nb = uABC[1];
	int nb1, nc1;
	int i;
	for (int Ip=0; Ip<AtomsPerBlock; Ip++)
	{
		for (int nc=ncini; nc<ncfin; nc++)//nc(or na or nb) is a neghbor in the direction of c (or a or b)-vector
		{
			nc1 = Na * Nb * nc;
			for (int nb=nbini; nb<nbfin; nb++)
			{
				nb1 = Na * nb;
				for (int na=naini; na<nafin; na++)
				{
					i = Ip + AtomsPerBlock * ( na + nb1 + nc1 );// index of spin "i"
					nx = inx[i];	ny = iny[i];	nz = inz[i];
					Hx = Heffx[i];	Hy = Heffy[i];	Hz = Heffz[i];
					Rx = rx[i];		Ry = ry[i];		Rz = rz[i];

					// deterministic terms of Landau–Lifshitz equation:
					Fx = Alpha_p * Hx + Alpha_d * (ny * Hz - nz * Hy);
					Fy = Alpha_p * Hy + Alpha_d * (nz * Hx - nx * Hz);
					Fz = Alpha_p * Hz + Alpha_d * (nx * Hy - ny * Hx);

					// Spin-torque term
					if (Cu!=0) {
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
					
					tnx[i] = nx - h * (ny * Fz - nz * Fy);
					tny[i] = ny - h * (nz * Fx - nx * Fz);
					tnz[i] = nz - h * (nx * Fy - ny * Fx);
				}
			}
		}
	}

	GetEffectiveField( 	tnx, tny, tnz, 
						NeighborPairs, AIdxBlock, NIdxBlock, NIdxGridA, NIdxGridB, NIdxGridC, SIdx,
						Jij, Bij, Dij, VDMx, VDMy, VDMz, VKu1, Ku1, VKu2, Ku2, Kc, VHf, Hf, Heffx, Heffy, Heffz, NOS,
						naini, nafin, nbini, nbfin, ncini, ncfin);

	//corrector step:
	for (int Ip=0; Ip<AtomsPerBlock; Ip++)
	{
		for (int nc=ncini; nc<ncfin; nc++)//nc(a,b)=neghbor in the direction of c(a,b)-vector
		{
			nc1 = Na * Nb * nc;
			for (int nb=nbini; nb<nbfin; nb++)
			{
				nb1 = Na * nb;
				for (int na=naini; na<nafin; na++)
				{
					i = Ip + AtomsPerBlock * ( na + nb1 + nc1 );	// index of spin "i"
					nx = tnx[i];	ny = tny[i];	nz = tnz[i];	// <-- compare this line to corresponding one in prediction step
					Hx = Heffx[i];	Hy = Heffy[i];	Hz = Heffz[i];	// <-- they are new values for heff
					Rx = rx[i];		Ry = ry[i];		Rz = rz[i];		// <-- they are the same values as in prediction step

					// deterministic terms of Landau–Lifshitz equation:
					Fx = Alpha_p * Hx + Alpha_d * (ny * Hz - nz * Hy);
					Fy = Alpha_p * Hy + Alpha_d * (nz * Hx - nx * Hz);
					Fz = Alpha_p * Hz + Alpha_d * (nx * Hy - ny * Hx);

					// Spin-torque term
					if (Cu!=0) {
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
					
					nx = 0.5 * ( nx + inx[i] - h * (ny * Fz - nz * Fy) );
					ny = 0.5 * ( ny + iny[i] - h * (nz * Fx - nx * Fz) );
					nz = 0.5 * ( nz + inz[i] - h * (nx * Fy - ny * Fx) );

					//normalize spin
					detMi = 1.0 / sqrt(nx*nx + ny*ny + nz*nz);
					inx[i] = nx * detMi;
					iny[i] = ny * detMi;
					inz[i] = nz * detMi;

					//find max torque
					detMi = Heffx[i]*inx[i]+Heffy[i]*iny[i]+Heffz[i]*inz[i];
					Hx = Heffx[i]-detMi*inx[i];
					Hy = Heffy[i]-detMi*iny[i];	
					Hz = Heffz[i]-detMi*inz[i];
					detMi = sqrt(Hx*Hx + Hy*Hy + Hz*Hz);
					if (detMi > Max_torque[thread]) Max_torque[thread] = detMi;
					bSx[i]=inx[i];
					bSy[i]=iny[i];
					bSz[i]=inz[i];
				}
			}
		}
	}	
}



void
StochasticLLG_RK23(	double* inx,		double* iny,		double* inz,		// input vector field
				double* tnx,		double* tny,		double* tnz,		// temporal storage
				double* heffx,	double* heffy,	double* heffz,	// effective field
				float* rx,		float* ry,		float* rz,		// random numbers 
				int nos,		float alpha, 	float h,		// number of spins, damping, time step
				float temperature, 	int thread,					
				int naini, 	int nafin,
				int nbini, 	int nbfin,
				int ncini, 	int ncfin)

{
	double rh=sqrt(h);// h = \delta t
	double nx, ny, nz;// components of the unit vector
	double Hx, Hy, Hz;// components of the effective field
	float Cx, Cy, Cz;// spin-torque term
	float Rx, Ry, Rz;// random variables - thermal fluctuations
	double Fx, Fy, Fz;// total matrix
	double detMi;	 // detMi = 1/detM
	double D = sqrt(2.0 * alpha / (1.0 + alpha * alpha) * temperature);
	double Alpha_d = alpha / (1.0 + alpha * alpha * Precession);
	double Alpha_p = 1.0f / (1.0 + alpha * alpha);

	//if (temperature>0) GetFluctuations( rx, ry, rz, nos );
	//electric DC current vector (VCu) and density (Cu)
	Cx = VCu[0] * Cu;
	Cy = VCu[1] * Cu;
	Cz = VCu[2] * Cu;
	GetEffectiveField( 	inx, iny, inz, 
						NeighborPairs, AIdxBlock, NIdxBlock, NIdxGridA, NIdxGridB, NIdxGridC, SIdx,
						Jij, Bij, Dij, VDMx, VDMy, VDMz, VKu1, Ku1, VKu2, Ku2, Kc, VHf, Hf, Heffx, Heffy, Heffz, NOS,
						naini, nafin, nbini, nbfin, ncini, ncfin);
	//k1:
	int Na = uABC[0];
	int Nb = uABC[1];
	int nb1, nc1;
	int i;
	for (int Ip=0; Ip<AtomsPerBlock; Ip++)
	{
		for (int nc=ncini; nc<ncfin; nc++)//nc(or na or nb) is a neghbor in the direction of c (or a or b)-vector
		{
			nc1 = Na * Nb * nc;
			for (int nb=nbini; nb<nbfin; nb++)
			{
				nb1 = Na * nb;
				for (int na=naini; na<nafin; na++)
				{
					i = Ip + AtomsPerBlock * ( na + nb1 + nc1 );// index of spin "i"
					nx = inx[i];	ny = iny[i];	nz = inz[i];
					Hx = Heffx[i];	Hy = Heffy[i];	Hz = Heffz[i];
					Rx = rx[i];		Ry = ry[i];		Rz = rz[i];

					// deterministic terms of Landau–Lifshitz equation:
					Fx = Alpha_p * Hx + Alpha_d * (ny * Hz - nz * Hy);
					Fy = Alpha_p * Hy + Alpha_d * (nz * Hx - nx * Hz);
					Fz = Alpha_p * Hz + Alpha_d * (nx * Hy - ny * Hx);

					// Spin-torque term
					if (Cu!=0) {
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
					
					tnx[i] = nx - 0.5*h * (ny * Fz - nz * Fy);
					tny[i] = ny - 0.5*h * (nz * Fx - nx * Fz);
					tnz[i] = nz - 0.5*h * (nx * Fy - ny * Fx);
				}
			}
		}
	}

	GetEffectiveField( 	tnx, tny, tnz, 
						NeighborPairs, AIdxBlock, NIdxBlock, NIdxGridA, NIdxGridB, NIdxGridC, SIdx,
						Jij, Bij, Dij, VDMx, VDMy, VDMz, VKu1, Ku1, VKu2, Ku2, Kc, VHf, Hf, Heffx, Heffy, Heffz, NOS,
						naini, nafin, nbini, nbfin, ncini, ncfin);

	//corrector step:
	for (int Ip=0; Ip<AtomsPerBlock; Ip++)
	{
		for (int nc=ncini; nc<ncfin; nc++)//nc(a,b)=neghbor in the direction of c(a,b)-vector
		{
			nc1 = Na * Nb * nc;
			for (int nb=nbini; nb<nbfin; nb++)
			{
				nb1 = Na * nb;
				for (int na=naini; na<nafin; na++)
				{
					i = Ip + AtomsPerBlock * ( na + nb1 + nc1 );	// index of spin "i"
					nx = tnx[i];	ny = tny[i];	nz = tnz[i];	// <-- compare this line to corresponding one in prediction step
					Hx = Heffx[i];	Hy = Heffy[i];	Hz = Heffz[i];	// <-- they are new values for heff
					Rx = rx[i];		Ry = ry[i];		Rz = rz[i];		// <-- they are the same values as in prediction step

					// deterministic terms of Landau–Lifshitz equation:
					Fx = Alpha_p * Hx + Alpha_d * (ny * Hz - nz * Hy);
					Fy = Alpha_p * Hy + Alpha_d * (nz * Hx - nx * Hz);
					Fz = Alpha_p * Hz + Alpha_d * (nx * Hy - ny * Hx);

					// Spin-torque term
					if (Cu!=0) {
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
					Cx = inx[i] - h * (ny * Fz - nz * Fy);
					Cy = iny[i] - h * (nz * Fx - nx * Fz);
					Cz = inz[i] - h * (nx * Fy - ny * Fx);

					//normalize spin
					detMi = 1.0 / sqrt(Cx*Cx + Cy*Cy + Cz*Cz);
					inx[i] = Cx * detMi;
					iny[i] = Cy * detMi;
					inz[i] = Cz * detMi;

					//find max torque
					detMi = Heffx[i]*inx[i]+Heffy[i]*iny[i]+Heffz[i]*inz[i];
					Hx = Heffx[i]-detMi*inx[i];
					Hy = Heffy[i]-detMi*iny[i];	
					Hz = Heffz[i]-detMi*inz[i];
					detMi = sqrt(Hx*Hx + Hy*Hy + Hz*Hz);
					if (detMi > Max_torque[thread]) Max_torque[thread] = detMi;
					bSx[i]=inx[i];
					bSy[i]=iny[i];
					bSz[i]=inz[i];
				}
			}
		}
	}	
}


void
StochasticLLG_RK45(	double* inx,		double* iny,		double* inz,		// input vector field
				double* tnx,		double* tny,		double* tnz,		// temporal storage
				double* heffx,	double* heffy,	double* heffz,	// effective field
				float* rx,		float* ry,		float* rz,		// random numbers 
				int nos,		float alpha, 	float h,		// number of spins, damping, time step
				float temperature, 	int thread,					
				int naini, 	int nafin,
				int nbini, 	int nbfin,
				int ncini, 	int ncfin)
{
	double rh=sqrt(h);// h = \delta t
	double nx, ny, nz;// components of the unit vector
	double Hx, Hy, Hz;// components of the effective field
	float Cx, Cy, Cz;// spin-torque term
	float Rx, Ry, Rz;// random variables - thermal fluctuations
	double Fx, Fy, Fz;// total matrix
	double detMi;	 // detMi = 1/detM
	double D = sqrt(2.0 * alpha / (1.0 + alpha * alpha) * temperature);
	double Alpha_d = alpha / (1.0 + alpha * alpha * Precession);
	double Alpha_p = 1.0f / (1.0 + alpha * alpha);

	//if (temperature>0) GetFluctuations( rx, ry, rz, nos );
	//electric DC current vector (VCu) and density (Cu)
	Cx = VCu[0] * Cu;
	Cy = VCu[1] * Cu;
	Cz = VCu[2] * Cu;
	GetEffectiveField( 	inx, iny, inz, 
						NeighborPairs, AIdxBlock, NIdxBlock, NIdxGridA, NIdxGridB, NIdxGridC, SIdx,
						Jij, Bij, Dij, VDMx, VDMy, VDMz, VKu1, Ku1, VKu2, Ku2, Kc, VHf, Hf, Heffx, Heffy, Heffz, NOS,
						naini, nafin, nbini, nbfin, ncini, ncfin);
	//Predictor step:
	int Na = uABC[0];
	int Nb = uABC[1];
	int nb1, nc1;
	int i;

	//k1:
	for (int Ip=0; Ip<AtomsPerBlock; Ip++)
	{
		for (int nc=ncini; nc<ncfin; nc++)//nc(or na or nb) is a neghbor in the direction of c (or a or b)-vector
		{
			nc1 = Na * Nb * nc;
			for (int nb=nbini; nb<nbfin; nb++)
			{
				nb1 = Na * nb;
				for (int na=naini; na<nafin; na++)
				{
					i = Ip + AtomsPerBlock * ( na + nb1 + nc1 );// index of spin "i"
					nx = inx[i];	ny = iny[i];	nz = inz[i];
					Hx = Heffx[i];	Hy = Heffy[i];	Hz = Heffz[i];
					Rx = rx[i];		Ry = ry[i];		Rz = rz[i];

					// deterministic terms of Landau–Lifshitz equation:
					Fx = Alpha_p * Hx + Alpha_d * (ny * Hz - nz * Hy);
					Fy = Alpha_p * Hy + Alpha_d * (nz * Hx - nx * Hz);
					Fz = Alpha_p * Hz + Alpha_d * (nx * Hy - ny * Hx);

					// Spin-torque term
					if (Cu!=0) {
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
					t2Sx[i] = Cx/6.0;
					t2Sy[i] = Cy/6.0;
					t2Sz[i] = Cz/6.0;
					//y_n+k1/2 will be used on the next step
					tnx[i] = nx + Cx*0.5;
					tny[i] = ny + Cy*0.5;
					tnz[i] = nz + Cz*0.5;				
				}
			}
		}
	}
	//Heff(y_n+k1/2):
	GetEffectiveField( 	tnx, tny, tnz, 
						NeighborPairs, AIdxBlock, NIdxBlock, NIdxGridA, NIdxGridB, NIdxGridC, SIdx,
						Jij, Bij, Dij, VDMx, VDMy, VDMz, VKu1, Ku1, VKu2, Ku2, Kc, VHf, Hf, Heffx, Heffy, Heffz, NOS,
						naini, nafin, nbini, nbfin, ncini, ncfin);

	//k2:
	for (int Ip=0; Ip<AtomsPerBlock; Ip++)
	{
		for (int nc=ncini; nc<ncfin; nc++)//nc(or na or nb) is a neghbor in the direction of c (or a or b)-vector
		{
			nc1 = Na * Nb * nc;
			for (int nb=nbini; nb<nbfin; nb++)
			{
				nb1 = Na * nb;
				for (int na=naini; na<nafin; na++)
				{
					i = Ip + AtomsPerBlock * ( na + nb1 + nc1 );// index of spin "i"
					nx = tnx[i];	ny = tny[i];	nz = tnz[i];
					Hx = Heffx[i];	Hy = Heffy[i];	Hz = Heffz[i];
					Rx = rx[i];		Ry = ry[i];		Rz = rz[i];

					// deterministic terms of Landau–Lifshitz equation:
					Fx = Alpha_p * Hx + Alpha_d * (ny * Hz - nz * Hy);
					Fy = Alpha_p * Hy + Alpha_d * (nz * Hx - nx * Hz);
					Fz = Alpha_p * Hz + Alpha_d * (nx * Hy - ny * Hx);

					// Spin-torque term
					if (Cu!=0) {
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
					t2Sx[i]+= Cx/3.0;
					t2Sy[i]+= Cy/3.0;
					t2Sz[i]+= Cz/3.0;
					//y_n+k2/2 will be used on the next step
					tnx[i] = inx[i] + Cx*0.5;
					tny[i] = iny[i] + Cy*0.5;
					tnz[i] = inz[i] + Cz*0.5;						
				}
			}
		}
	}
	//Heff(y_n+k2/2):
	GetEffectiveField( 	tnx, tny, tnz, 
						NeighborPairs, AIdxBlock, NIdxBlock, NIdxGridA, NIdxGridB, NIdxGridC, SIdx,
						Jij, Bij, Dij, VDMx, VDMy, VDMz, VKu1, Ku1, VKu2, Ku2, Kc, VHf, Hf, Heffx, Heffy, Heffz, NOS,
						naini, nafin, nbini, nbfin, ncini, ncfin);	
	//k3:
	for (int Ip=0; Ip<AtomsPerBlock; Ip++)
	{
		for (int nc=ncini; nc<ncfin; nc++)//nc(or na or nb) is a neghbor in the direction of c (or a or b)-vector
		{
			nc1 = Na * Nb * nc;
			for (int nb=nbini; nb<nbfin; nb++)
			{
				nb1 = Na * nb;
				for (int na=naini; na<nafin; na++)
				{
					i = Ip + AtomsPerBlock * ( na + nb1 + nc1 );// index of spin "i"
					nx = tnx[i];	ny = tny[i];	nz = tnz[i];
					Hx = Heffx[i];	Hy = Heffy[i];	Hz = Heffz[i];
					Rx = rx[i];		Ry = ry[i];		Rz = rz[i];

					// deterministic terms of Landau–Lifshitz equation:
					Fx = Alpha_p * Hx + Alpha_d * (ny * Hz - nz * Hy);
					Fy = Alpha_p * Hy + Alpha_d * (nz * Hx - nx * Hz);
					Fz = Alpha_p * Hz + Alpha_d * (nx * Hy - ny * Hx);

					// Spin-torque term
					if (Cu!=0) {
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
					t2Sx[i]+= Cx/3.0;
					t2Sy[i]+= Cy/3.0;
					t2Sz[i]+= Cz/3.0;
					//y_n+k3 will be used on the next step
					tnx[i] = inx[i] + Cx;
					tny[i] = iny[i] + Cy;
					tnz[i] = inz[i] + Cz;							
				}
			}
		}
	}
	//Heff(y_n+k3):
	GetEffectiveField( 	tnx, tny, tnz, 
						NeighborPairs, AIdxBlock, NIdxBlock, NIdxGridA, NIdxGridB, NIdxGridC, SIdx,
						Jij, Bij, Dij, VDMx, VDMy, VDMz, VKu1, Ku1, VKu2, Ku2, Kc, VHf, Hf, Heffx, Heffy, Heffz, NOS,
						naini, nafin, nbini, nbfin, ncini, ncfin);	
	//k4:
	for (int Ip=0; Ip<AtomsPerBlock; Ip++)
	{
		for (int nc=ncini; nc<ncfin; nc++)//nc(or na or nb) is a neghbor in the direction of c (or a or b)-vector
		{
			nc1 = Na * Nb * nc;
			for (int nb=nbini; nb<nbfin; nb++)
			{
				nb1 = Na * nb;
				for (int na=naini; na<nafin; na++)
				{
					i = Ip + AtomsPerBlock * ( na + nb1 + nc1 );// index of spin "i"
					nx = tnx[i];	ny = tny[i];	nz = tnz[i];
					Hx = Heffx[i];	Hy = Heffy[i];	Hz = Heffz[i];
					Rx = rx[i];		Ry = ry[i];		Rz = rz[i];

					// deterministic terms of Landau–Lifshitz equation:
					Fx = Alpha_p * Hx + Alpha_d * (ny * Hz - nz * Hy);
					Fy = Alpha_p * Hy + Alpha_d * (nz * Hx - nx * Hz);
					Fz = Alpha_p * Hz + Alpha_d * (nx * Hy - ny * Hx);

					// Spin-torque term
					if (Cu!=0) {
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
					t2Sx[i]+= Cx/6.0;
					t2Sy[i]+= Cy/6.0;
					t2Sz[i]+= Cz/6.0;
					//y_{n+1}=y_n+k1/6+k2/3+k3/3+k4/6 - final step:
					inx[i]+= t2Sx[i];
					iny[i]+= t2Sy[i];
					inz[i]+= t2Sz[i];
					nx = inx[i];	ny = iny[i];	nz = inz[i];
					//normalize spin
					detMi = 1.0 / sqrt(nx*nx + ny*ny + nz*nz);
					inx[i] = nx * detMi;
					iny[i] = ny * detMi;
					inz[i] = nz * detMi;

					//find max torque
					detMi = Heffx[i]*inx[i]+Heffy[i]*iny[i]+Heffz[i]*inz[i];
					Hx = Heffx[i]-detMi*inx[i];
					Hy = Heffy[i]-detMi*iny[i];	
					Hz = Heffz[i]-detMi*inz[i];
					detMi = sqrt(Hx*Hx + Hy*Hy + Hz*Hz);
					if (detMi > Max_torque[thread]) Max_torque[thread] = detMi;	
					bSx[i]=inx[i];
					bSy[i]=iny[i];
					bSz[i]=inz[i];					
				}
			}
		}
	}
}


void
SimpleMinimizer(double* inx,		double* iny,		double* inz,		// input vector field
				double* tnx,		double* tny,		double* tnz,		// temporal storage
				double* heffx,	double* heffy,	double* heffz,	// effective field
				float* rx,		float* ry,		float* rz,		// random numbers 
				int nos,		float alpha, 	float h,		// number of spins, damping, time step
				float temperature)
{	

}


float GetACfield(int type)
{
	float R=0;
	float temp;
	switch (type){
		case SIN_FIELD:
			R = AC_FIELD_ON*Hac*sin(Omega_dc*t_step*ITERATION);
		break;
		case GAUSSIAN_FIELD:
			temp = (t_step*ITERATION-t_offset)/GPulseWidth;
			R = AC_FIELD_ON*Hac*exp(-0.5*temp*temp );
		break;
	}
		Bac[0] = R * VHac[0];
		Bac[1] = R * VHac[1];
		Bac[2] = R * VHac[2];
	return R;
}
	

/* this function is run by the distinct thread */
void *CALC_THREAD(void *void_ptr)
{
    int threadindex = *((int *) void_ptr);
    //printf("threadindex =%d\n", threadindex );
    int dNa=0;
    int dNb=0;
    int dNc=0;

    if (uABC[0]%THREADS_NUMBER==0){
    		dNa = uABC[0]/THREADS_NUMBER;
    }else{	dNa = (int)uABC[0]/THREADS_NUMBER+1;}

    if (uABC[1]%THREADS_NUMBER==0){
    		dNb = uABC[1]/THREADS_NUMBER;
    }else{	dNb = (int)uABC[1]/THREADS_NUMBER+1;}

    if (uABC[2]%THREADS_NUMBER==0){
    		dNc = uABC[2]/THREADS_NUMBER;
    }else{	dNc = (int)uABC[2]/THREADS_NUMBER+1;}

    int naini=0;
    int nafin=0;
    int nbini=0; 
    int nbfin=0;
    int ncini=0;
    int ncfin=0;

    if (uABC[0]>=uABC[1]&&uABC[0]>=uABC[2]){      //a-axis is the longest side of the box
    	naini = dNa*threadindex; 
    	if (dNa*(threadindex+1)<uABC[0]){
    		nafin = dNa*(threadindex+1);
    	}else{
    		nafin = uABC[0];
    	}
    	nbini = 0; nbfin = uABC[1];
    	ncini = 0; ncfin = uABC[2];
    }else if (uABC[2]>=uABC[0]&&uABC[2]>=uABC[1]){//c-axis is the longest side of the box
    	ncini = dNc*threadindex; 
    	if ( (dNc*(threadindex+1))<uABC[2]){
    		ncfin = dNc*(threadindex+1);
    	}else{
    		ncfin = uABC[2];
    	}
    	naini = 0; nafin = uABC[0];	
    	nbini = 0; nbfin = uABC[1];
    }else if (uABC[1]>=uABC[0]&&uABC[1]>=uABC[2]){//b-axis is the longest side of the box
    	nbini = dNb*threadindex; 
    	if (dNb*(threadindex+1)<uABC[1]){
    		nbfin = dNb*(threadindex+1);
    	}else{
    		nbfin = uABC[1];
    	}
    	naini = 0; nafin = uABC[0];	
    	ncini = 0; ncfin = uABC[2];    	
    }
    // printf("thread[%d]\n",threadindex );
    // printf("dNa = %d, dNb = %d, dNc = %d\n",dNa, dNb, dNc);
    // printf("naini = %d, nafin = %d\n",naini, nafin);
    // printf("nbini = %d, nbfin = %d\n",nbini, nbfin);
    // printf("ncini = %d, ncfin = %d\n",ncini, ncfin);
    // printf("*************************************\n");

	while(true)
	{
		while(ENGINE_MUTEX != DO_IT){ usleep(SleepTime);}
		HacTime = GetACfield(WhichACField);
		if (threadindex==0 && Temperature > 0) GetFluctuations(RNx, RNy, RNz, NOS );
		//StochasticLLG( 	Sx, 	Sy, 	Sz, 
		//StochasticLLG_Heun( 	Sx, 	Sy, 	Sz, 
		//StochasticLLG_RK23( 	Sx, 	Sy, 	Sz,
		StochasticLLG_RK45( 	Sx, 	Sy, 	Sz, 
						tSx,	tSy, 	tSz, 
						Heffx, 	Heffy, 	Heffz, 
						RNx, 	RNy, 	RNz, 
						NOS, 	damping,t_step, 
						Temperature, threadindex,
						naini, 	nafin,
						nbini, 	nbfin,
						ncini, 	ncfin);
		//printf("Tread[%d]\n",threadindex );
		if (threadindex==THREADS_NUMBER-1){ 
			
			//first thread opens the first (in) door in the next (second) thread
			sem_post(sem_in[(threadindex+1)%THREADS_NUMBER]);
			// first (in)door will be open from the last thread (first sem_post)
			sem_wait(sem_in[threadindex]);

			MAX_TORQUE=0;
			for (int i=0;i<THREADS_NUMBER;i++){
				if (Max_torque[i] > MAX_TORQUE) MAX_TORQUE = Max_torque[i];
				Max_torque[i] = 0;
			}
			
			ITERATION++;
			if (ITERATION==Max_Numb_Iteration){
			    Play=0;
			    currentIteration=ITERATION;
			    for (int i=0;i<NOS;i++){
					bSx[i]=Sx[i];
					bSy[i]=Sy[i];
					bSz[i]=Sz[i];
				}
				pthread_mutex_lock(&culc_mutex);
		            ENGINE_MUTEX=WAIT;
		            SleepTime=3000;
				pthread_mutex_unlock(&culc_mutex);
			}
			//normalize all spins every 100 iterations
			// if (ITERATION%100==0) 
			// {
			// 	for (int i=0;i<NOS;i++)
			// 	{
			// 		if (Kind[i]!=0){
			// 		double absS = 1.0f/sqrt(Sx[i]*Sx[i]+Sy[i]*Sy[i]+Sz[i]*Sz[i]);
			// 		Sx[i] = Sx[i] * absS;
			// 		Sy[i] = Sy[i] * absS;
			// 		Sz[i] = Sz[i] * absS;
			// 		}		
			// 	}
			// }

			//save to file if recording is on
			if (Record!=0 && ITERATION%rec_iteration == 0){
				// outputEtotal = GetTotalEnergy( 	bSx, bSy, bSz, 
				// 	NeighborPairs, AIdxBlock, NIdxBlock, NIdxGridA, NIdxGridB, NIdxGridC, SIdx,
				// 	Jij, Bij, Dij, VDMx, VDMy, VDMz, VKu, Ku, Kc, VHf, Hf, Etot0, outputMtotal, NOS );
				outputEtotal = GetTotalEnergyMoment( bSx, bSy, bSz, Heffx, 	Heffy, 	Heffz, Etot0, outputMtotal, NOS);
				// BigDataBank[0][recordsCounter] = (float)ITERATION;
				// BigDataBank[1][recordsCounter] = outputMtotal[0]*iNOS;
				// BigDataBank[2][recordsCounter] = outputMtotal[1]*iNOS;
				// BigDataBank[3][recordsCounter] = outputMtotal[2]*iNOS;
				// BigDataBank[4][recordsCounter] = outputEtotal;
				//metka text LLG
				BigDataBank[0][recordsCounter] = ITERATION*t_step;
				BigDataBank[1][recordsCounter] = bSx[0];
				BigDataBank[2][recordsCounter] = bSy[0];
				BigDataBank[3][recordsCounter] = bSz[1];
				BigDataBank[4][recordsCounter] = bSx[1];
				BigDataBank[5][recordsCounter] = bSy[1];
				BigDataBank[6][recordsCounter] = bSz[1];				recordsCounter++;
				if (recordsCounter==100){
					if (outFile!=NULL){
						for (int i=0; i<recordsCounter; i++){
							//snprintf(BuferString,80,"%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,\n",BigDataBank[0][i],BigDataBank[0][i]*t_step,BigDataBank[1][i],BigDataBank[2][i],BigDataBank[3][i],BigDataBank[4][i]);
							//metka test LLG
							snprintf(BuferString,200,"%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,\n",BigDataBank[0][i],BigDataBank[1][i],BigDataBank[2][i],BigDataBank[3][i],BigDataBank[4][i],BigDataBank[5][i],BigDataBank[6][i]);
							fputs(BuferString,outFile);  							
						}
					}
					printf("%s\n", "Recording to file table.csv is done!");
					recordsCounter=0;
				}
			}

			if (DATA_TRANSFER_MUTEX==WAIT_DATA){
				for (int i=0;i<NOS;i++){
					bSx[i]=Sx[i];
					bSy[i]=Sy[i];
					bSz[i]=Sz[i];
				}
				pthread_mutex_lock(&show_mutex);
					DATA_TRANSFER_MUTEX=TAKE_DATA;
					currentIteration=ITERATION;
				pthread_mutex_unlock(&show_mutex);	
			}


			// now it opens the second (out) door in the next (second) thread
			sem_post(sem_out[(threadindex+1)%THREADS_NUMBER]);
			// second (out)door will be open from the last thread (second sem_post)
			sem_wait(sem_out[threadindex]);

		}else{
			
			//all other calculation threads
			sem_wait(sem_in[threadindex]);
			// first button which open the first door in the next (second) thread
			sem_post(sem_in[(threadindex+1)%THREADS_NUMBER]);
			// second door will be open from the last thread (second button)
			sem_wait(sem_out[threadindex]);
			// second button which open the second door in the next (second) thread
			sem_post(sem_out[(threadindex+1)%THREADS_NUMBER]);
			
		}

}
fclose (outFile);
/* the function must return something - NULL will do */
return NULL;
}