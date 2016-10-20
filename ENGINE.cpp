void
GetEffectiveField(	double* sx, double* sy, double* sz, 
					int numNeighbors, int* aidxBlock, int* nidxBlock, int* nidxGridA, int* nidxGridB, int* nidxGridC, int* shellIdx,
					float* Jij, float* Bij, float* Dij, float* VDMx, float* VDMy, float* VDMz, float* vku, float ku, float kc, float* VHfield, float Hfield,
					double* heffx, double* heffy, double* heffz, int N)
{
	double tmp0;
	//single spin interactions (or potentila terms): Zeeman and Anizotropy:
	for (int i=0; i<N; i++)
	{
	//H-field (Zeeman energy):
	heffx[i] = Hfield*VHfield[0]+AC_FIELD_ON*HacTime*VHac[0];
	heffy[i] = Hfield*VHfield[1]+AC_FIELD_ON*HacTime*VHac[1];
	heffz[i] = Hfield*VHfield[2]+AC_FIELD_ON*HacTime*VHac[2];
	//uniaxial anisotropy:
	tmp0 = sx[i]*vku[0] + sy[i]*vku[1] + sz[i]*vku[2];
	heffx[i]+= 2 * ku * vku[0] * tmp0;
	heffy[i]+= 2 * ku * vku[1] * tmp0;
	heffz[i]+= 2 * ku * vku[2] * tmp0;
	//cubic anisotropy:
	heffx[i]+= 4 * kc * sx[i]*sx[i]*sx[i];
	heffy[i]+= 4 * kc * sy[i]*sy[i]*sy[i];
	heffz[i]+= 4 * kc * sz[i]*sz[i]*sz[i];	
	}

	// pairwise spin interactions
	int Ip, I, J, K, L;
	int S;
	float Je, Bq, dx, dy, dz, DM;
	int i,j;
	int na1, Na = ABC[0];
	int nb1, Nb = ABC[1];
	int nc1, Nc = ABC[2];
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
		Je= Jij[S];//Comapre to the expression for the total energy
		//Je= Jij[S]*Jexc[ni];
		Bq= Bij[S];//where there is a factor 1/2 which commes from the double summation
		DM= Dij[S];//Here one has not the factor 1/2 because Ei=Heff*Si
		dx= VDMx[ni];
		dy= VDMy[ni];
		dz= VDMz[ni];
		for (int nc=0; nc<Nc; nc++)//nc(a,b)=neghbor in the direction of c(a,b)-vector
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
					nc1 = Na * Nb * ((Nc + nc + L)%Nc);					
					j = I  + AtomsPerBlock * ( na1 + nb1 + nc1 );// index of neighbouring spin "j"
					//Symmetric Heisenberg exchange:
					// heffx[i]+= Je * sx[j] * bc_f;
					// heffy[i]+= Je * sy[j] * bc_f;
					// heffz[i]+= Je * sz[j] * bc_f;
					// //bi-quadratic exchange:
					// tmp0 = (sx[i]*sx[j] + sy[i]*sy[j] + sz[i]*sz[j]) * bc_f;
					// heffx[i]+= 2 * Bq * sx[j] * tmp0;
					// heffy[i]+= 2 * Bq * sy[j] * tmp0;
					// heffz[i]+= 2 * Bq * sz[j] * tmp0;
					// //Dzyaloshinskii-Moriya interaction (antisymmetric exchange):
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
		tmp0 = tmp0 + Etot[i];
	}
	return tmp0;
}

double
GetTotalEnergyFerro(double sx, double sy, double sz, 
					int numNeighbors, int* aidxBlock, int* nidxBlock, int* nidxGridA, int* nidxGridB, int* nidxGridC, int* shellIdx,
					float* Jij, float* Bij, float* Dij, float* VDMx, float* VDMy, float* VDMz, float* vku, float ku, float kc, float* VHfield, float Hfield,
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
	//uniaxial anisotropy:
	tmp0 = sx*vku[0] + sy*vku[1] + sz*vku[2];
	Etot[i]-= ku * tmp0 * tmp0;
	//cubic anisotropy:
	Etot[i]-= kc * (sx*sx*sx*sx + sy*sy*sy*sy + sz*sz*sz*sz);	
	}
	// pairwise spin interactions
	int Ip, I, J, K, L;
	int S;
	double Je, Bq, dx, dy, dz, DM;
	int i,j;
	int na1, Na = ABC[0];
	int nb1, Nb = ABC[1];
	int nc1, Nc = ABC[2];
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
					float* Jij, float* Bij, float* Dij, float* VDMx, float* VDMy, float* VDMz, float* vku, float ku, float kc, float* VHfield, float Hfield,
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
	tmp0 = sx[i]*vku[0] + sy[i]*vku[1] + sz[i]*vku[2];
	Etot[i]-= ku * tmp0 * tmp0;
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
	int na1, Na = ABC[0];
	int nb1, Nb = ABC[1];
	int nc1, Nc = ABC[2];
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
				float temperature)
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

	if (temperature>0) GetFluctuations( rx, ry, rz, nos );
	//electric DC current vector (VCu) and density (Cu)
	Cx = VCu[0] * Cu;
	Cy = VCu[1] * Cu;
	Cz = VCu[2] * Cu;
	GetEffectiveField( 	inx, iny, inz, 
						NeighborPairs, AIdxBlock, NIdxBlock, NIdxGridA, NIdxGridB, NIdxGridC, SIdx,
						Jij, Bij, Dij, VDMx, VDMy, VDMz, VKu, Ku, Kc, VHf, Hf, Heffx, Heffy, Heffz, NOS );
	//prediction step of midpoint solver:
	for (int i=0;i<nos;i++)
	{
		nx = inx[i];	ny = iny[i];	nz = inz[i];
		Hx = Heffx[i];	Hy = Heffy[i];	Hz = Heffz[i];
		Rx = rx[i];		Ry = ry[i];		Rz = rz[i];

		// deterministic terms of Landau–Lifshitz equation:
		Ax = 0.5f * h * ( - Hx - alpha * (ny * Hz - nz * Hy) );
		Ay = 0.5f * h * ( - Hy - alpha * (nz * Hx - nx * Hz) );
		Az = 0.5f * h * ( - Hz - alpha * (nx * Hy - ny * Hx) );

		// Spin-torque term
		Ax = Ax + 0.5f * h * ( - alpha * Cx + (ny * Cz - nz * Cy) ); //pay attention to the signe and factors
		Ay = Ay + 0.5f * h * ( - alpha * Cy + (nz * Cx - nx * Cz) );
		Az = Az + 0.5f * h * ( - alpha * Cz + (nx * Cy - ny * Cx) );

		// stochastic terms of Landau–Lifshitz equation:
		Ax = Ax + 0.5f * rh * D * ( - Rx - alpha * (ny * Rz - nz * Ry) );
		Ay = Ay + 0.5f * rh * D * ( - Ry - alpha * (nz * Rx - nx * Rz) );
		Az = Az + 0.5f * rh * D * ( - Rz - alpha * (nx * Ry - ny * Rx) );

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

	GetEffectiveField( 	tnx, tny, tnz, 
						NeighborPairs, AIdxBlock, NIdxBlock, NIdxGridA, NIdxGridB, NIdxGridC, SIdx,
						Jij, Bij, Dij, VDMx, VDMy, VDMz, VKu, Ku, Kc, VHf, Hf, Heffx, Heffy, Heffz, NOS );

	//final step of midpoint solver:
	for (int i=0;i<nos;i++)
		{
			nx = tnx[i];	ny = tny[i];	nz = tnz[i];	// <-- compare this line to corresponding one in prediction step
			Hx = Heffx[i];	Hy = Heffy[i];	Hz = Heffz[i];	// <-- they are new values for heff
			Rx = rx[i];		Ry = ry[i];		Rz = rz[i];		// <-- they are the same values as in prediction step
			// deterministic terms of Landau–Lifshitz equation:
			Ax = 0.5 * h * ( - Hx - alpha * (ny * Hz - nz * Hy) );
			Ay = 0.5 * h * ( - Hy - alpha * (nz * Hx - nx * Hz) );
			Az = 0.5 * h * ( - Hz - alpha * (nx * Hy - ny * Hx) );
			// Spin-torque term
			Ax = Ax + 0.5f * h * ( - alpha * Cx + (ny * Cz - nz * Cy) ); //pay attention to the signe and factors
			Ay = Ay + 0.5f * h * ( - alpha * Cy + (nz * Cx - nx * Cz) );
			Az = Az + 0.5f * h * ( - alpha * Cz + (nx * Cy - ny * Cx) );
			// stochastic terms of Landau–Lifshitz equation:
			Ax = Ax + 0.5 * rh * D * ( - Rx - alpha * (ny * Rz - nz * Ry) );
			Ay = Ay + 0.5 * rh * D * ( - Ry - alpha * (nz * Rx - nx * Rz) );
			Az = Az + 0.5 * rh * D * ( - Rz - alpha * (nx * Ry - ny * Rx) );

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
		
		}	

	if (ITERATION%100!=0) 
	{
		for (int i=0;i<nos;i++)
		{
		detMi = 1.f/sqrt(inx[i]*inx[i]+iny[i]*iny[i]+inz[i]*inz[i]);
		inx[i] = inx[i] * detMi;
		iny[i] = iny[i] * detMi;
		inz[i] = inz[i] * detMi;		
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


float GetACfield()
{
	return AC_FIELD_ON*Hac*sin(Omega_dc*t_step*ITERATION);
}

/* this function is run by the distinct thread */
void *
CALC_THREAD(void *void_ptr)
{
	char	shortBufer[80];
 	outFile = fopen ("sxsysz.csv","w");
	if (outFile!=NULL) {fputs ("iter,sx,sy,sz,e_tot,\n",outFile);}

	float Etotal;
	double m[3];
	while(true)
	{
		while(FLAG_CALC != DO_IT){ 
			usleep(10);
			//nanosleep (&tw, NULL);
		}
		HacTime = GetACfield();
		if (damping<10){
			StochasticLLG( Sx, Sy, Sz, tSx, tSy, tSz, Heffx, Heffy, Heffz, RNx, RNy, RNz, NOS, damping, t_step, Temperature);
		}else{
			//here should be energy minimization function
			StochasticLLG( Sx, Sy, Sz, tSx, tSy, tSz, Heffx, Heffy, Heffz, RNx, RNy, RNz, NOS, 100, t_step, Temperature);
			//SimpleMinimizer( Sx, Sy, Sz, tSx, tSy, tSz, Heffx, Heffy, Heffz, RNx, RNy, RNz, NOS, damping, t_step, Temperature);
		}

		if (FLAG_SHOW==READY)
		{
			for (int i=0;i<NOS;i++)
			{
				bSx[i]=Sx[i];
				bSy[i]=Sy[i];
				bSz[i]=Sz[i];
			}

			EnterCriticalSection(&show_mutex);
				FLAG_SHOW=TAKE_DATA;
				currentIteration=ITERATION;
			LeaveCriticalSection(&show_mutex);	
		}

		if (Record!=0 && ITERATION%rec_iteration == 0){
			Etotal = GetTotalEnergyMoment(Sx, Sy, Sz, Heffx, Heffy, Heffz, Etot, m, NOS);
			// Etotal = GetTotalEnergy( 	Sx, Sy, Sz, 
			// 				NeighborPairs, AIdxBlock, NIdxBlock, NIdxGridA, NIdxGridB, NIdxGridC, SIdx,
			// 	 			Jij, Bij, Dij, VDMx, VDMy, VDMz, VKu, Ku, Kc, VHf, Hf, Etot, Mtot, NOS );
			  if (outFile!=NULL)
			  {
				snprintf(shortBufer,80,"%d,%2.5f,%2.5f,%2.5f,%2.5f,\n",ITERATION,m[0]/NOS,m[1]/NOS,m[2]/NOS,Etotal);
			    fputs (shortBufer,outFile);  		
			  }
			//printf("%d \t %f \t %f \t %f \n",  ITERATION, Mtot[0],Mtot[1],Mtot[2] );
		}
		ITERATION++;
}
fclose (outFile);
/* the function must return something - NULL will do */
return NULL;
}
