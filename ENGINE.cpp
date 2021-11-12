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
					if(nc==0||nc==Nc-1){
					heffx[i] /= 2;
					heffy[i] /= 2;
					heffz[i] /= 2;
					}

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

//series for function:
//     a*k1 + a^2*k2
//  ------------------ = Sum(a^(n)*coef[n-1],n=1,6)
//	1 + a*q1 + a^2*q2
//
void fun2(double k1, double k2, double q1, double q2, double coef[6]){
	double q12 = q1*q1, c1 = q12*q12-3*q12*q2+q2*q2;
	coef[0] = k1;
	coef[1] = k2-k1*q1;
	coef[2] = -k2*q1+k1*(q12-q2);
	coef[3] = k2*q12-k1*q1*q12-k2*q2+2*k1*q1*q2;
	coef[4] = -k2*(q1*q12-2*q1*q2)+k1*c1;
	coef[5] = k2*c1-k1*q1*(q12*q12-4*q12*q2+3*q2*q2);
}

//series for function:
//     a*k1 + a^2*k2 + a^3*k3 + a^4*k4
//  ------------------------------------ = Sum(a^(n)*coef[n-1],n=1,6)
//	1 + a*q1 + a^2*q2 + a^3*q3 + a^4*q4
//
void fun4(double k1, double k2, double k3, double k4,
 double q1, double q2, double q3, double q4, double* coef){
 	double c1,c2,c3,c4, q12 = q1*q1, q13 = q12*q1;
 	c1 = q12*q12-3*q12*q2+q2*q2+2*q1*q3-q4;
 	c2 = q12*q13-4*q13*q2+3*q1*q2*q2+3*q12*q3-2*q2*q3-2*q1*q4;
 	c3 = q12-q2;
 	c4 = q13-2*q1*q2+q3;
	coef[0] = k1;
	coef[1] = k2-k1*q1;
	coef[2] = k3-k2*q1+k1*c3;
	coef[3] = k4-k3*q1+k2*c3-k1*c4;
	coef[4] = -k4*q1+k3*c3-k2*c4+k1*c1;
	coef[5] = k4*c3-k3*c4+k2*c1-k1*c2;
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

		Etot[i] = 1.0-Hx[i]*sx[i] - Hy[i]*sy[i] - Hz[i]*sz[i];
		// metka test stochastic LLG
		double vtemp[3];
		// opposite to the rotation of the vector of external field
		RotateVector(sx[i],sy[i],sz[i],0,0,1,-VHphi,vtemp); //rotate about y by theta of the external field
		RotateVector(vtemp[0],vtemp[1],vtemp[2],0,1,0,-VHtheta,vtemp); //rotate about z by phi of the external field
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

double
GetTotalEnergyMomentE0(	double* sx, double* sy, double* sz, double* Hx, double* Hy, double* Hz, double* Etot, double* Mtot, int N)
{
	double tmp0 = 0;
	Mtot[0] = 0;
	Mtot[1] = 0;
	Mtot[2] = 0;
	double LD = 2*PI/Dij[0], HD = Dij[0]*Dij[0];
	//single spin interactions (or potentila terms): Zeeman and Anizotropy:

	for (int i=0; i<N; i++)
	{

		// Etot[i] = -Hx[i]*sx[i] - Hy[i]*sy[i] - Hz[i]*sz[i];
		Etot[i] = 4*PI*PI*(Hf/HD)*(1-sx[i]*sin(VHtheta*PI/180)*cos(VHphi*PI/180)-sy[i]*sin(VHtheta*PI/180)*sin(VHphi*PI/180)-sz[i]*cos(VHtheta*PI/180))/(LD*LD*LD);
		// metka test stochastic LLG
		double vtemp[3];
		// opposite to the rotation of the vector of external field
		RotateVector(sx[i],sy[i],sz[i],0,0,1,-VHphi,vtemp); //rotate about y by theta of the external field
		RotateVector(vtemp[0],vtemp[1],vtemp[2],0,1,0,-VHtheta,vtemp); //rotate about z by phi of the external field
		Mtot[0] = Mtot[0] + vtemp[0];
		Mtot[1] = Mtot[1] + vtemp[1];
		Mtot[2] = Mtot[2] + vtemp[2];		
		// Mtot[0] = Mtot[0] + sx[i];
		// Mtot[1] = Mtot[1] + sy[i];
		// Mtot[2] = Mtot[2] + sz[i];
		tmp0 = tmp0 + Etot[i];
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
	int Ip, J, K, L;
	int S;
	double Je, Bq, dx, dy, dz, DM;
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
					int i = Ip + AtomsPerBlock * ( na1 + nb1 + nc1 );// index of spin "i"
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
	for (int i=0;i<NOS;i++)
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

double
GetTotalEnergyE0(	double* sx, double* sy, double* sz, 
					int numNeighbors, int* aidxBlock, int* nidxBlock, int* nidxGridA, int* nidxGridB, int* nidxGridC, int* shellIdx,
					float* Jij, float* Bij, float* Dij, float* VDMx, float* VDMy, float* VDMz, float* vku1, float ku1, float* vku2, float ku2, float kc, float* VHfield, float Hfield,
					double* Etot,double* Mtot, int N)
{
	double tmp0=0;
	Mtot[0] = 0;
	Mtot[1] = 0;
	Mtot[2] = 0;
	double LD = 2*PI/Dij[0], HD = Dij[0]*Dij[0];


	int Nx = uABC[0], Ny = uABC[1], Nz = uABC[2];
	//single spin interactions (or potentila terms): Zeeman and Anizotropy:
	tmp0=0;
	// if(Nz == 1 && Nx > 1 && Ny > 1){
		for(int i = 0; i < Nx; i++){
			for(int j = 0; j < Ny; j++){
				int n = i + j*Nx, n1 = (i-1+Nx)%Nx+j*Nx, n2 = i + ((j-1+Ny)%Ny)*Nx;
				Etot[n] = 4*PI*PI*(Hf/HD)*(1-sx[n]*sin(VHtheta*PI/180)*cos(VHphi*PI/180)-sy[n]*sin(VHtheta*PI/180)*sin(VHphi*PI/180)-sz[n]*cos(VHtheta*PI/180))/(LD*LD) -\
				(sx[n]*sx[n1]+sy[n]*sy[n1]+sz[n]*sz[n1]-1) + (1./LD)*2.*PI*(sy[n]*sz[n1]-sz[n]*sy[n1]) -\
				(sx[n]*sx[n2]+sy[n]*sy[n2]+sz[n]*sz[n2]-1) + (1./LD)*2.*PI*(sz[n]*sx[n2]-sx[n]*sz[n2]);
				tmp0 = tmp0 + Etot[n];
			}
		}
	// }
	// if(Nz > 1 && Nx > 1 && Ny > 1){

	// }

	
	

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
GetFluctuations( float* rx, float* ry, float* rz, int N ){
	for (int i=0; i<N; i++){
		float r[3];
		float A=sqrt(2*abs(log(t_step)));
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
	double Alpha_d = alpha;// / (1.0 + alpha * alpha * Precession);
	double Alpha_p = 1.0f;// / (1.0 + alpha * alpha);

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
				float temperature, float Xi, float Curr_u, 	int thread,					
				int naini, 	int nafin,
				int nbini, 	int nbfin,
				int ncini, 	int ncfin)

{
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
	double Alpha_d = alpha / (1.0 + alpha * alpha * Precession);
	double Alpha_p = 1.0f / (1.0 + alpha * alpha);

	double xi = Xi; 
	double u = Curr_u;
	double c1 = (xi-alpha)*u, c2=(1+xi*alpha)*u/alpha;
	// c2 = (alpha>1e-20)? (1+xi*alpha)*u/alpha : 0;

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

					//Li-Zhang
					int ip = Ip + AtomsPerBlock * ( (na +1)%Na + nb1 + nc1 );
					int im = Ip + AtomsPerBlock * ( (na-1+Na)%Na + nb1 + nc1 );
					int ipp = Ip + AtomsPerBlock * ( (na +2)%Na + nb1 + nc1 );
					int imm = Ip + AtomsPerBlock * ( (na-2+Na)%Na + nb1 + nc1 );
				
					Hx1 = Hx + c1*(inx[ip]-inx[im])/2;
					Hy1 = Hy + c1*(iny[ip]-iny[im])/2;
					Hz1 = Hz + c1*(inz[ip]-inz[im])/2;

					Hx2 = Hx + c2*(inx[ip]-inx[im])/2;
					Hy2 = Hy + c2*(iny[ip]-iny[im])/2;
					Hz2 = Hz + c2*(inz[ip]-inz[im])/2;
					// Hx1 = Hx + c1*(inx[imm]/12 - 2*inx[im]/3 + 2*inx[ip]/3 -inx[ipp]/12);
					// Hy1 = Hy + c1*(iny[imm]/12 - 2*iny[im]/3 + 2*iny[ip]/3 -iny[ipp]/12);
					// Hz1 = Hz + c1*(inz[imm]/12 - 2*inz[im]/3 + 2*inz[ip]/3 -inz[ipp]/12);

					// Hx2 = Hx + c2*(inx[imm]/12 - 2*inx[im]/3 + 2*inx[ip]/3 -inx[ipp]/12);
					// Hy2 = Hy + c2*(iny[imm]/12 - 2*iny[im]/3 + 2*iny[ip]/3 -iny[ipp]/12);
					// Hz2 = Hz + c2*(inz[imm]/12 - 2*inz[im]/3 + 2*inz[ip]/3 -inz[ipp]/12);

					Rx = rx[i];		Ry = ry[i];		Rz = rz[i];

					// deterministic terms of Landau–Lifshitz equation:
					// Fx = Alpha_p * Hx + Alpha_d * (ny * Hz - nz * Hy);
					// Fy = Alpha_p * Hy + Alpha_d * (nz * Hx - nx * Hz);
					// Fz = Alpha_p * Hz + Alpha_d * (nx * Hy - ny * Hx);
					Fx = Alpha_p * Hx1 + Alpha_d * (ny * Hz2 - nz * Hy2);
					Fy = Alpha_p * Hy1 + Alpha_d * (nz * Hx2 - nx * Hz2);
					Fz = Alpha_p * Hz1 + Alpha_d * (nx * Hy2 - ny * Hx2);

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
					int ip = Ip + AtomsPerBlock * ( (na +1)%Na + nb1 + nc1 );
					int im = Ip + AtomsPerBlock * ( (na-1+Na)%Na + nb1 + nc1 );
					int ipp = Ip + AtomsPerBlock * ( (na +2)%Na + nb1 + nc1 );
					int imm = Ip + AtomsPerBlock * ( (na-2+Na)%Na + nb1 + nc1 );
					nx = tnx[i];	ny = tny[i];	nz = tnz[i];	// <-- compare this line to corresponding one in prediction step
					Hx = Heffx[i];	Hy = Heffy[i];	Hz = Heffz[i];	// <-- they are new values for heff
					Hx1 = Hx + c1*(tnx[ip]-tnx[im])/2;
					Hy1 = Hy + c1*(tny[ip]-tny[im])/2;
					Hz1 = Hz + c1*(tnz[ip]-tnz[im])/2;

					Hx2 = Hx + c2*(tnx[ip]-tnx[im])/2;
					Hy2 = Hy + c2*(tny[ip]-tny[im])/2;
					Hz2 = Hz + c2*(tnz[ip]-tnz[im])/2;
					Rx = rx[i];		Ry = ry[i];		Rz = rz[i];		// <-- they are the same values as in prediction step

					// deterministic terms of Landau–Lifshitz equation:
					// Fx = Alpha_p * Hx + Alpha_d * (ny * Hz - nz * Hy);
					// Fy = Alpha_p * Hy + Alpha_d * (nz * Hx - nx * Hz);
					// Fz = Alpha_p * Hz + Alpha_d * (nx * Hy - ny * Hx);
					Fx = Alpha_p * Hx1 + Alpha_d * (ny * Hz2 - nz * Hy2);
					Fy = Alpha_p * Hy1 + Alpha_d * (nz * Hx2 - nx * Hz2);
					Fz = Alpha_p * Hz1 + Alpha_d * (nx * Hy2 - ny * Hx2);

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
				float temperature, float Xi, float Curr_u,	int thread,					
				int naini, 	int nafin,
				int nbini, 	int nbfin,
				int ncini, 	int ncfin)
{
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

	double xi = Xi; 
	double u = Curr_u;
	double c1 = (xi-alpha)*u, c2=(1+xi*alpha)*u/alpha;

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

					int ip = Ip + AtomsPerBlock * ( (na +1)%Na + nb1 + nc1 );
					int im = Ip + AtomsPerBlock * ( (na-1+Na)%Na + nb1 + nc1 );
					int ipp = Ip + AtomsPerBlock * ( (na +2)%Na + nb1 + nc1 );
					int imm = Ip + AtomsPerBlock * ( (na-2+Na)%Na + nb1 + nc1 );
					// Hx1 = Hx + c1*(inx[imm]/12 - 2*inx[im]/3 + 2*inx[ip]/3 -inx[ipp]/12);
					// Hy1 = Hy + c1*(iny[imm]/12 - 2*iny[im]/3 + 2*iny[ip]/3 -iny[ipp]/12);
					// Hz1 = Hz + c1*(inz[imm]/12 - 2*inz[im]/3 + 2*inz[ip]/3 -inz[ipp]/12);

					// Hx2 = Hx + c2*(inx[imm]/12 - 2*inx[im]/3 + 2*inx[ip]/3 -inx[ipp]/12);
					// Hy2 = Hy + c2*(iny[imm]/12 - 2*iny[im]/3 + 2*iny[ip]/3 -iny[ipp]/12);
					// Hz2 = Hz + c2*(inz[imm]/12 - 2*inz[im]/3 + 2*inz[ip]/3 -inz[ipp]/12);
					Hx1 = Hx + c1*(inx[ip]-inx[im])/2;
					Hy1 = Hy + c1*(iny[ip]-iny[im])/2;
					Hz1 = Hz + c1*(inz[ip]-inz[im])/2;

					Hx2 = Hx + c2*(inx[ip]-inx[im])/2;
					Hy2 = Hy + c2*(iny[ip]-iny[im])/2;
					Hz2 = Hz + c2*(inz[ip]-inz[im])/2;

					// deterministic terms of Landau–Lifshitz equation:
					// Fx = Alpha_p * Hx + Alpha_d * (ny * Hz - nz * Hy);
					// Fy = Alpha_p * Hy + Alpha_d * (nz * Hx - nx * Hz);
					// Fz = Alpha_p * Hz + Alpha_d * (nx * Hy - ny * Hx);
					Fx = Alpha_p * Hx1 + Alpha_d * (ny * Hz2 - nz * Hy2);
					Fy = Alpha_p * Hy1 + Alpha_d * (nz * Hx2 - nx * Hz2);
					Fz = Alpha_p * Hz1 + Alpha_d * (nx * Hy2 - ny * Hx2);

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

					int ip = Ip + AtomsPerBlock * ( (na +1)%Na + nb1 + nc1 );
					int im = Ip + AtomsPerBlock * ( (na-1+Na)%Na + nb1 + nc1 );
					int ipp = Ip + AtomsPerBlock * ( (na +2)%Na + nb1 + nc1 );
					int imm = Ip + AtomsPerBlock * ( (na-2+Na)%Na + nb1 + nc1 );
					// Hx1 = Hx + c1*(tnx[imm]/12 - 2*tnx[im]/3 + 2*tnx[ip]/3 -tnx[ipp]/12);
					// Hy1 = Hy + c1*(tny[imm]/12 - 2*tny[im]/3 + 2*tny[ip]/3 -tny[ipp]/12);
					// Hz1 = Hz + c1*(tnz[imm]/12 - 2*tnz[im]/3 + 2*tnz[ip]/3 -tnz[ipp]/12);

					// Hx2 = Hx + c2*(tnx[imm]/12 - 2*tnx[im]/3 + 2*tnx[ip]/3 -tnx[ipp]/12);
					// Hy2 = Hy + c2*(tny[imm]/12 - 2*tny[im]/3 + 2*tny[ip]/3 -tny[ipp]/12);
					// Hz2 = Hz + c2*(tnz[imm]/12 - 2*tnz[im]/3 + 2*tnz[ip]/3 -tnz[ipp]/12);

					Hx1 = Hx + c1*(tnx[ip]-tnx[im])/2;
					Hy1 = Hy + c1*(tny[ip]-tny[im])/2;
					Hz1 = Hz + c1*(tnz[ip]-tnz[im])/2;

					Hx2 = Hx + c2*(tnx[ip]-tnx[im])/2;
					Hy2 = Hy + c2*(tny[ip]-tny[im])/2;
					Hz2 = Hz + c2*(tnz[ip]-tnz[im])/2;

					// deterministic terms of Landau–Lifshitz equation:
					// Fx = Alpha_p * Hx + Alpha_d * (ny * Hz - nz * Hy);
					// Fy = Alpha_p * Hy + Alpha_d * (nz * Hx - nx * Hz);
					// Fz = Alpha_p * Hz + Alpha_d * (nx * Hy - ny * Hx);
					Fx = Alpha_p * Hx1 + Alpha_d * (ny * Hz2 - nz * Hy2);
					Fy = Alpha_p * Hy1 + Alpha_d * (nz * Hx2 - nx * Hz2);
					Fz = Alpha_p * Hz1 + Alpha_d * (nx * Hy2 - ny * Hx2);


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

					int ip = Ip + AtomsPerBlock * ( (na +1)%Na + nb1 + nc1 );
					int im = Ip + AtomsPerBlock * ( (na-1+Na)%Na + nb1 + nc1 );
					int ipp = Ip + AtomsPerBlock * ( (na +2)%Na + nb1 + nc1 );
					int imm = Ip + AtomsPerBlock * ( (na-2+Na)%Na + nb1 + nc1 );
					// Hx1 = Hx + c1*(tnx[imm]/12 - 2*tnx[im]/3 + 2*tnx[ip]/3 -tnx[ipp]/12);
					// Hy1 = Hy + c1*(tny[imm]/12 - 2*tny[im]/3 + 2*tny[ip]/3 -tny[ipp]/12);
					// Hz1 = Hz + c1*(tnz[imm]/12 - 2*tnz[im]/3 + 2*tnz[ip]/3 -tnz[ipp]/12);

					// Hx2 = Hx + c2*(tnx[imm]/12 - 2*tnx[im]/3 + 2*tnx[ip]/3 -tnx[ipp]/12);
					// Hy2 = Hy + c2*(tny[imm]/12 - 2*tny[im]/3 + 2*tny[ip]/3 -tny[ipp]/12);
					// Hz2 = Hz + c2*(tnz[imm]/12 - 2*tnz[im]/3 + 2*tnz[ip]/3 -tnz[ipp]/12);
					Hx1 = Hx + c1*(tnx[ip]-tnx[im])/2;
					Hy1 = Hy + c1*(tny[ip]-tny[im])/2;
					Hz1 = Hz + c1*(tnz[ip]-tnz[im])/2;

					Hx2 = Hx + c2*(tnx[ip]-tnx[im])/2;
					Hy2 = Hy + c2*(tny[ip]-tny[im])/2;
					Hz2 = Hz + c2*(tnz[ip]-tnz[im])/2;

					// deterministic terms of Landau–Lifshitz equation:
					// Fx = Alpha_p * Hx + Alpha_d * (ny * Hz - nz * Hy);
					// Fy = Alpha_p * Hy + Alpha_d * (nz * Hx - nx * Hz);
					// Fz = Alpha_p * Hz + Alpha_d * (nx * Hy - ny * Hx);
					Fx = Alpha_p * Hx1 + Alpha_d * (ny * Hz2 - nz * Hy2);
					Fy = Alpha_p * Hy1 + Alpha_d * (nz * Hx2 - nx * Hz2);
					Fz = Alpha_p * Hz1 + Alpha_d * (nx * Hy2 - ny * Hx2);

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

					int ip = Ip + AtomsPerBlock * ( (na +1)%Na + nb1 + nc1 );
					int im = Ip + AtomsPerBlock * ( (na-1+Na)%Na + nb1 + nc1 );
					int ipp = Ip + AtomsPerBlock * ( (na +2)%Na + nb1 + nc1 );
					int imm = Ip + AtomsPerBlock * ( (na-2+Na)%Na + nb1 + nc1 );
					// Hx1 = Hx + c1*(tnx[imm]/12 - 2*tnx[im]/3 + 2*tnx[ip]/3 -tnx[ipp]/12);
					// Hy1 = Hy + c1*(tny[imm]/12 - 2*tny[im]/3 + 2*tny[ip]/3 -tny[ipp]/12);
					// Hz1 = Hz + c1*(tnz[imm]/12 - 2*tnz[im]/3 + 2*tnz[ip]/3 -tnz[ipp]/12);

					// Hx2 = Hx + c2*(tnx[imm]/12 - 2*tnx[im]/3 + 2*tnx[ip]/3 -tnx[ipp]/12);
					// Hy2 = Hy + c2*(tny[imm]/12 - 2*tny[im]/3 + 2*tny[ip]/3 -tny[ipp]/12);
					// Hz2 = Hz + c2*(tnz[imm]/12 - 2*tnz[im]/3 + 2*tnz[ip]/3 -tnz[ipp]/12);
					Hx1 = Hx + c1*(tnx[ip]-tnx[im])/2;
					Hy1 = Hy + c1*(tny[ip]-tny[im])/2;
					Hz1 = Hz + c1*(tnz[ip]-tnz[im])/2;

					Hx2 = Hx + c2*(tnx[ip]-tnx[im])/2;
					Hy2 = Hy + c2*(tny[ip]-tny[im])/2;
					Hz2 = Hz + c2*(tnz[ip]-tnz[im])/2;
					// deterministic terms of Landau–Lifshitz equation:
					// Fx = Alpha_p * Hx + Alpha_d * (ny * Hz - nz * Hy);
					// Fy = Alpha_p * Hy + Alpha_d * (nz * Hx - nx * Hz);
					// Fz = Alpha_p * Hz + Alpha_d * (nx * Hy - ny * Hx);
					Fx = Alpha_p * Hx1 + Alpha_d * (ny * Hz2 - nz * Hy2);
					Fy = Alpha_p * Hy1 + Alpha_d * (nz * Hx2 - nx * Hz2);
					Fz = Alpha_p * Hz1 + Alpha_d * (nx * Hy2 - ny * Hx2);

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
Relax(	double* inx,		double* iny,		double* inz,		// input vector field
				double* tnx,		double* tny,		double* tnz,		// temporal storage
				double* heffx,	double* heffy,	double* heffz,	// effective field
				float* rx,		float* ry,		float* rz,		// random numbers 
				int nos,		float alpha, 	float h,		// number of spins, damping, time step
				float temperature, 	int thread,					
				int naini, 	int nafin,
				int nbini, 	int nbfin,
				int ncini, 	int ncfin, bool* proj)
{
	
		GetEffectiveField( 	inx, iny, inz, 
					NeighborPairs, AIdxBlock, NIdxBlock, NIdxGridA, NIdxGridB, NIdxGridC, SIdx,
					Jij, Bij, Dij, VDMx, VDMy, VDMz, VKu1, Ku1, VKu2, Ku2, Kc, VHf, Hf, Heffx, Heffy, Heffz, NOS,
					naini, nafin, nbini, nbfin, ncini, ncfin);
	
		

			double Hx, Hy, Hz, temp;// components of the effective field
			int Na = uABC[0];
			int Nb = uABC[1];
			int nb1, nc1;
			int i,s;

			double ALPHA,g1,g2, d1,d2;
			ALPHA = damping;

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
							Hx = Heffx[i];	Hy = Heffy[i];	Hz = Heffz[i];
							//find max torque
							temp = Heffx[i]*inx[i]+Heffy[i]*iny[i]+Heffz[i]*inz[i];
							Hx = Heffx[i]-temp*inx[i];
							Hy = Heffy[i]-temp*iny[i];	
							Hz = Heffz[i]-temp*inz[i];
							temp = sqrt(Hx*Hx + Hy*Hy + Hz*Hz);
							if (temp > Max_torque[thread]) Max_torque[thread] = temp;
							// constant step descent
							proj[i] = (inz[i]>0)? true : false;
							s = (proj[i])? 1 : -1;
							g1 = inx[i]/(1+s*inz[i]);
							g2 = iny[i]/(1+s*inz[i]);
							
							Hx = Heffx[i];	Hy = Heffy[i];	Hz = Heffz[i];
							d1 = (Hx*(iny[i]*iny[i]+s*inz[i]*(1+s*inz[i])) - Hy*inx[i]*iny[i] - Hz*inx[i]*(s+inz[i]));
							d2 = (-Hx*inx[i]*iny[i] + Hy*(inx[i]*inx[i]+s*inz[i]*(1+s*inz[i])) - Hz*iny[i]*(s+inz[i]));

							g1 += ALPHA*d1;
							g2 += ALPHA*d2;

							inx[i] = 2*g1/(1+g1*g1+g2*g2);
							iny[i] = 2*g2/(1+g1*g1+g2*g2);
							inz[i] = s*(1-g1*g1-g2*g2)/(1+g1*g1+g2*g2);

							bSx[i]=inx[i];
							bSy[i]=iny[i];
							bSz[i]=inz[i];					
						}
					}
				}
			}
	
	
}



double GetACfield(int type)
{
	double R=0;
	double temp;
	switch (type){
		case SIN_FIELD:
			R = AC_FIELD_ON*Hac*sin(Omega_dc*t_step*ITERATION);
		break;

		case GAUSSIAN_FIELD:
			temp = (t_step*ITERATION-t_offset)/GPulseWidth;
			R = AC_FIELD_ON*Hac*exp(-0.5*temp*temp );
		break;

		case SINC_FIELD:
			if (t_step*ITERATION<=GPulseWidth){
				if (t_step*ITERATION==t_offset){
					R = AC_FIELD_ON*Hac*1;
				}else{
					R = AC_FIELD_ON*Hac*sin((Omega_dc*(t_step*ITERATION-t_offset)))/(Omega_dc*(t_step*ITERATION-t_offset));
				}				
			}else{ R = 0.0; }
		break;

		case CIRCULAR_FIELD:
			if(AC_FIELD_ON!=0){
				VHac[0] = cos(Omega_dc*t_step*ITERATION);
				VHac[1] = sin(Omega_dc*t_step*ITERATION);
				VHac[2] = 0.0f;
				R = AC_FIELD_ON*Hac;			
			}
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
    // printf("threadindex =%d\n", threadindex );
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
    	// nbini = 1; nbfin = uABC[1];//metka
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

	while(true)
	{
		while(ENGINE_MUTEX != DO_IT){ usleep(SleepTime);}
		HacTime = GetACfield(WhichACField);
		
		switch (WhichIntegrationScheme){
			case HEUN: 
				// SimpleMinimizer(Sx,Sy,Sz,tSx,tSy,tSz,Heffx,Heffy,Heffz,RNx,RNy,RNz,NOS,damping,t_step,Temperature, threadindex,naini,nafin,nbini,nbfin,ncini,ncfin);
				StochasticLLG_Heun(Sx,Sy,Sz,tSx,tSy,tSz,Heffx,Heffy,Heffz,RNx,RNy,RNz,NOS,damping,t_step,Temperature, threadindex,naini,nafin,nbini,nbfin,ncini,ncfin);
			break;

			case SIB: 
				StochasticLLG(Sx,Sy,Sz,tSx,tSy,tSz,Heffx,Heffy,Heffz,RNx,RNy,RNz,NOS,damping,t_step,Temperature, threadindex,naini,nafin,nbini,nbfin,ncini,ncfin);
			break;

			case RK23: 
				StochasticLLG_RK23(Sx,Sy,Sz,tSx,tSy,tSz,Heffx,Heffy,Heffz,RNx,RNy,RNz,NOS,damping,t_step,Temperature,Xi,Curr_u, threadindex,naini,nafin,nbini,nbfin,ncini,ncfin);
			break;

			case RK45: 
				StochasticLLG_RK45(Sx,Sy,Sz,tSx,tSy,tSz,Heffx,Heffy,Heffz,RNx,RNy,RNz,NOS,damping,t_step,Temperature,Xi,Curr_u, threadindex,naini,nafin,nbini,nbfin,ncini,ncfin);
			break;

			case RELAX: 	
				Relax(Sx,Sy,Sz,tSx,tSy,tSz,Heffx,Heffy,Heffz,RNx,RNy,RNz,NOS,damping,t_step,Temperature, threadindex,naini,nafin,nbini,nbfin,ncini,ncfin, Proj);
			break;
		}
		if (threadindex==0 && Temperature > 0) GetFluctuations(RNx, RNy, RNz, NOS );

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
			if (WhichIntegrationScheme != RELAX && ITERATION%100==0) 
			{
				// printf("ich bin hier!\n");
				for (int i=0;i<NOS;i++)
				{
					if (Kind[i]!=0){
					double absS = 1.0f/sqrt(Sx[i]*Sx[i]+Sy[i]*Sy[i]+Sz[i]*Sz[i]);
					Sx[i] = Sx[i] * absS;
					Sy[i] = Sy[i] * absS;
					Sz[i] = Sz[i] * absS;
					}		
				}
			}



			//save to file if recording is on
			if (Record!=0 && ITERATION%rec_iteration == 0){

				outputEtotal = GetTotalEnergyMoment( bSx, bSy, bSz, Heffx, 	Heffy, 	Heffz, Etot0, outputMtotal, NOS);
				BigDataBank[0][recordsCounter] = (float)ITERATION;
				BigDataBank[1][recordsCounter] = outputMtotal[0]*iNOS;
				BigDataBank[2][recordsCounter] = outputMtotal[1]*iNOS;
				BigDataBank[3][recordsCounter] = outputMtotal[2]*iNOS;
				BigDataBank[4][recordsCounter] = outputEtotal;
				//metka test LLG
				// BigDataBank[0][recordsCounter] = ITERATION*t_step;
				// BigDataBank[1][recordsCounter] = bSx[0];
				// BigDataBank[2][recordsCounter] = bSy[0];
				// BigDataBank[3][recordsCounter] = bSz[1];
				// BigDataBank[4][recordsCounter] = bSx[1];
				// BigDataBank[5][recordsCounter] = bSy[1];
				// BigDataBank[6][recordsCounter] = bSz[1];				
				recordsCounter++;
				if (recordsCounter==100){
					if (outFile!=NULL){
						for (int i=0; i<recordsCounter; i++){
							snprintf(BuferString,80,"%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,\n",BigDataBank[0][i],BigDataBank[0][i]*t_step,BigDataBank[1][i],BigDataBank[2][i],BigDataBank[3][i],BigDataBank[4][i]);
							//metka test LLG
							// snprintf(BuferString,200,"%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,\n",BigDataBank[0][i],BigDataBank[1][i],BigDataBank[2][i],BigDataBank[3][i],BigDataBank[4][i],BigDataBank[5][i],BigDataBank[6][i]);
							fputs(BuferString,outFile);  							
						}
					}
					printf("%s\n", "Recording to file table.csv is done!");
					recordsCounter=0;
				}
			}

			//save mode snapshot
			if (AC_MODE_REC*AC_FIELD_ON!=0){
				float phase = Omega_dc*t_step*ITERATION*iTPI;
				phase = phase - floor(phase);
				int Im=-1;
				float tolerance=0.5*t_step/Period_dc;

				for (int i=0; i<=num_images; i++){
					if ( ABS(phase - i/(float)num_images)<tolerance){
						Im = i%num_images;
					}
				}
				if (Im!=-1){
						current_rec_num_mode++;
						for (int i=0; i<NOS; i++){
							// for delta m:
							dImage_x[Im][i]=dImage_x[Im][i]+(Sx[i]-t3Sx[i]);
							dImage_y[Im][i]=dImage_y[Im][i]+(Sy[i]-t3Sy[i]);
							dImage_z[Im][i]=dImage_z[Im][i]+(Sz[i]-t3Sz[i]);	
							// for m:
							Image_x[Im][i] = Image_x[Im][i] + Sx[i];
							Image_y[Im][i] = Image_y[Im][i] + Sy[i];
							Image_z[Im][i] = Image_z[Im][i] + Sz[i];						
						}
						printf("tolerance %1.8f \n", tolerance);
						// printf("ABS(phase - i/(float)num_images)= %1.8f \n", ABS(phase - Im/(float)num_images));
						printf("Phase %1.8f (%d) (%1.8f)\n", phase,Im,Im/(float)num_images);
						printf("Image %d\n", current_rec_num_mode);
						if (current_rec_num_mode==num_images*rec_num_mode){	
							Play=0;
							AC_MODE_REC=0;
							AC_FIELD_ON=0;
							current_rec_num_mode=0;
							for (int j=0; j<num_images; j++){
								for (int i=0; i<NOS; i++){
								// for dm:
									dImage_x[j][i]=dImage_x[j][i]/rec_num_mode;
									dImage_y[j][i]=dImage_y[j][i]/rec_num_mode;
									dImage_z[j][i]=dImage_z[j][i]/rec_num_mode;
								// for m:
									Image_x[j][i]=Image_x[j][i]/rec_num_mode;
									Image_y[j][i]=Image_y[j][i]/rec_num_mode;
									Image_z[j][i]=Image_z[j][i]/rec_num_mode;
									float Norm=sqrt(Image_x[j][i]*Image_x[j][i]+Image_y[j][i]*Image_y[j][i]+Image_z[j][i]*Image_z[j][i]);					
									Image_x[j][i]=Image_x[j][i]/Norm;
									Image_y[j][i]=Image_y[j][i]/Norm;
									Image_z[j][i]=Image_z[j][i]/Norm;

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
									snprintf(vtk_filename,64,"phase%d.vtk",j);
									Save_VTK_6(Image_x[j], Image_y[j], Image_z[j], dImage_x[j], dImage_y[j], dImage_z[j],0, vtk_filename);
								/*dTheta dPhi*/
								for (int i=0; i<NOS; i++){
								// get theta and phi for the equilibrium state:
									float T=acos(t3Sz[i]);
									float F=atan2(t3Sy[i],t3Sx[i]);
								// get spin i
									float s[3], r[3];
									mat4x4 My, Mz, M;
									s[0]=(float)Image_x[j][i];
									s[1]=(float)Image_y[j][i];
									s[2]=(float)Image_z[j][i];
									mat4x4_identity(My);
									mat4x4_rotate_Y(My, My, T);
									mat4x4_identity(Mz);
									mat4x4_rotate_Z(Mz, Mz, F);
									// mat4x4_mul(M, Mz, My);
			
									mat4x4_mul_vec3(r, My, s);
									mat4x4_mul_vec3(s, Mz, r);
									dImage_x[j][i]=s[0];
									dImage_y[j][i]=s[1];
									dImage_z[j][i]=s[2];	
								}
									snprintf(vtk_filename,64,"dTdF%d.vtk",j);
									Save_VTK(dImage_x[j], dImage_y[j], dImage_z[j],0, vtk_filename);
							}
					}
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
