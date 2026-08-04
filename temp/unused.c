/*
 * Archived solver functions with no live callers.
 *
 * This file is intentionally excluded from the build. The E0 functions depend
 * on declarations from Magnoom's unity translation unit and are preserved here
 * only as historical implementations.
 */

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
GetTotalEnergyMomentE0(	magnoom_ctx *ctx, double* sx, double* sy, double* sz, double* Hx, double* Hy, double* Hz, double* Etot, double* Mtot, int N)
{
	double tmp0 = 0;
	Mtot[0] = 0;
	Mtot[1] = 0;
	Mtot[2] = 0;
	double LD = 2*PI/ctx->Dij[0], HD = ctx->Dij[0]*ctx->Dij[0];
	//single spin interactions (or potentila terms): Zeeman and Anizotropy:

	for (int i=0; i<N; i++)
	{

		// Etot[i] = -Hx[i]*sx[i] - Hy[i]*sy[i] - Hz[i]*sz[i];
		Etot[i] = 4*PI*PI*(ctx->BextDCMagnitude/HD)*(1-sx[i]*sin(ctx->BextDCTheta*PI/180)*cos(ctx->BextDCPhi*PI/180)-sy[i]*sin(ctx->BextDCTheta*PI/180)*sin(ctx->BextDCPhi*PI/180)-sz[i]*cos(ctx->BextDCTheta*PI/180))/(LD*LD*LD);
		// metka test stochastic LLG
		double vtemp[3];
		// opposite to the rotation of the vector of external field
		RotateVector(sx[i],sy[i],sz[i],0,0,1,-ctx->BextDCPhi,vtemp); //rotate about y by theta of the external field
		RotateVector(vtemp[0],vtemp[1],vtemp[2],0,1,0,-ctx->BextDCTheta,vtemp); //rotate about z by phi of the external field
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
GetTotalEnergyE0(	magnoom_ctx *ctx, double* sx, double* sy, double* sz,
					int numNeighbors, int* aidxBlock, int* nidxBlock, int* nidxGridA, int* nidxGridB, int* nidxGridC, int* shellIdx,
					float* Jij, float* Bij, float* Dij, float* VDMx, float* VDMy, float* VDMz, float* vku1, float ku1, float* vku2, float ku2, float kc, float* BextDCDirection, float BextDCMagnitude,
					double* Etot,double* Mtot, int N)
{
	double tmp0=0;
	Mtot[0] = 0;
	Mtot[1] = 0;
	Mtot[2] = 0;
	double LD = 2*PI/Dij[0], HD = Dij[0]*Dij[0];


	int Nx = ctx->uABC[0], Ny = ctx->uABC[1], Nz = ctx->uABC[2];
	//single spin interactions (or potentila terms): Zeeman and Anizotropy:
	tmp0=0;
	// if(Nz == 1 && Nx > 1 && Ny > 1){
		for(int i = 0; i < Nx; i++){
			for(int j = 0; j < Ny; j++){
				int n = i + j*Nx, n1 = (i-1+Nx)%Nx+j*Nx, n2 = i + ((j-1+Ny)%Ny)*Nx;
				Etot[n] = 4*PI*PI*(ctx->BextDCMagnitude/HD)*(1-sx[n]*sin(ctx->BextDCTheta*PI/180)*cos(ctx->BextDCPhi*PI/180)-sy[n]*sin(ctx->BextDCTheta*PI/180)*sin(ctx->BextDCPhi*PI/180)-sz[n]*cos(ctx->BextDCTheta*PI/180))/(LD*LD) -\
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
