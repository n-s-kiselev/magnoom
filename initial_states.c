void CreatSkyrmion(magnoom_ctx *ctx, float * px, float * py, float * pz, double * s, float Sk_R, float tx, float ty)
{
	float T = 0.f;
	float F = 0.f;
	float r = 0.f;
	float rx=0.0f;
	float ry=0.0f;

	for (int n=0; n<ctx->NOS; n++)
	{
		rx=px[n]-tx;
		ry=py[n]-ty;
		r = sqrt(rx*rx+ry*ry);
		// if (r<Sk_R)//metka
		{
			T= PI*exp(-0.5*r/Sk_R);//<-- defines skyrmion profile you may put periodical function to get target like skyrmions
			F= atan2(ry,rx)+PI*0.5;//<--chiral (bloch) skyrmion |Q|=1
			VEC_X(ctx->S,n) = sin(T)*cos(F)*ctx->Kind[n];
			VEC_Y(ctx->S,n) = sin(T)*sin(F)*ctx->Kind[n];
			VEC_Z(ctx->S,n) = cos(T)*ctx->Kind[n];
			//metka test  color code:	
			// float inv_abs_S=1.0f/sqrt(VEC_X(ctx->S,n)*VEC_X(ctx->S,n)+VEC_Y(ctx->S,n)*VEC_Y(ctx->S,n)+VEC_Z(ctx->S,n)*VEC_Z(ctx->S,n));
			// VEC_X(ctx->S,n) = VEC_X(ctx->S,n)*inv_abs_S/(0.1*r+1.0);
			// VEC_Y(ctx->S,n) = VEC_Y(ctx->S,n)*inv_abs_S/(0.1*r+1.0);
			// VEC_Z(ctx->S,n) = VEC_Z(ctx->S,n)*inv_abs_S/(0.1*r+1.0);
		}
	}
}

void GetSkyrmion(magnoom_ctx *ctx, float * px, float * py, float * pz, double * s, float Sk_R, float tx, float ty)
{
	float T = 0.f;
	float F = 0.f;
	float r = 0.f;
	float rx=0.0f;
	float ry=0.0f;



	for (int n=0; n<ctx->NOS; n++)
	{
		rx=px[n]-tx;
		ry=py[n]-ty;
		r = sqrt(rx*rx+ry*ry);
		if (r<Sk_R){
			T = PI*(1-r/Sk_R);
			F = atan2(ry,rx)+PI*0.5;
			VEC_X(s,n) = sin(T)*cos(F)*ctx->Kind[n];
			VEC_Y(s,n) = sin(T)*sin(F)*ctx->Kind[n];
			VEC_Z(s,n) = cos(T)*ctx->Kind[n];
		}
	}
}

void GetAntiskyrmion(magnoom_ctx *ctx, float * px, float * py, float * pz, double * s, float Sk_R, float tx, float ty, bool orient)
{
	float T = 0.f;
	float F = 0.f;
	float r = 0.f;
	float rx=0.0f;
	float ry=0.0f;



	for (int n=0; n<ctx->NOS; n++)
	{	
		
		if(orient){
			rx=(px[n]-tx);
			ry=(py[n]-ty)*2;
			F = -atan2(ry,rx)-PI*0.5;	
		}else{
			rx=(px[n]-tx)*2;
			ry=(py[n]-ty);
			F = -atan2(rx,-ry)-PI;
		}
		
		r = sqrt(rx*rx+ry*ry);
		if (r<Sk_R){
			T = PI*(1-r/Sk_R);
			
			VEC_X(s,n) = sin(T)*cos(F)*ctx->Kind[n];
			VEC_Y(s,n) = sin(T)*sin(F)*ctx->Kind[n];
			VEC_Z(s,n) = cos(T)*ctx->Kind[n];
		}
	}
}

void TiltSpinsToX(magnoom_ctx *ctx, float * px, float * py, float * pz, double * s, float Sk_R, float tx, float ty)
{
	float T = 0.f;
	float F = 0.f;
	float r = 0.f;
	float rx=0.0f;
	float ry=0.0f;

	for (int n=0; n<ctx->NOS; n++)
	{
		rx=px[n]-tx;
		ry=py[n]-ty;
		r = sqrt(rx*rx+ry*ry);
		if (r<Sk_R)
		{
			T= PI/2;
			F= 0;
			VEC_X(ctx->S,n) = sin(T)*cos(F)*ctx->Kind[n];
			VEC_Y(ctx->S,n) = sin(T)*sin(F)*ctx->Kind[n];
			VEC_Z(ctx->S,n) = cos(T)*ctx->Kind[n];
		}
	}
}

void CreatSkyrmionSoliton(magnoom_ctx *ctx, float * px, float * py, float * pz, double * s, float Sk_R, float tx, float ty)
{
	float T = 0.f;
	float F = 0.f;
	float A = 0.f;
	float a = ctx->chDir[0];//sqrt(alpha)
	float v = ctx->chDir[1];
	float r = 0.f;
	float rx= 0.f;
	float ry= 0.f;
	float rz= 0.f;
	for (int n=0; n<ctx->NOS; n++){
		rx=px[n];
		ry=py[n];
		rz=pz[n];

		F = 2*a*tanh(a*v*rz*ctx->chDir[2])+v*rz*ctx->chDir[2];
		A = 10*2*a/(fabs(v)*cosh(a*v*rz*ctx->chDir[2]));
		rx= px[n]-A*cos(F);
		ry= py[n]-A*sin(F);
		r = sqrt(rx*rx+ry*ry);
		if (r<Sk_R){
			T = PI*exp(-2*r/Sk_R);
			F = atan2(ry,rx)+PI*0.5;
			VEC_X(ctx->S,n) = sin(T)*cos(F)*ctx->Kind[n];
			VEC_Y(ctx->S,n) = sin(T)*sin(F)*ctx->Kind[n];
			VEC_Z(ctx->S,n) = cos(T)*ctx->Kind[n];			
		}else{
			VEC_X(ctx->S,n) = 0.0f;
			VEC_Y(ctx->S,n) = 0.0f;
			VEC_Z(ctx->S,n) = 1.0f*ctx->Kind[n];
		}
	}
}


void CreatGlobule(magnoom_ctx *ctx, float * px, float * py, float * pz, double * s, float Sk_R, float tx, float ty, float tz)
{
	float T = 0.f;
	float F = 0.f;
	float rx=0.0f;
	float ry=0.0f;
	float rz=0.0f;
	float r1,r2;

	for (int n=0; n<ctx->NOS; n++)
	{
		rx=px[n]-tx;
		ry=py[n]-ty;
		rz=pz[n]-tz;
		r1 = sqrt(rx*rx+ry*ry);
		r2 = sqrt(rx*rx+ry*ry+rz*rz);
		if (r2<ctx->chSize*0.5)
		{
			T=PI*exp(-2*r1/ctx->chSize);//<-- defines skyrmion profile you may put periodical function to get target like skyrmions
			F= atan2(ry,rx)+PI*0.5;//<--chiral (bloch) skyrmion |Q|=1
			VEC_X(ctx->S,n) = sin(T)*cos(F)*ctx->Kind[n];
			VEC_Y(ctx->S,n) = sin(T)*sin(F)*ctx->Kind[n];
			VEC_Z(ctx->S,n) = cos(T)*ctx->Kind[n];			
		}
	}
}

void Rx(float T,double* sx, double* sy, double* sz)
{
	double sxp,syp,szp;
	sxp = *sx;
	syp = *sy * cos(T) - *sz * sin(T);
	szp = *sy * sin(T) + *sz * cos(T);
	*sx = sxp;
	*sy = syp;
	*sz = szp;
}

void Ry(float T, double* sx, double* sy, double* sz)
{
	double sxp,syp,szp;
	sxp = *sx * cos(T) + *sz * sin(T);
	syp = *sy;
	szp =-*sx * sin(T) + *sz * cos(T);
	*sx = sxp;
	*sy = syp;
	*sz = szp;
}

void Rz(float T, double* sx, double* sy, double* sz)
{
	double sxp,syp,szp;
	sxp = *sx * cos(T) - *sy * sin(T);
	syp = *sx * sin(T) - *sy * cos(T);
	szp = *sz;
	*sx = sxp;
	*sy = syp;
	*sz = szp;
}

void CreatBobber(magnoom_ctx *ctx, float * px, float * py, float * pz, double * s, float Bob_R ,float tx, float ty, int Loc_B)
{
	float top = 0.5*ctx->Box[2][2];
	float R;
	float T = 0.f;
	float F = 0.f;
	float r = 0.f;
	//float Bob_R = chSize;
	float Bob_D = 2*Bob_R;
	float rx=0.0f;
	float ry=0.0f;
	if(ctx->chSize>0)
	{
		if (Bob_D>top) { Bob_D = 0.75*top;}

		for (int n=0; n<ctx->NOS; n++)
		{		
			rx=px[n]-tx;
			ry=py[n]-ty;			
			R = Bob_R*(1-(top-Loc_B*pz[n])/(Bob_D));//If the bobber set in positive Loc_B=1 side else bobber set in negative Loc_B=-1 
			r = sqrt(rx*rx+ry*ry);
			if ( (r<R) && (fabs(pz[n])>(top-Bob_D)) )
			{
			T=PI*(R-r)/R;//<-- defines skyrmion profile you may put periodical function to get target like skyrmions
			F= 1*(-atan2(rx,ry));//<--chiral (bloch) skyrmion |Q|=1
			VEC_X(ctx->S,n) = -sin(T)*cos(F)*ctx->Kind[n];
			VEC_Y(ctx->S,n) = -sin(T)*sin(F)*ctx->Kind[n];
			VEC_Z(ctx->S,n) = cos(T)*ctx->Kind[n];	
			}
		}	
	}
}

void CreatHorisontalBobber(magnoom_ctx *ctx, float * px, float * py, float * pz, double * s, float Bob_R){
	for(int i = 0; i < ctx->uABC[0]; i++){
		for(int j = 0; j < ctx->uABC[1]; j++){
			for(int k = 0; k < ctx->uABC[2]; k++){
				int n = i + j*ctx->uABC[0] + k*ctx->uABC[0]*ctx->uABC[1];
				VEC_X(ctx->S,n) = 0;
				VEC_Y(ctx->S,n) = 0;
				VEC_Z(ctx->S,n) = 1;	
				if(px[n]<0 && (py[n]*py[n]+pz[n]*pz[n]) < Bob_R*Bob_R){
					float r = sqrt(py[n]*py[n]+pz[n]*pz[n]);
					float T = PI*(Bob_R-r)/Bob_R;
					float F = atan2(py[n],pz[n]) + 0.5*PI;
					VEC_X(ctx->S,n) = sin(T)*cos(F);
					VEC_Y(ctx->S,n) = sin(T)*sin(F);
					VEC_Z(ctx->S,n) = cos(T);	
				}
			}
		}
	}
}

void InitSpinComponents(magnoom_ctx *ctx, float * px, float * py, float * pz, double * s, int id)
{
	float T = 0.f;
	float F = 0.f;
	// T, F are theta and phi for spherical spase coordinates 
	float r = 0.f;

	switch (id)
	{

	case 0:; //random spins
		double rnd[3];
		for (int n=0; n<ctx->NOS; n++)
		// for (int n=0; n<ctx->uABC[0]; n++)//metka
		{	
			rnd[0] = 2.0 * (0.5 - rand() / (double)RAND_MAX);
			rnd[1] = 2.0 * (0.5 - rand() / (double)RAND_MAX);
			rnd[2] = 2.0 * (0.5 - rand() / (double)RAND_MAX);
			(void)Unit(rnd,rnd);
			VEC_X(ctx->S,n) = rnd[0]*ctx->Kind[n];
			VEC_Y(ctx->S,n) = rnd[1]*ctx->Kind[n];
			VEC_Z(ctx->S,n) = rnd[2]*ctx->Kind[n];	
		}	
		//test
			// VEC_X(ctx->S,0) = 1.0;
			// VEC_Y(ctx->S,0) = 0.0;
			// VEC_Z(ctx->S,0) = 0.0;
			// VEC_X(ctx->S,1) = 0.0;
			// VEC_Y(ctx->S,1) = 0.0;
			// VEC_Z(ctx->S,1) = 1.0;	
		// Stefan		
		// if(chSize>0)
		// {
		// 	if( rand() / (double)RAND_MAX>0.5){
		// 		TiltSpinsToX(ctx, px, py, pz, s, chSize, 0.5*ctx->uABC[0], 0);				
		// 	}

		// 	if( 2.0 * (0.5 - rand() / (double)RAND_MAX)>0){
		// 		TiltSpinsToX(ctx, px, py, pz, s, chSize,-0.5*ctx->uABC[0], 0);				
		// 	}

		// 	// if( 2.0 * (0.5 - rand() / (double)RAND_MAX)>0){
		// 	// 	TiltSpinsToX(ctx, px, py, pz, s, chSize, 0.5*ctx->uABC[0], -0.5*ctx->uABC[1]);				
		// 	// }

		// 	// if( 2.0 * (0.5 - rand() / (double)RAND_MAX)>0){
		// 	// 	TiltSpinsToX(ctx, px, py, pz, s, chSize, 0.5*ctx->uABC[0],  0.5*ctx->uABC[1]);				
		// 	// }
		// }
	break;

	case 1: //homogeneous
		for (int n=0; n<ctx->NOS; n++)
		{	
			VEC_X(ctx->S,n) = ctx->chDir[0]*ctx->Kind[n];
			VEC_Y(ctx->S,n) = ctx->chDir[1]*ctx->Kind[n];
			VEC_Z(ctx->S,n) = ctx->chDir[2]*ctx->Kind[n];	

			//tilted FM
			// VEC_X(ctx->S,n) = sin(BextDCTheta*PI/180)*cos(BextDCPhi*PI/180);
			// VEC_Y(ctx->S,n) = sin(BextDCTheta*PI/180)*sin(BextDCPhi*PI/180);
			// VEC_Z(ctx->S,n) = cos(BextDCTheta*PI/180);

			//cone
			// VEC_X(ctx->S,n) = sin(acos(BextDCMagnitude/(ctx->Dij[0]*ctx->Dij[0])))*cos(pz[n]*2*PI/64);
			// VEC_Y(ctx->S,n) = sin(acos(BextDCMagnitude/(ctx->Dij[0]*ctx->Dij[0])))*sin(pz[n]*2*PI/64);
			// VEC_Z(ctx->S,n) = BextDCMagnitude/(ctx->Dij[0]*ctx->Dij[0]);
			// int kz = n%(ctx->uABC[0]*ctx->uABC[1]);
			// VEC_X(ctx->S,n) = sin(acos(0.7975))*cos(pz[n]*2*PI/128);
			// VEC_Y(ctx->S,n) = sin(acos(0.7975))*sin(pz[n]*2*PI/128);
			// VEC_Z(ctx->S,n) = 0.7975;	
		}

	break;

	case 2: // skyrmion Q=1
		if(ctx->chSize>0)
		{
			CreatSkyrmion(ctx, px, py, pz, s, ctx->chSize,0,0);	
		}
	break;

	case 3: // skyrmion Q=2
		if(ctx->chSize>0)
		{
			// for (int n=0; n<ctx->NOS; n++)
			// {	
			// r = sqrt(px[n]*px[n]+py[n]*py[n]+0*pz[n]*pz[n]);
			// //r = sqrt(px[n]*px[n]+py[n]*py[n]);
			// if (r<chSize){
			// T=PI*exp(-2*r/chSize);//<-- defines skyrmion profile you may put periodical function to get target like skyrmions
			// }
			// else
			// {
			// 	T=0.f;
			// }
			// //try to modify phi angle as shown below:
			// //F= 1*(-atan2(px[n],py[n]));//<--chiral (bloch) skyrmion |Q|=1
			// //F= 1*(-atan2(px[n],py[n]))+PI*0.5;//<--chiral (neel) skyrmion |Q|=1
			// //F=-1*(-atan2(px[n],py[n]))+PI*0.5;//<-- achiral skyrmion |Q|=1
			// F=-2*(-atan2(px[n],py[n]))+PI*0.5;//<-- achiral skyrmion |Q|=2
			// //F=-3*(-atan2(px[n],py[n]))+PI*0.5;//<-- achiral skyrmion |Q|=3
			// VEC_X(ctx->S,n) = sin(T)*cos(F)*ctx->Kind[n];
			// VEC_Y(ctx->S,n) = sin(T)*sin(F)*ctx->Kind[n];
			// VEC_Z(ctx->S,n) = cos(T)*ctx->Kind[n];	
			// }
			//metka eto antiskyrmion
				float T = 0.f;
				float F = 0.f;
				float r = 0.f;
				float rx=0.0f;
				float ry=0.0f;

				for (int n=0; n<ctx->NOS; n++)
				{
					rx=px[n];
					ry=py[n];
					r = sqrt(rx*rx+ry*ry);
					if (r<ctx->chSize)
					{
						T= PI*exp(-2*r/ctx->chSize);//<-- defines skyrmion profile you may put periodical function to get target like skyrmions
						F= atan2(ry,rx)+PI*0.5;//<--chiral (bloch) skyrmion |Q|=1
						VEC_Y(ctx->S,n) = sin(T)*cos(F)*ctx->Kind[n];
						VEC_X(ctx->S,n) = sin(T)*sin(F)*ctx->Kind[n];
						VEC_Z(ctx->S,n) = cos(T)*ctx->Kind[n];
					}
				}
		}
	break;

	case 4: // skyrmion Q=3
		if(ctx->chSize>0)
		{
			for (int n=0; n<ctx->NOS; n++)
			{	
			r = sqrt(px[n]*px[n]+py[n]*py[n]+0*pz[n]*pz[n]);
			//r = sqrt(px[n]*px[n]+py[n]*py[n]);
			if (r<ctx->chSize){
			T=PI*exp(-2*r/ctx->chSize);//<-- defines skyrmion profile you may put periodical function to get target like skyrmions
			}
			else
			{
				T=0.f;
			}
			//try to modify phi angle as shown below:
			//F= 1*(-atan2(px[n],py[n]));//<--chiral (bloch) skyrmion |Q|=1
			//F= 1*(-atan2(px[n],py[n]))+PI*0.5;//<--chiral (neel) skyrmion |Q|=1
			//F=-1*(-atan2(px[n],py[n]))+PI*0.5;//<-- achiral skyrmion |Q|=1
			//F=-2*(-atan2(px[n],py[n]))+PI*0.5;//<-- achiral skyrmion |Q|=2
			F=-7*(-atan2(px[n],py[n]))+PI*0.5;//<-- achiral skyrmion |Q|=3
			VEC_X(ctx->S,n) = sin(T)*cos(F)*ctx->Kind[n];
			VEC_Y(ctx->S,n) = sin(T)*sin(F)*ctx->Kind[n];
			VEC_Z(ctx->S,n) = cos(T)*ctx->Kind[n];	
			}
		}	
	break;

	case 5: // bobber_top
		if(ctx->chSize>0)
		
		{
			CreatBobber(ctx, px, py, pz, s, ctx->chSize, 0, 0, 1);	
		}
	break;

	case 6: // bobber_bottom
		if(ctx->chSize>0)
		
		{
			// CreatBobber(ctx, px, py, pz, s, chSize, 0, 0, -1);	
			CreatHorisontalBobber(ctx, px, py, pz, s, ctx->chSize);			
		}
	break;

	case 7: // bobber_lattice
		if(ctx->chSize>0 && !(ctx->chDir[0]==0&&ctx->chDir[1]==0) )
		{

			float PerBobLat=ctx->chSize;
			float aBobLat=2*PerBobLat/(sqrt(3)); 
			float Bob_Radius=0.8*0.5*aBobLat; 
			float TranLatX,TranLatX1;
			float TranLatY,TranLatY1;
			float CoChDir=ctx->chDir[0]/(sqrt(ctx->chDir[0]*ctx->chDir[0]+ctx->chDir[1]*ctx->chDir[1]));//cos of rotation angle 
			float SiChDir=ctx->chDir[1]/(sqrt(ctx->chDir[0]*ctx->chDir[0]+ctx->chDir[1]*ctx->chDir[1]));//sin of rotation angle
			float RotTranLatX;
			float RotTranLatY;
		
			for (int n=0; n<ctx->NOS; n++)
			{	
				VEC_X(ctx->S,n) = 0.f*ctx->Kind[n];
				VEC_Y(ctx->S,n) = 0.f*ctx->Kind[n];
				VEC_Z(ctx->S,n) = 1.f*ctx->Kind[n];	
			}
			
			int Ntr=10;
			for (int kx=-Ntr; kx<=Ntr; kx++)
			{
				for (int ky=-Ntr; ky<=Ntr; ky++)
				{
					TranLatX=aBobLat*kx+aBobLat*0.5*ky;//Tranlation of Bobber Lattice (along X) without rotate
					TranLatY=PerBobLat*ky;//Tranlation of Bobber Lattice (along Y) without rotate

					RotTranLatX=CoChDir*TranLatX-SiChDir*TranLatY;//Tranlation of Bobber Lattice (along Y) with rotate in XY
					RotTranLatY=SiChDir*TranLatX+CoChDir*TranLatY;//Tranlation of  Lattice (along Y) with rotate in XY
					CreatBobber(ctx, px, py, pz, s, Bob_Radius, RotTranLatX, RotTranLatY, 1);

					TranLatX1=TranLatX+aBobLat/2;//
					TranLatY1=TranLatY+PerBobLat/3;//

					RotTranLatX=CoChDir*TranLatX1-SiChDir*TranLatY1;//Tranlation of Bobber Lattice (along Y) with rotate in XY
					RotTranLatY=SiChDir*TranLatX1+CoChDir*TranLatY1;//Tranlation of  Lattice (along Y) with rotate in XY
					CreatBobber(ctx, px, py, pz, s, Bob_Radius, RotTranLatX, RotTranLatY, -1);
				}
			}	
		}
	break;

	case 8: // bobber_lattice_top
		if(ctx->chSize>0 && !(ctx->chDir[0]==0&&ctx->chDir[1]==0) )
		{

			float PerBobLat=ctx->chSize;
			float aBobLat=2*PerBobLat/(sqrt(3)); 
			float Bob_Radius=0.8*0.5*aBobLat; 
			float TranLatX;
			float TranLatY;
			float CoChDir=ctx->chDir[0]/(sqrt(ctx->chDir[0]*ctx->chDir[0]+ctx->chDir[1]*ctx->chDir[1]));//cos of rotation angle 
			float SiChDir=ctx->chDir[1]/(sqrt(ctx->chDir[0]*ctx->chDir[0]+ctx->chDir[1]*ctx->chDir[1]));//sin of rotation angle
			float RotTranLatX;
			float RotTranLatY;
		
			for (int n=0; n<ctx->NOS; n++)
			{	
				if (pz[n]>0)
				{
				VEC_X(ctx->S,n) = 0.f*ctx->Kind[n];
				VEC_Y(ctx->S,n) = 0.f*ctx->Kind[n];
				VEC_Z(ctx->S,n) = 1.f*ctx->Kind[n];
				}	
			}
			
			int Ntr=10;
			for (int kx=-Ntr; kx<=Ntr; kx++)
			{
				for (int ky=-Ntr; ky<=Ntr; ky++)
				{
					TranLatX=aBobLat*kx+aBobLat*0.5*ky;//Tranlation of Bobber Lattice (along X) without rotate
					TranLatY=PerBobLat*ky;//Tranlation of Bobber Lattice (along Y) without rotate
					RotTranLatX=CoChDir*TranLatX-SiChDir*TranLatY;//Tranlation of Bobber Lattice (along Y) with rotate in XY
					RotTranLatY=SiChDir*TranLatX+CoChDir*TranLatY;//Tranlation of  Lattice (along Y) with rotate in XY
					CreatBobber(ctx, px, py, pz, s, Bob_Radius, RotTranLatX, RotTranLatY, 1);
				}
			}	
		}
	break;

	case 9: // bobber_lattice_bottom
		if(ctx->chSize>0 && !(ctx->chDir[0]==0&&ctx->chDir[1]==0) )
		{

			float PerBobLat=ctx->chSize;
			float aBobLat=2*PerBobLat/(sqrt(3)); 
			float Bob_Radius=0.8*0.5*aBobLat; 
			float TranLatX;
			float TranLatY;
			float CoChDir=ctx->chDir[0]/(sqrt(ctx->chDir[0]*ctx->chDir[0]+ctx->chDir[1]*ctx->chDir[1]));//cos of rotation angle 
			float SiChDir=ctx->chDir[1]/(sqrt(ctx->chDir[0]*ctx->chDir[0]+ctx->chDir[1]*ctx->chDir[1]));//sin of rotation angle
			float RotTranLatX;
			float RotTranLatY;
		
			for (int n=0; n<ctx->NOS; n++)
			{	
				if (pz[n]<0)
				{
				VEC_X(ctx->S,n) = 0.f*ctx->Kind[n];
				VEC_Y(ctx->S,n) = 0.f*ctx->Kind[n];
				VEC_Z(ctx->S,n) = 1.f*ctx->Kind[n];
				}	
			}
			
			int Ntr=10;
			for (int kx=-Ntr; kx<=Ntr; kx++)
			{
				for (int ky=-Ntr; ky<=Ntr; ky++)
				{
					TranLatX=aBobLat*kx+aBobLat*0.5*ky;//Tranlation of Bobber Lattice (along X) without rotate
					TranLatY=PerBobLat*ky;//Tranlation of Bobber Lattice (along Y) without rotate
					RotTranLatX=CoChDir*TranLatX-SiChDir*TranLatY;//Tranlation of Bobber Lattice (along Y) with rotate in XY
					RotTranLatY=SiChDir*TranLatX+CoChDir*TranLatY;//Tranlation of  Lattice (along Y) with rotate in XY
					CreatBobber(ctx, px, py, pz, s, Bob_Radius, RotTranLatX, RotTranLatY, -1);
				}
			}	
		}	
	break;

	case 10: // hopfion H=1
		// if(chSize>0)
		// {   mat4x4 M;
		// 	mat4x4_identity(M);
		// 	// mat4x4_rotate_X(M, M, 1*PI/3);
		// 	// mat4x4_rotate_Z(M, M, 1*PI/4);
		// 	// mat4x4_rotate_X(M, M, PI/2);
		// 	float tmp;
		// 	vec3 v1,v2;
		// 	for (int n=0; n<ctx->NOS; n++)
		// 	{	
		// 		v1[0]=px[n];
		// 		v1[1]=py[n];
		// 		v1[2]=pz[n];
		// 		v2[0]=v1[0];
		// 		v2[1]=v1[1];
		// 		v2[2]=v1[2];				
		// 		r = sqrt(px[n]*px[n]+py[n]*py[n]+pz[n]*pz[n]);
		// 		if (r==0){
		// 			T = 0;
		// 		}else{
		// 			T = v2[2]/r; // cos of the angle with respect to the main axis of toroid [0,0,1]
		// 		}
		// 		T = acos(T);
		// 		t = r*r/chSize;//metka initially t = r/chSize;
		// 		t = 1.0 + 4.22/(t*t);
		// 		tmp = PI*(1.0-1.0/sqrt(t));
		// 		t = sin(tmp)*sin(T);
		// 		t = acos(1.0-2.0*t*t);
		// 		F = atan2(v2[1],v2[0]);
		// 		// f = F + atan2( 1.0/(tan(tmp)),cos(T) );
		// 		f = F - atan2( 1.0/(tan(tmp)),cos(T) );//metka +/-
		// 		VEC_X(ctx->S,n) = sin(t)*cos(f)*ctx->Kind[n];
		// 		VEC_Y(ctx->S,n) = sin(t)*sin(f)*ctx->Kind[n];
		// 		VEC_Z(ctx->S,n) = cos(t)*ctx->Kind[n];
		// 	}
		// }
		if(ctx->chSize>0)
		{
			for (int n=0; n<ctx->NOS; n++)
			// {	
			// 	r = sqrt(px[n]*px[n]+py[n]*py[n]+pz[n]*pz[n]);
			// 	if (r==0){
			// 		T = 0;
			// 	}else{
			// 		T = pz[n]/r; // angle with respect to the main axis of toroid [0,0,1]
			// 	}
			// 	T = acos(T);
			// 	F = atan2(py[n],px[n]);

			// 	t = r*r/chSize;//metka initially t = r/chSize;
			// 	t = 1.0 + 4.22/(t*t);
			// 	float tmp = PI*(1.0-1.0/sqrt(t));
			// 			// Filipp anzats 
			// 			// t = 2*r/chSize;
			// 			// t = t*t;
			// 			// tmp = PI/(1.0-t/exp(t));
			// 	t = sin(tmp)*sin(T);
			// 	t = acos(1.0-2.0*t*t);
			// 	f = -F + atan2( 1.0/(tan(tmp)),cos(T) );
			// 	// f = -F + atan( 1.0/( tan(tmp)*cos(T) ) );
			// 	VEC_X(ctx->S,n) = sin(t)*cos(f)*ctx->Kind[n];
			// 	VEC_Y(ctx->S,n) = sin(t)*sin(f)*ctx->Kind[n];
			// 	VEC_Z(ctx->S,n) = cos(t)*ctx->Kind[n];
			// }
			// metka Vlad anzats
			{	
				float rx=px[n]/ctx->chSize;
				float ry=py[n]/ctx->chSize;
				float rz=pz[n]/ctx->chSize;
				r = sqrt(rx*rx+ry*ry+rz*rz);
				float G = 2*atan(exp(-2*r)/r);
				VEC_X(ctx->S,n) = (1-cos(2*G))*rz*rx/(r*r) + sin(2*G)*ry/r;
				VEC_Y(ctx->S,n) = -(1-cos(2*G))*rz*ry/(r*r) + sin(2*G)*rx/r;
				VEC_Z(ctx->S,n) = (1-cos(2*G))*rz*rz/(r*r) + cos(2*G);
			}				
		}
	break;

	case 11: // spiral
		if(fabs(ctx->chDir[0])+fabs(ctx->chDir[1])+fabs(ctx->chDir[2])!=0)
		{	

		float angle;//in degrees
		double tmp[3];
		vec3_norm(ctx->chDir, ctx->chDir);
		//metka spiralization
		for (int i=0; i<ctx->NOS; i++)
		{	
			angle=360*(px[i]*ctx->chDir[0]+py[i]*ctx->chDir[1]+pz[i]*ctx->chDir[2])/ctx->chSize;
			RotateVector(VEC_X(ctx->S,i), VEC_Y(ctx->S,i), VEC_Z(ctx->S,i), ctx->chDir[0], ctx->chDir[1], ctx->chDir[2], angle, tmp);
			VEC_X(ctx->bS,i) = VEC_X(ctx->S,i) = tmp[0];
			VEC_Y(ctx->bS,i) = VEC_Y(ctx->S,i) = tmp[1];
			VEC_Z(ctx->bS,i) = VEC_Z(ctx->S,i) = tmp[2];
		}	
		// float Psp = 1;
		// for (int i=0; i<ctx->NOS; i++)
		// {	
		// 	float r = sqrt(px[i]*px[i]+py[i]*py[i])/64;
		// 	T = 2*PI*Psp*(px[i]+py[i])/(64*sqrt(2));
		// 	bSx[i] = VEC_X(ctx->S,i) = -sin(T)/sqrt(2);
		// 	bSy[i] = VEC_Y(ctx->S,i) = sin(T)/sqrt(2);
		// 	bSz[i] = VEC_Z(ctx->S,i) = cos(T);
		// }	

		// for (int i=0; i<ctx->NOS; i++)
		// {	
		// 	float r = (px[i]+0.5*ctx->uABC[0])/64;
		// 	T = 2*PI*(1-r);
		// 	bSx[i] = VEC_X(ctx->S,i) = 0;
		// 	bSy[i] = VEC_Y(ctx->S,i) = sin(T);
		// 	bSz[i] = VEC_Z(ctx->S,i) = cos(T);
		// }	


	//metka
			/*
			float tmpv[3];//temporal vector
			float tmpv2[3];//temporal vector
			//float F = 0.f;//initial phase
			float prp[3];// vector perpendicular to the chDir
			(void)Unitf(chDir,chDir);//normalize chDir vector

			tmpv[0] = 0;//x
			tmpv[1] = 0;//y
			tmpv[2] = 1;//z

			prp[0] = 0; prp[1] = 1; prp[2] = 0;
			if (1-Dotf(chDir, prp)<1e-6f){
				tmpv[0] = 1;//x
				tmpv[1] = 0;//y
				tmpv[2] = 0;//z
			}

			prp[0] = 0; prp[1] = 0; prp[2] = 1;
			if (1-Dotf(chDir, prp)<1e-6f){
				tmpv[0] = 0;//x
				tmpv[1] = 1;//y
				tmpv[2] = 0;//z
			}

			Crossf(tmpv, chDir, prp);
			(void)Unitf(prp,prp);

			for (int n=0; n<ctx->NOS; n++)
			{	
				tmpv[0] = px[n];
				tmpv[1] = py[n];
				tmpv[2] = pz[n];
				r = Dotf(chDir, tmpv); //distance
				T = TPI*r/chSize; //rotation angle
				tmpv[0] = cos(T);
				tmpv[1] = sin(T);
				tmpv[2] = 0.0f;
				NewBasisCartesian(tmpv, chDir, tmpv2);
				VEC_X(ctx->S,n) = tmpv2[0]*ctx->Kind[n];
				VEC_Y(ctx->S,n) = tmpv2[1]*ctx->Kind[n];
				VEC_Z(ctx->S,n) = tmpv2[2]*ctx->Kind[n];
			}
			*/
		}
	break;

	case 12: // skyrmion Lattice
		if(ctx->chSize>0)
		{
			for (int i=0; i<ctx->NOS; i++)
			{	
				VEC_X(ctx->bS,i) = VEC_X(ctx->S,i) = 0;
				VEC_Y(ctx->bS,i) = VEC_Y(ctx->S,i) = 0;
				VEC_Z(ctx->bS,i) = VEC_Z(ctx->S,i) = 1;
			}	

			// GetSkyrmion(ctx, px, py, pz, s, chSize, 0, 0);
			// GetSkyrmion(ctx, px, py, pz, s, chSize, -0.5*ctx->uABC[0], -0.5*ctx->uABC[1]);
			// GetSkyrmion(ctx, px, py, pz, s, chSize, 0.5*ctx->uABC[0], 0.5*ctx->uABC[1]);
			// GetSkyrmion(ctx, px, py, pz, s, chSize, 0.5*ctx->uABC[0], -0.5*ctx->uABC[1]);
			// GetSkyrmion(ctx, px, py, pz, s, chSize, -0.5*ctx->uABC[0], 0.5*ctx->uABC[1]);

			// ask lattice

			GetAntiskyrmion(ctx, px, py, pz, s, ctx->chSize, 0, 0,0);
			GetAntiskyrmion(ctx, px, py, pz, s, ctx->chSize, -0.5*ctx->uABC[0], -0.5*ctx->uABC[1],1);
			GetAntiskyrmion(ctx, px, py, pz, s, ctx->chSize, 0.5*ctx->uABC[0], 0.5*ctx->uABC[1],1);
			GetAntiskyrmion(ctx, px, py, pz, s, ctx->chSize, 0.5*ctx->uABC[0], -0.5*ctx->uABC[1],1);
			GetAntiskyrmion(ctx, px, py, pz, s, ctx->chSize, -0.5*ctx->uABC[0], 0.5*ctx->uABC[1],1);

			
			float T = ctx->chDir[0];
			float nx,ny,nz;
			for (int n=0; n<ctx->NOS; n++)
			{	

				int kz = (int)n/(ctx->uABC[0]*ctx->uABC[1]);
				float F = 2*PI*kz/ctx->chDir[1];
				nx = VEC_X(s,n)*cos(F) - VEC_Y(s,n)*sin(F);
				ny = VEC_X(s,n)*sin(F) + VEC_Y(s,n)*cos(F);

				VEC_X(s,n) = nx; VEC_Y(s,n)=ny;

				nx = VEC_X(s,n)*cos(T) - VEC_Z(s,n)*sin(T);
				nz = VEC_X(s,n)*sin(T) + VEC_Z(s,n)*cos(T);

				VEC_X(s,n) = nx; VEC_Z(s,n)=nz;

				nx = VEC_X(s,n)*cos(F) - VEC_Y(s,n)*sin(F);
				ny = VEC_X(s,n)*sin(F) + VEC_Y(s,n)*cos(F);

				VEC_X(s,n)=nx; VEC_Y(s,n)=ny;

				// VEC_X(ctx->S,n) = sin(acos(BextDCMagnitude/(ctx->Dij[0]*ctx->Dij[0])))*cos(pz[n]*2*PI/128);
				// VEC_Y(ctx->S,n) = sin(acos(BextDCMagnitude/(ctx->Dij[0]*ctx->Dij[0])))*sin(pz[n]*2*PI/128);
				// VEC_Z(ctx->S,n) = BextDCMagnitude/(ctx->Dij[0]*ctx->Dij[0]);
				
			
			}	

		
			// float PerSkirmLat=chSize;
			// float aSkirmLat=2*PerSkirmLat/(sqrt(3)); 
			// float Sk_Radius=0.8*0.5*aSkirmLat; 
			// float TranLatX;
			// float TranLatY;
			// float CoChDir=chDir[0]/(sqrt(chDir[0]*chDir[0]+chDir[1]*chDir[1]));//cos of rotation angle 
			// float SiChDir=chDir[1]/(sqrt(chDir[0]*chDir[0]+chDir[1]*chDir[1]));//sin of rotation angle
			// float RotTranLatX;
			// float RotTranLatY;
		
			// for (int n=0; n<ctx->NOS; n++)
			// {	
			// 	VEC_X(ctx->S,n) = 0.f*ctx->Kind[n];
			// 	VEC_Y(ctx->S,n) = 0.f*ctx->Kind[n];
			// 	VEC_Z(ctx->S,n) = 1.f*ctx->Kind[n];	
			// }
			
			// int Ntr=10;
			// for (int kx=-Ntr; kx<=Ntr; kx++)
			// {
			// 	for (int ky=-Ntr; ky<=Ntr; ky++)
			// 	{
			// 		TranLatX=aSkirmLat*kx+aSkirmLat*0.5*ky;//Tranlation of skirmion Lattice (along X) without rotate
			// 		TranLatY=PerSkirmLat*ky;//Tranlation of skirmion Lattice (along Y) without rotate
			// 		RotTranLatX=CoChDir*TranLatX-SiChDir*TranLatY;//Tranlation of skirmion Lattice (along Y) with rotate in XY
			// 		RotTranLatY=SiChDir*TranLatX+CoChDir*TranLatY;//Tranlation of skirmion Lattice (along Y) with rotate in XY
			// 		CreatSkyrmion(ctx, px, py, pz, s, Sk_Radius, RotTranLatX, RotTranLatY);
			// 	}
			// }	
		}
	break;

	case 13: // globula
		{
			for (int n=0; n<ctx->NOS; n++)
			{	
				VEC_X(ctx->S,n) = 0.f*ctx->Kind[n];
				VEC_Y(ctx->S,n) = 0.f*ctx->Kind[n];
				VEC_Z(ctx->S,n) = 1.f*ctx->Kind[n];	
			}
			CreatGlobule(ctx, px, py, pz, s, ctx->chSize, 0, 0, -1);

			//hex Lx/Ly=sqrt(3)
			// CreatGlobule(ctx, px, py, pz, s, chSize,-ctx->uABC[0]/2,-ctx->uABC[1]/2,-ctx->uABC[2]/2);
			// CreatGlobule(ctx, px, py, pz, s, chSize, ctx->uABC[0]/2,-ctx->uABC[1]/2,-ctx->uABC[2]/2);
			// CreatGlobule(ctx, px, py, pz, s, chSize,-ctx->uABC[0]/2, ctx->uABC[1]/2,-ctx->uABC[2]/2);
			// CreatGlobule(ctx, px, py, pz, s, chSize, ctx->uABC[0]/2, ctx->uABC[1]/2,-ctx->uABC[2]/2);

			// CreatGlobule(ctx, px, py, pz, s, chSize,         0,         0,-ctx->uABC[2]/2);

			// CreatGlobule(ctx, px, py, pz, s, chSize,-ctx->uABC[0]/2,-ctx->uABC[1]/2, ctx->uABC[2]/2);
			// CreatGlobule(ctx, px, py, pz, s, chSize, ctx->uABC[0]/2,-ctx->uABC[1]/2, ctx->uABC[2]/2);
			// CreatGlobule(ctx, px, py, pz, s, chSize,-ctx->uABC[0]/2, ctx->uABC[1]/2, ctx->uABC[2]/2);
			// CreatGlobule(ctx, px, py, pz, s, chSize, ctx->uABC[0]/2, ctx->uABC[1]/2, ctx->uABC[2]/2);

			// CreatGlobule(ctx, px, py, pz, s, chSize,         0,         0, ctx->uABC[2]/2);

			// CreatGlobule(ctx, px, py, pz, s, chSize,         0,-ctx->uABC[0]*sqrt(3)/3, 0);
			// CreatGlobule(ctx, px, py, pz, s, chSize,-ctx->uABC[0]/2, ctx->uABC[0]*sqrt(3)/6, 0);
			// CreatGlobule(ctx, px, py, pz, s, chSize, ctx->uABC[0]/2, ctx->uABC[0]*sqrt(3)/6, 0);

			//bcc
			// CreatGlobule(ctx, px, py, pz, s, chSize,-ctx->uABC[0]/2,-ctx->uABC[1]/2,-ctx->uABC[2]/2);
			// CreatGlobule(ctx, px, py, pz, s, chSize, ctx->uABC[0]/2,-ctx->uABC[1]/2,-ctx->uABC[2]/2);
			// CreatGlobule(ctx, px, py, pz, s, chSize,-ctx->uABC[0]/2, ctx->uABC[1]/2,-ctx->uABC[2]/2);
			// CreatGlobule(ctx, px, py, pz, s, chSize, ctx->uABC[0]/2, ctx->uABC[1]/2,-ctx->uABC[2]/2);

			// CreatGlobule(ctx, px, py, pz, s, chSize,-ctx->uABC[0]/2,-ctx->uABC[1]/2, ctx->uABC[2]/2);
			// CreatGlobule(ctx, px, py, pz, s, chSize, ctx->uABC[0]/2,-ctx->uABC[1]/2, ctx->uABC[2]/2);
			// CreatGlobule(ctx, px, py, pz, s, chSize,-ctx->uABC[0]/2, ctx->uABC[1]/2, ctx->uABC[2]/2);
			// CreatGlobule(ctx, px, py, pz, s, chSize, ctx->uABC[0]/2, ctx->uABC[1]/2, ctx->uABC[2]/2);

			// CreatGlobule(ctx, px, py, pz, s, chSize,         0,         0,         0);


			//fcc	
			// CreatGlobule(ctx, px, py, pz, s, chSize,-ctx->uABC[0]/2,-ctx->uABC[1]/2,-ctx->uABC[2]/2);
			// CreatGlobule(ctx, px, py, pz, s, chSize, ctx->uABC[0]/2,-ctx->uABC[1]/2,-ctx->uABC[2]/2);
			// CreatGlobule(ctx, px, py, pz, s, chSize,-ctx->uABC[0]/2, ctx->uABC[1]/2,-ctx->uABC[2]/2);
			// CreatGlobule(ctx, px, py, pz, s, chSize, ctx->uABC[0]/2, ctx->uABC[1]/2,-ctx->uABC[2]/2);

			// CreatGlobule(ctx, px, py, pz, s, chSize,-ctx->uABC[0]/2,-ctx->uABC[1]/2, ctx->uABC[2]/2);
			// CreatGlobule(ctx, px, py, pz, s, chSize, ctx->uABC[0]/2,-ctx->uABC[1]/2, ctx->uABC[2]/2);
			// CreatGlobule(ctx, px, py, pz, s, chSize,-ctx->uABC[0]/2, ctx->uABC[1]/2, ctx->uABC[2]/2);
			// CreatGlobule(ctx, px, py, pz, s, chSize, ctx->uABC[0]/2, ctx->uABC[1]/2, ctx->uABC[2]/2);

			// CreatGlobule(ctx, px, py, pz, s, chSize,-ctx->uABC[0]/2,         0,         0);
			// CreatGlobule(ctx, px, py, pz, s, chSize, ctx->uABC[0]/2,         0,         0);
			// CreatGlobule(ctx, px, py, pz, s, chSize,         0, ctx->uABC[1]/2,         0);
			// CreatGlobule(ctx, px, py, pz, s, chSize,         0,-ctx->uABC[1]/2,         0);

			// CreatGlobule(ctx, px, py, pz, s, chSize,         0,         0, ctx->uABC[2]/2);
			// CreatGlobule(ctx, px, py, pz, s, chSize,         0,         0,-ctx->uABC[2]/2);

			// float F0 = PI/2;
			// float P = TPI*ctx->Jij[0]/ctx->Dij[0];
			// float KZ = 0.f;
			// float Theta = 0;
			// float b=ctx->Jij[0]/(ctx->Dij[0]*sin(F0));
			// for (int n=0; n<ctx->NOS; n++)
			// {
			// 	rx=px[n];
			// 	ry=py[n];
			// 	rz=pz[n]-(ctx->uABC[2]*0.5-chSize);
			// 	r1 = sqrt(rx*rx+ry*ry);
			// 	r2 = sqrt(rx*rx+ry*ry+rz*rz);
			// 	KZ = 2*PI*rz/P;

			// 	//F = atan2(ry,rx)+F0;
			// 	//T = 2*atan(1/(r2/chSize+1)/tan(acos(rz/r2)/2));
			// 	if(rz>0){
			// 		F = atan2(ry,rx) + F0 + chDir[2]*rz*exp((rz-chSize)/chSize)/chSize;			
			// 	}else{
			// 		F = atan2(ry,rx)+F0;				
			// 	}

			// 	T = 2*atan(chDir[0]/(r2/b+1)/tan(acos(rz/r2)/2));
			// 	Theta = acos(BextDCMagnitude*ctx->Jij[0]/(ctx->Dij[0]*ctx->Dij[0])) * (1-T/PI) * exp(-T*chDir[1]);

			// 	rx = sin(T)*cos(F);
			// 	ry = sin(T)*sin(F);
			// 	rz = cos(T);
			// 	VEC_X(ctx->S,n) = ((cos(KZ)*cos(KZ)+cos(Theta)*sin(KZ)*sin(KZ))*rx + sin(2*KZ)*sin(Theta/2)*sin(Theta/2)*ry + sin(KZ)*sin(Theta)*rz)*ctx->Kind[n];
			// 	VEC_Y(ctx->S,n) = (sin(2*KZ)*sin(Theta/2)*sin(Theta/2)*rx+(cos(KZ)*cos(KZ)*cos(Theta)+sin(KZ)*sin(KZ))*ry-cos(KZ)*sin(Theta)*rz)*ctx->Kind[n];
			// 	VEC_Z(ctx->S,n) = (-sin(KZ)*sin(Theta)*rx+cos(KZ)*sin(Theta)*ry+cos(Theta)*rz)*ctx->Kind[n];

			// 	// VEC_X(ctx->S,n) = sin(T)*cos(F);
			// 	// VEC_Y(ctx->S,n) = sin(T)*sin(F);
			// 	// VEC_Z(ctx->S,n) = cos(T);
			// }
		}
	break;

	case 14: // 3D lattice of Bloch points
		if(ctx->chSize>0)
		{
			float T1 = 0.f;
			// float T2 = 0.f;
			float F1 = 0.f;
			// float F2 = 0.f;


			/*for (int n=0; n<ctx->NOS; n++)
			{	
				rx = px[n];
				ry = py[n];
				rz = pz[n];
				T1 = (rx) * TPI / chSize;
				T2 = (ry) * TPI / chSize;
				F1 = (rx) * TPI / chSize;
				F2 = (ry) * TPI / chSize;
				VEC_X(ctx->S,n) = 0;
				VEC_Y(ctx->S,n) = 0;
				VEC_Z(ctx->S,n) = 1;
				//kanazawa
				VEC_X(ctx->S,n) = sin(ry * TPI / chSize)+cos(rz * TPI / chSize);
				VEC_Y(ctx->S,n) = sin(rz * TPI / chSize)+cos(rx * TPI / chSize);
				VEC_Z(ctx->S,n) = sin(rx * TPI / chSize)+cos(ry * TPI / chSize);

				rx=sqrt(VEC_X(ctx->S,n)*VEC_X(ctx->S,n)+VEC_Y(ctx->S,n)*VEC_Y(ctx->S,n)+VEC_Z(ctx->S,n)*VEC_Z(ctx->S,n));
				VEC_X(ctx->S,n) = VEC_X(ctx->S,n)/rx;
				VEC_Y(ctx->S,n) = VEC_Y(ctx->S,n)/rx;
				VEC_Z(ctx->S,n) = VEC_Z(ctx->S,n)/rx;
				// Rx(T1,&VEC_X(ctx->S,n),&VEC_Y(ctx->S,n),&VEC_Z(ctx->S,n));
				// Ry(T2,&VEC_X(ctx->S,n),&VEC_Y(ctx->S,n),&VEC_Z(ctx->S,n));
				// Rz(rz * TPI / chSize,&VEC_X(ctx->S,n),&VEC_Y(ctx->S,n),&VEC_Z(ctx->S,n));

				// VEC_X(ctx->S,n) = (sin(T1)+sin(T2)) * (cos(F1-F2)) * ctx->Kind[n];
				// VEC_Y(ctx->S,n) = (sin(T1)+sin(T2)) * (sin(F1+F2)) * ctx->Kind[n];
				// VEC_Z(ctx->S,n) = (cos(T1)-cos(T2)) * ctx->Kind[n];
			}*/
			float Lx=px[ctx->uABC[0]-1]-px[0];
			float Ly=py[ctx->uABC[1]*ctx->uABC[0]-1]-py[0];
			for (int n=0; n<ctx->NOS; n++)
			{	
				float rx = px[n]-px[0];
				float ry = py[n]-py[0];
				// float rz = pz[n]-pz[0];
				F1 = rx * TPI / Lx;
				T1 = ry * PI / Ly;
				VEC_X(ctx->S,n) = 0;
				VEC_Y(ctx->S,n) = 0;
				VEC_Z(ctx->S,n) = 1;
				VEC_X(ctx->S,n) = sin(T1)*cos(F1)*ctx->Kind[n];
				VEC_Y(ctx->S,n) = sin(T1)*sin(F1)*ctx->Kind[n];
				VEC_Z(ctx->S,n) = cos(T1)*ctx->Kind[n];
			}			
		}
		printf("spin0x=%f",px[0]);
	break;
	
	case 15:; // Normalize all spins
	double inv_abs_S;
	for (int n=0; n<ctx->NOS; n++)
			{	
				if (ctx->Kind[n]!=0){
					inv_abs_S=1.0f/sqrt(VEC_X(ctx->S,n)*VEC_X(ctx->S,n)+VEC_Y(ctx->S,n)*VEC_Y(ctx->S,n)+VEC_Z(ctx->S,n)*VEC_Z(ctx->S,n));
					VEC_X(ctx->S,n) = VEC_X(ctx->S,n)*inv_abs_S;
					VEC_Y(ctx->S,n) = VEC_Y(ctx->S,n)*inv_abs_S;
					VEC_Z(ctx->S,n) = VEC_Z(ctx->S,n)*inv_abs_S;
				}	
			}
	break;
	}	
}
