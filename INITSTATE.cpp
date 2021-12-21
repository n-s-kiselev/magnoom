void CreatSkyrmion(float * px, float * py, float * pz, double * sx, double * sy, double * sz, float Sk_R, float tx, float ty)
{
	float T = 0.f;
	float F = 0.f;
	float r = 0.f;
	float rx=0.0f;
	float ry=0.0f;

	for (int n=0; n<NOS; n++)
	{
		rx=px[n]-tx;
		ry=py[n]-ty;
		r = sqrt(rx*rx+ry*ry);
		// if (r<Sk_R)//metka
		{
			T= PI*exp(-0.5*r/Sk_R);//<-- defines skyrmion profile you may put periodical function to get target like skyrmions
			F= atan2(ry,rx)+PI*0.5;//<--chiral (bloch) skyrmion |Q|=1
			Sx[n] = sin(T)*cos(F)*Kind[n];
			Sy[n] = sin(T)*sin(F)*Kind[n];
			Sz[n] = cos(T)*Kind[n];
			//metka test  color code:	
			// float inv_abs_S=1.0f/sqrt(Sx[n]*Sx[n]+Sy[n]*Sy[n]+Sz[n]*Sz[n]);
			// Sx[n] = Sx[n]*inv_abs_S/(0.1*r+1.0);
			// Sy[n] = Sy[n]*inv_abs_S/(0.1*r+1.0);
			// Sz[n] = Sz[n]*inv_abs_S/(0.1*r+1.0);
		}
	}
}

void GetSkyrmion(float * px, float * py, float * pz, double * sx, double * sy, double * sz, float Sk_R, float tx, float ty)
{
	float T = 0.f;
	float F = 0.f;
	float r = 0.f;
	float rx=0.0f;
	float ry=0.0f;



	for (int n=0; n<NOS; n++)
	{
		rx=px[n]-tx;
		ry=py[n]-ty;
		r = sqrt(rx*rx+ry*ry);
		if (r<Sk_R){
			T = PI*(1-r/Sk_R);
			F = atan2(ry,rx)+PI*0.5;
			sx[n] = sin(T)*cos(F)*Kind[n];
			sy[n] = sin(T)*sin(F)*Kind[n];
			sz[n] = cos(T)*Kind[n];
		}
	}
}

void GetAntiskyrmion(float * px, float * py, float * pz, double * sx, double * sy, double * sz, float Sk_R, float tx, float ty, bool orient)
{
	float T = 0.f;
	float F = 0.f;
	float r = 0.f;
	float rx=0.0f;
	float ry=0.0f;



	for (int n=0; n<NOS; n++)
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
			
			sx[n] = sin(T)*cos(F)*Kind[n];
			sy[n] = sin(T)*sin(F)*Kind[n];
			sz[n] = cos(T)*Kind[n];
		}
	}
}

void TiltSpinsToX(float * px, float * py, float * pz, double * sx, double * sy, double * sz, float Sk_R, float tx, float ty)
{
	float T = 0.f;
	float F = 0.f;
	float r = 0.f;
	float rx=0.0f;
	float ry=0.0f;

	for (int n=0; n<NOS; n++)
	{
		rx=px[n]-tx;
		ry=py[n]-ty;
		r = sqrt(rx*rx+ry*ry);
		if (r<Sk_R)
		{
			T= PI/2;
			F= 0;
			Sx[n] = sin(T)*cos(F)*Kind[n];
			Sy[n] = sin(T)*sin(F)*Kind[n];
			Sz[n] = cos(T)*Kind[n];
		}
	}
}

void CreatSkyrmionSoliton(float * px, float * py, float * pz, double * sx, double * sy, double * sz, float Sk_R, float tx, float ty)
{
	float T = 0.f;
	float F = 0.f;
	float A = 0.f;
	float a = chDir[0];//sqrt(alpha)
	float v = chDir[1];
	float r = 0.f;
	float rx= 0.f;
	float ry= 0.f;
	float rz= 0.f;
	for (int n=0; n<NOS; n++){
		rx=px[n];
		ry=py[n];
		rz=pz[n];

		F = 2*a*tanh(a*v*rz*chDir[2])+v*rz*chDir[2];
		A = 10*2*a/(fabs(v)*cosh(a*v*rz*chDir[2]));
		rx= px[n]-A*cos(F);
		ry= py[n]-A*sin(F);
		r = sqrt(rx*rx+ry*ry);
		if (r<Sk_R){
			T = PI*exp(-2*r/Sk_R);
			F = atan2(ry,rx)+PI*0.5;
			Sx[n] = sin(T)*cos(F)*Kind[n];
			Sy[n] = sin(T)*sin(F)*Kind[n];
			Sz[n] = cos(T)*Kind[n];			
		}else{
			Sx[n] = 0.0f;
			Sy[n] = 0.0f;
			Sz[n] = 1.0f*Kind[n];
		}
	}
}


void CreatGlobule(float * px, float * py, float * pz, double * sx, double * sy, double * sz, float Sk_R, float tx, float ty, float tz)
{
	float T = 0.f;
	float F = 0.f;
	float rx=0.0f;
	float ry=0.0f;
	float rz=0.0f;
	float r1,r2;

	for (int n=0; n<NOS; n++)
	{
		rx=px[n]-tx;
		ry=py[n]-ty;
		rz=pz[n]-tz;
		r1 = sqrt(rx*rx+ry*ry);
		r2 = sqrt(rx*rx+ry*ry+rz*rz);
		if (r2<chSize*0.5)
		{
			T=PI*exp(-2*r1/chSize);//<-- defines skyrmion profile you may put periodical function to get target like skyrmions
			F= atan2(ry,rx)+PI*0.5;//<--chiral (bloch) skyrmion |Q|=1
			Sx[n] = sin(T)*cos(F)*Kind[n];
			Sy[n] = sin(T)*sin(F)*Kind[n];
			Sz[n] = cos(T)*Kind[n];			
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

void CreatBobber(float * px, float * py, float * pz, double * sx, double * sy, double * sz, float Bob_R ,float tx, float ty, int Loc_B)
{
	float top = 0.5*Box[2][2];
	float R;
	float T = 0.f;
	float F = 0.f;
	float r = 0.f;
	//float Bob_R = chSize;
	float Bob_D = 2*Bob_R;
	float rx=0.0f;
	float ry=0.0f;
	if(chSize>0)
	{
		if (Bob_D>top) { Bob_D = 0.75*top;}

		for (int n=0; n<NOS; n++)
		{		
			rx=px[n]-tx;
			ry=py[n]-ty;			
			R = Bob_R*(1-(top-Loc_B*pz[n])/(Bob_D));//If the bobber set in positive Loc_B=1 side else bobber set in negative Loc_B=-1 
			r = sqrt(rx*rx+ry*ry);
			if ( (r<R) && (fabs(pz[n])>(top-Bob_D)) )
			{
			T=PI*(R-r)/R;//<-- defines skyrmion profile you may put periodical function to get target like skyrmions
			F= 1*(-atan2(rx,ry));//<--chiral (bloch) skyrmion |Q|=1
			Sx[n] = -sin(T)*cos(F)*Kind[n];
			Sy[n] = -sin(T)*sin(F)*Kind[n];
			Sz[n] = cos(T)*Kind[n];	
			}
		}	
	}
}

void CreatHorisontalBobber(float * px, float * py, float * pz, double * sx, double * sy, double * sz, float Bob_R){
	for(int i = 0; i < uABC[0]; i++){
		for(int j = 0; j < uABC[1]; j++){
			for(int k = 0; k < uABC[2]; k++){
				int n = i + j*uABC[0] + k*uABC[0]*uABC[1];
				Sx[n] = 0;
				Sy[n] = 0;
				Sz[n] = 1;	
				if(px[n]<0 && (py[n]*py[n]+pz[n]*pz[n]) < Bob_R*Bob_R){
					float r = sqrt(py[n]*py[n]+pz[n]*pz[n]);
					float T = PI*(Bob_R-r)/Bob_R;
					float F = atan2(py[n],pz[n]) + 0.5*PI;
					Sx[n] = sin(T)*cos(F);
					Sy[n] = sin(T)*sin(F);
					Sz[n] = cos(T);	
				}
			}
		}
	}
}

void InitSpinComponents(float * px, float * py, float * pz, double * sx, double * sy, double * sz, int id)
{
	float T = 0.f;
	float F = 0.f;
	// T, F are theta and phi for spherical spase coordinates 
	float t = 0.f;
	float f = 0.f;
	// t, f are theta and phi for spin components 
	float r = 0.f;

	switch (id)
	{

	case 0: //random spins
		double rnd[3];
		for (int n=0; n<NOS; n++)
		// for (int n=0; n<uABC[0]; n++)//metka
		{	
			rnd[0] = 2.0 * (0.5 - rand() / (double)RAND_MAX);
			rnd[1] = 2.0 * (0.5 - rand() / (double)RAND_MAX);
			rnd[2] = 2.0 * (0.5 - rand() / (double)RAND_MAX);
			(void)Unit(rnd,rnd);
			Sx[n] = rnd[0]*Kind[n];
			Sy[n] = rnd[1]*Kind[n];
			Sz[n] = rnd[2]*Kind[n];	
		}	
		//test
			// Sx[0] = 1.0;
			// Sy[0] = 0.0;
			// Sz[0] = 0.0;
			// Sx[1] = 0.0;
			// Sy[1] = 0.0;
			// Sz[1] = 1.0;	
		// Stefan		
		// if(chSize>0)
		// {
		// 	if( rand() / (double)RAND_MAX>0.5){
		// 		TiltSpinsToX(px, py, pz, sx, sy, sz, chSize, 0.5*uABC[0], 0);				
		// 	}

		// 	if( 2.0 * (0.5 - rand() / (double)RAND_MAX)>0){
		// 		TiltSpinsToX(px, py, pz, sx, sy, sz, chSize,-0.5*uABC[0], 0);				
		// 	}

		// 	// if( 2.0 * (0.5 - rand() / (double)RAND_MAX)>0){
		// 	// 	TiltSpinsToX(px, py, pz, sx, sy, sz, chSize, 0.5*uABC[0], -0.5*uABC[1]);				
		// 	// }

		// 	// if( 2.0 * (0.5 - rand() / (double)RAND_MAX)>0){
		// 	// 	TiltSpinsToX(px, py, pz, sx, sy, sz, chSize, 0.5*uABC[0],  0.5*uABC[1]);				
		// 	// }
		// }
	break;

	case 1: //homogeneous
		for (int n=0; n<NOS; n++)
		{	
			Sx[n] = chDir[0]*Kind[n];
			Sy[n] = chDir[1]*Kind[n];
			Sz[n] = chDir[2]*Kind[n];	

			//tilted FM
			// Sx[n] = sin(VHtheta*PI/180)*cos(VHphi*PI/180);
			// Sy[n] = sin(VHtheta*PI/180)*sin(VHphi*PI/180);
			// Sz[n] = cos(VHtheta*PI/180);	

			//cone
			// Sx[n] = sin(acos(Hf/(Dij[0]*Dij[0])))*cos(pz[n]*2*PI/64);
			// Sy[n] = sin(acos(Hf/(Dij[0]*Dij[0])))*sin(pz[n]*2*PI/64);
			// Sz[n] = Hf/(Dij[0]*Dij[0]);	
			// int kz = n%(uABC[0]*uABC[1]);
			// Sx[n] = sin(acos(0.7975))*cos(pz[n]*2*PI/128);
			// Sy[n] = sin(acos(0.7975))*sin(pz[n]*2*PI/128);
			// Sz[n] = 0.7975;	
		}

	break;

	case 2: // skyrmion Q=1
		if(chSize>0)
		{
			CreatSkyrmion(px, py, pz, sx, sy, sz, chSize,0,0);	
		}
	break;

	case 3: // skyrmion Q=2
		if(chSize>0)
		{
			// for (int n=0; n<NOS; n++)
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
			// Sx[n] = sin(T)*cos(F)*Kind[n];
			// Sy[n] = sin(T)*sin(F)*Kind[n];
			// Sz[n] = cos(T)*Kind[n];	
			// }
			//metka eto antiskyrmion
				float T = 0.f;
				float F = 0.f;
				float r = 0.f;
				float rx=0.0f;
				float ry=0.0f;

				for (int n=0; n<NOS; n++)
				{
					rx=px[n];
					ry=py[n];
					r = sqrt(rx*rx+ry*ry);
					if (r<chSize)
					{
						T= PI*exp(-2*r/chSize);//<-- defines skyrmion profile you may put periodical function to get target like skyrmions
						F= atan2(ry,rx)+PI*0.5;//<--chiral (bloch) skyrmion |Q|=1
						Sy[n] = sin(T)*cos(F)*Kind[n];
						Sx[n] = sin(T)*sin(F)*Kind[n];
						Sz[n] = cos(T)*Kind[n];
					}
				}
		}
	break;

	case 4: // skyrmion Q=3
		if(chSize>0)
		{
			for (int n=0; n<NOS; n++)
			{	
			r = sqrt(px[n]*px[n]+py[n]*py[n]+0*pz[n]*pz[n]);
			//r = sqrt(px[n]*px[n]+py[n]*py[n]);
			if (r<chSize){
			T=PI*exp(-2*r/chSize);//<-- defines skyrmion profile you may put periodical function to get target like skyrmions
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
			Sx[n] = sin(T)*cos(F)*Kind[n];
			Sy[n] = sin(T)*sin(F)*Kind[n];
			Sz[n] = cos(T)*Kind[n];	
			}
		}	
	break;

	case 5: // bobber_top
		if(chSize>0)
		
		{
			CreatBobber(px, py, pz, sx, sy, sz, chSize, 0, 0, 1);	
		}
	break;

	case 6: // bobber_bottom
		if(chSize>0)
		
		{
			// CreatBobber(px, py, pz, sx, sy, sz, chSize, 0, 0, -1);	
			CreatHorisontalBobber(px, py, pz, sx, sy, sz, chSize);			
		}
	break;

	case 7: // bobber_lattice
		if(chSize>0 && !(chDir[0]==0&&chDir[1]==0) )
		{

			float PerBobLat=chSize;
			float aBobLat=2*PerBobLat/(sqrt(3)); 
			float Bob_Radius=0.8*0.5*aBobLat; 
			float TranLatX,TranLatX1;
			float TranLatY,TranLatY1;
			float CoChDir=chDir[0]/(sqrt(chDir[0]*chDir[0]+chDir[1]*chDir[1]));//cos of rotation angle 
			float SiChDir=chDir[1]/(sqrt(chDir[0]*chDir[0]+chDir[1]*chDir[1]));//sin of rotation angle
			float RotTranLatX;
			float RotTranLatY;
		
			for (int n=0; n<NOS; n++)
			{	
				Sx[n] = 0.f*Kind[n];
				Sy[n] = 0.f*Kind[n];
				Sz[n] = 1.f*Kind[n];	
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
					CreatBobber(px, py, pz, sx, sy, sz, Bob_Radius, RotTranLatX, RotTranLatY, 1);

					TranLatX1=TranLatX+aBobLat/2;//
					TranLatY1=TranLatY+PerBobLat/3;//

					RotTranLatX=CoChDir*TranLatX1-SiChDir*TranLatY1;//Tranlation of Bobber Lattice (along Y) with rotate in XY
					RotTranLatY=SiChDir*TranLatX1+CoChDir*TranLatY1;//Tranlation of  Lattice (along Y) with rotate in XY
					CreatBobber(px, py, pz, sx, sy, sz, Bob_Radius, RotTranLatX, RotTranLatY, -1);
				}
			}	
		}
	break;

	case 8: // bobber_lattice_top
		if(chSize>0 && !(chDir[0]==0&&chDir[1]==0) )
		{

			float PerBobLat=chSize;
			float aBobLat=2*PerBobLat/(sqrt(3)); 
			float Bob_Radius=0.8*0.5*aBobLat; 
			float TranLatX;
			float TranLatY;
			float CoChDir=chDir[0]/(sqrt(chDir[0]*chDir[0]+chDir[1]*chDir[1]));//cos of rotation angle 
			float SiChDir=chDir[1]/(sqrt(chDir[0]*chDir[0]+chDir[1]*chDir[1]));//sin of rotation angle
			float RotTranLatX;
			float RotTranLatY;
		
			for (int n=0; n<NOS; n++)
			{	
				if (pz[n]>0)
				{
				Sx[n] = 0.f*Kind[n];
				Sy[n] = 0.f*Kind[n];
				Sz[n] = 1.f*Kind[n];
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
					CreatBobber(px, py, pz, sx, sy, sz, Bob_Radius, RotTranLatX, RotTranLatY, 1);
				}
			}	
		}
	break;

	case 9: // bobber_lattice_bottom
		if(chSize>0 && !(chDir[0]==0&&chDir[1]==0) )
		{

			float PerBobLat=chSize;
			float aBobLat=2*PerBobLat/(sqrt(3)); 
			float Bob_Radius=0.8*0.5*aBobLat; 
			float TranLatX;
			float TranLatY;
			float CoChDir=chDir[0]/(sqrt(chDir[0]*chDir[0]+chDir[1]*chDir[1]));//cos of rotation angle 
			float SiChDir=chDir[1]/(sqrt(chDir[0]*chDir[0]+chDir[1]*chDir[1]));//sin of rotation angle
			float RotTranLatX;
			float RotTranLatY;
		
			for (int n=0; n<NOS; n++)
			{	
				if (pz[n]<0)
				{
				Sx[n] = 0.f*Kind[n];
				Sy[n] = 0.f*Kind[n];
				Sz[n] = 1.f*Kind[n];
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
					CreatBobber(px, py, pz, sx, sy, sz, Bob_Radius, RotTranLatX, RotTranLatY, -1);
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
		// 	for (int n=0; n<NOS; n++)
		// 	{	
		// 		v1[0]=px[n];
		// 		v1[1]=py[n];
		// 		v1[2]=pz[n];
		// 		// mat4x4_mul_vec3(v2, M, v1);
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
		// 		Sx[n] = sin(t)*cos(f)*Kind[n];
		// 		Sy[n] = sin(t)*sin(f)*Kind[n];
		// 		Sz[n] = cos(t)*Kind[n];
		// 	}
		// }
		if(chSize>0)
		{
			float tmp;
			for (int n=0; n<NOS; n++)
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
			// 	tmp = PI*(1.0-1.0/sqrt(t));
			// 			// Filipp anzats 
			// 			// t = 2*r/chSize;
			// 			// t = t*t;
			// 			// tmp = PI/(1.0-t/exp(t));
			// 	t = sin(tmp)*sin(T);
			// 	t = acos(1.0-2.0*t*t);
			// 	f = -F + atan2( 1.0/(tan(tmp)),cos(T) );
			// 	// f = -F + atan( 1.0/( tan(tmp)*cos(T) ) );
			// 	Sx[n] = sin(t)*cos(f)*Kind[n];
			// 	Sy[n] = sin(t)*sin(f)*Kind[n];
			// 	Sz[n] = cos(t)*Kind[n];
			// }
			// metka Vlad anzats
			{	
				float rx=px[n]/chSize;
				float ry=py[n]/chSize;
				float rz=pz[n]/chSize;
				r = sqrt(rx*rx+ry*ry+rz*rz);
				float G = 2*atan(exp(-2*r)/r);
				Sx[n] = (1-cos(2*G))*rz*rx/(r*r) + sin(2*G)*ry/r;
				Sy[n] = -(1-cos(2*G))*rz*ry/(r*r) + sin(2*G)*rx/r;
				Sz[n] = (1-cos(2*G))*rz*rz/(r*r) + cos(2*G);
			}				
		}
	break;

	case 11: // spiral
		if(fabs(chDir[0])+fabs(chDir[1])+fabs(chDir[2])!=0)
		{	

		float angle;//in degrees
		double tmp[3];
		vec3_norm(chDir, chDir);
		//metka spiralization
		for (int i=0; i<NOS; i++)
		{	
			angle=360*(px[i]*chDir[0]+py[i]*chDir[1]+pz[i]*chDir[2])/chSize;
			RotateVector(Sx[i], Sy[i], Sz[i], chDir[0], chDir[1], chDir[2], angle, tmp);
			bSx[i] = Sx[i] = tmp[0];
			bSy[i] = Sy[i] = tmp[1];
			bSz[i] = Sz[i] = tmp[2];
		}	
		// float Psp = 1;
		// for (int i=0; i<NOS; i++)
		// {	
		// 	float r = sqrt(px[i]*px[i]+py[i]*py[i])/64;
		// 	T = 2*PI*Psp*(px[i]+py[i])/(64*sqrt(2));
		// 	bSx[i] = Sx[i] = -sin(T)/sqrt(2);
		// 	bSy[i] = Sy[i] = sin(T)/sqrt(2);
		// 	bSz[i] = Sz[i] = cos(T);
		// }	

		// for (int i=0; i<NOS; i++)
		// {	
		// 	float r = (px[i]+0.5*uABC[0])/64;
		// 	T = 2*PI*(1-r);
		// 	bSx[i] = Sx[i] = 0;
		// 	bSy[i] = Sy[i] = sin(T);
		// 	bSz[i] = Sz[i] = cos(T);
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

			for (int n=0; n<NOS; n++)
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
				Sx[n] = tmpv2[0]*Kind[n];
				Sy[n] = tmpv2[1]*Kind[n];
				Sz[n] = tmpv2[2]*Kind[n];
			}
			*/
		}
	break;

	case 12: // skyrmion Lattice
		if(chSize>0)
		{
			for (int i=0; i<NOS; i++)
			{	
				bSx[i] = Sx[i] = 0;
				bSy[i] = Sy[i] = 0;
				bSz[i] = Sz[i] = 1;
			}	

			// GetSkyrmion(px, py, pz, sx, sy, sz, chSize, 0, 0);
			// GetSkyrmion(px, py, pz, sx, sy, sz, chSize, -0.5*uABC[0], -0.5*uABC[1]);
			// GetSkyrmion(px, py, pz, sx, sy, sz, chSize, 0.5*uABC[0], 0.5*uABC[1]);
			// GetSkyrmion(px, py, pz, sx, sy, sz, chSize, 0.5*uABC[0], -0.5*uABC[1]);
			// GetSkyrmion(px, py, pz, sx, sy, sz, chSize, -0.5*uABC[0], 0.5*uABC[1]);

			// ask lattice

			GetAntiskyrmion(px, py, pz, sx, sy, sz, chSize, 0, 0,0);
			GetAntiskyrmion(px, py, pz, sx, sy, sz, chSize, -0.5*uABC[0], -0.5*uABC[1],1);
			GetAntiskyrmion(px, py, pz, sx, sy, sz, chSize, 0.5*uABC[0], 0.5*uABC[1],1);
			GetAntiskyrmion(px, py, pz, sx, sy, sz, chSize, 0.5*uABC[0], -0.5*uABC[1],1);
			GetAntiskyrmion(px, py, pz, sx, sy, sz, chSize, -0.5*uABC[0], 0.5*uABC[1],1);

			
			float T = chDir[0];
			float nx,ny,nz;
			for (int n=0; n<NOS; n++)
			{	

				int kz = (int)n/(uABC[0]*uABC[1]);
				float F = 2*PI*kz/chDir[1];
				nx = sx[n]*cos(F) - sy[n]*sin(F);
				ny = sx[n]*sin(F) + sy[n]*cos(F);

				sx[n] = nx; sy[n]=ny;

				nx = sx[n]*cos(T) - sz[n]*sin(T);
				nz = sx[n]*sin(T) + sz[n]*cos(T);

				sx[n] = nx; sz[n]=nz;

				nx = sx[n]*cos(F) - sy[n]*sin(F);
				ny = sx[n]*sin(F) + sy[n]*cos(F);

				sx[n]=nx; sy[n]=ny;

				// Sx[n] = sin(acos(Hf/(Dij[0]*Dij[0])))*cos(pz[n]*2*PI/128);
				// Sy[n] = sin(acos(Hf/(Dij[0]*Dij[0])))*sin(pz[n]*2*PI/128);
				// Sz[n] = Hf/(Dij[0]*Dij[0]);	
				
			
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
		
			// for (int n=0; n<NOS; n++)
			// {	
			// 	Sx[n] = 0.f*Kind[n];
			// 	Sy[n] = 0.f*Kind[n];
			// 	Sz[n] = 1.f*Kind[n];	
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
			// 		CreatSkyrmion(px, py, pz, sx, sy, sz, Sk_Radius, RotTranLatX, RotTranLatY);
			// 	}
			// }	
		}
	break;

	case 13: // globula
		{
			for (int n=0; n<NOS; n++)
			{	
				Sx[n] = 0.f*Kind[n];
				Sy[n] = 0.f*Kind[n];
				Sz[n] = 1.f*Kind[n];	
			}
			CreatGlobule(px, py, pz, sx, sy, sz, chSize, 0, 0, -1);

			//hex Lx/Ly=sqrt(3)
			// CreatGlobule(px, py, pz, sx, sy, sz, chSize,-uABC[0]/2,-uABC[1]/2,-uABC[2]/2);
			// CreatGlobule(px, py, pz, sx, sy, sz, chSize, uABC[0]/2,-uABC[1]/2,-uABC[2]/2);
			// CreatGlobule(px, py, pz, sx, sy, sz, chSize,-uABC[0]/2, uABC[1]/2,-uABC[2]/2);
			// CreatGlobule(px, py, pz, sx, sy, sz, chSize, uABC[0]/2, uABC[1]/2,-uABC[2]/2);

			// CreatGlobule(px, py, pz, sx, sy, sz, chSize,         0,         0,-uABC[2]/2);

			// CreatGlobule(px, py, pz, sx, sy, sz, chSize,-uABC[0]/2,-uABC[1]/2, uABC[2]/2);
			// CreatGlobule(px, py, pz, sx, sy, sz, chSize, uABC[0]/2,-uABC[1]/2, uABC[2]/2);
			// CreatGlobule(px, py, pz, sx, sy, sz, chSize,-uABC[0]/2, uABC[1]/2, uABC[2]/2);
			// CreatGlobule(px, py, pz, sx, sy, sz, chSize, uABC[0]/2, uABC[1]/2, uABC[2]/2);

			// CreatGlobule(px, py, pz, sx, sy, sz, chSize,         0,         0, uABC[2]/2);

			// CreatGlobule(px, py, pz, sx, sy, sz, chSize,         0,-uABC[0]*sqrt(3)/3, 0);
			// CreatGlobule(px, py, pz, sx, sy, sz, chSize,-uABC[0]/2, uABC[0]*sqrt(3)/6, 0);
			// CreatGlobule(px, py, pz, sx, sy, sz, chSize, uABC[0]/2, uABC[0]*sqrt(3)/6, 0);

			//bcc
			// CreatGlobule(px, py, pz, sx, sy, sz, chSize,-uABC[0]/2,-uABC[1]/2,-uABC[2]/2);
			// CreatGlobule(px, py, pz, sx, sy, sz, chSize, uABC[0]/2,-uABC[1]/2,-uABC[2]/2);
			// CreatGlobule(px, py, pz, sx, sy, sz, chSize,-uABC[0]/2, uABC[1]/2,-uABC[2]/2);
			// CreatGlobule(px, py, pz, sx, sy, sz, chSize, uABC[0]/2, uABC[1]/2,-uABC[2]/2);

			// CreatGlobule(px, py, pz, sx, sy, sz, chSize,-uABC[0]/2,-uABC[1]/2, uABC[2]/2);
			// CreatGlobule(px, py, pz, sx, sy, sz, chSize, uABC[0]/2,-uABC[1]/2, uABC[2]/2);
			// CreatGlobule(px, py, pz, sx, sy, sz, chSize,-uABC[0]/2, uABC[1]/2, uABC[2]/2);
			// CreatGlobule(px, py, pz, sx, sy, sz, chSize, uABC[0]/2, uABC[1]/2, uABC[2]/2);

			// CreatGlobule(px, py, pz, sx, sy, sz, chSize,         0,         0,         0);


			//fcc	
			// CreatGlobule(px, py, pz, sx, sy, sz, chSize,-uABC[0]/2,-uABC[1]/2,-uABC[2]/2);
			// CreatGlobule(px, py, pz, sx, sy, sz, chSize, uABC[0]/2,-uABC[1]/2,-uABC[2]/2);
			// CreatGlobule(px, py, pz, sx, sy, sz, chSize,-uABC[0]/2, uABC[1]/2,-uABC[2]/2);
			// CreatGlobule(px, py, pz, sx, sy, sz, chSize, uABC[0]/2, uABC[1]/2,-uABC[2]/2);

			// CreatGlobule(px, py, pz, sx, sy, sz, chSize,-uABC[0]/2,-uABC[1]/2, uABC[2]/2);
			// CreatGlobule(px, py, pz, sx, sy, sz, chSize, uABC[0]/2,-uABC[1]/2, uABC[2]/2);
			// CreatGlobule(px, py, pz, sx, sy, sz, chSize,-uABC[0]/2, uABC[1]/2, uABC[2]/2);
			// CreatGlobule(px, py, pz, sx, sy, sz, chSize, uABC[0]/2, uABC[1]/2, uABC[2]/2);

			// CreatGlobule(px, py, pz, sx, sy, sz, chSize,-uABC[0]/2,         0,         0);
			// CreatGlobule(px, py, pz, sx, sy, sz, chSize, uABC[0]/2,         0,         0);
			// CreatGlobule(px, py, pz, sx, sy, sz, chSize,         0, uABC[1]/2,         0);
			// CreatGlobule(px, py, pz, sx, sy, sz, chSize,         0,-uABC[1]/2,         0);

			// CreatGlobule(px, py, pz, sx, sy, sz, chSize,         0,         0, uABC[2]/2);
			// CreatGlobule(px, py, pz, sx, sy, sz, chSize,         0,         0,-uABC[2]/2);

			// float F0 = PI/2;
			// float P = TPI*Jij[0]/Dij[0];
			// float KZ = 0.f;
			// float Theta = 0;
			// float b=Jij[0]/(Dij[0]*sin(F0));
			// for (int n=0; n<NOS; n++)
			// {
			// 	rx=px[n];
			// 	ry=py[n];
			// 	rz=pz[n]-(uABC[2]*0.5-chSize);
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
			// 	Theta = acos(Hf*Jij[0]/(Dij[0]*Dij[0])) * (1-T/PI) * exp(-T*chDir[1]);

			// 	rx = sin(T)*cos(F);
			// 	ry = sin(T)*sin(F);
			// 	rz = cos(T);
			// 	Sx[n] = ((cos(KZ)*cos(KZ)+cos(Theta)*sin(KZ)*sin(KZ))*rx + sin(2*KZ)*sin(Theta/2)*sin(Theta/2)*ry + sin(KZ)*sin(Theta)*rz)*Kind[n];
			// 	Sy[n] = (sin(2*KZ)*sin(Theta/2)*sin(Theta/2)*rx+(cos(KZ)*cos(KZ)*cos(Theta)+sin(KZ)*sin(KZ))*ry-cos(KZ)*sin(Theta)*rz)*Kind[n];
			// 	Sz[n] = (-sin(KZ)*sin(Theta)*rx+cos(KZ)*sin(Theta)*ry+cos(Theta)*rz)*Kind[n];

			// 	// Sx[n] = sin(T)*cos(F);
			// 	// Sy[n] = sin(T)*sin(F);
			// 	// Sz[n] = cos(T);
			// }
		}
	break;

	case 14: // 3D lattice of Bloch points
		if(chSize>0)
		{
			float T1 = 0.f;
			// float T2 = 0.f;
			float F1 = 0.f;
			// float F2 = 0.f;


			/*for (int n=0; n<NOS; n++)
			{	
				rx = px[n];
				ry = py[n];
				rz = pz[n];
				T1 = (rx) * TPI / chSize;
				T2 = (ry) * TPI / chSize;
				F1 = (rx) * TPI / chSize;
				F2 = (ry) * TPI / chSize;
				Sx[n] = 0;
				Sy[n] = 0;
				Sz[n] = 1;
				//kanazawa
				Sx[n] = sin(ry * TPI / chSize)+cos(rz * TPI / chSize);
				Sy[n] = sin(rz * TPI / chSize)+cos(rx * TPI / chSize);
				Sz[n] = sin(rx * TPI / chSize)+cos(ry * TPI / chSize);

				rx=sqrt(Sx[n]*Sx[n]+Sy[n]*Sy[n]+Sz[n]*Sz[n]);
				Sx[n] = Sx[n]/rx;
				Sy[n] = Sy[n]/rx;
				Sz[n] = Sz[n]/rx;
				// Rx(T1,&Sx[n],&Sy[n],&Sz[n]);
				// Ry(T2,&Sx[n],&Sy[n],&Sz[n]);
				// Rz(rz * TPI / chSize,&Sx[n],&Sy[n],&Sz[n]);

				// Sx[n] = (sin(T1)+sin(T2)) * (cos(F1-F2)) * Kind[n];
				// Sy[n] = (sin(T1)+sin(T2)) * (sin(F1+F2)) * Kind[n];
				// Sz[n] = (cos(T1)-cos(T2)) * Kind[n];
			}*/
			float Lx=px[uABC[0]-1]-px[0];
			float Ly=py[uABC[1]*uABC[0]-1]-py[0];
			for (int n=0; n<NOS; n++)
			{	
				float rx = px[n]-px[0];
				float ry = py[n]-py[0];
				// float rz = pz[n]-pz[0];
				F1 = rx * TPI / Lx;
				T1 = ry * PI / Ly;
				Sx[n] = 0;
				Sy[n] = 0;
				Sz[n] = 1;
				Sx[n] = sin(T1)*cos(F1)*Kind[n];
				Sy[n] = sin(T1)*sin(F1)*Kind[n];
				Sz[n] = cos(T1)*Kind[n];
			}			
		}
		printf("spin0x=%f",px[0]);
	break;
	
	case 15: // Normalize all spins
	double inv_abs_S;
	for (int n=0; n<NOS; n++)
			{	
				if (Kind[n]!=0){
					inv_abs_S=1.0f/sqrt(Sx[n]*Sx[n]+Sy[n]*Sy[n]+Sz[n]*Sz[n]);
					Sx[n] = Sx[n]*inv_abs_S;
					Sy[n] = Sy[n]*inv_abs_S;
					Sz[n] = Sz[n]*inv_abs_S;
				}	
			}
	break;
	}	
}
