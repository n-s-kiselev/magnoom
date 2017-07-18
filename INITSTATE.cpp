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
		if (r<Sk_R)
		{
			T=PI*exp(-2*r/Sk_R);//<-- defines skyrmion profile you may put periodical function to get target like skyrmions
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
	if(chSize!=0)
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
		{	
			rnd[0] = 2.0 * (0.5 - rand() / (double)RAND_MAX);
			rnd[1] = 2.0 * (0.5 - rand() / (double)RAND_MAX);
			rnd[2] = 2.0 * (0.5 - rand() / (double)RAND_MAX);
			(void)Unit(rnd,rnd);
			Sx[n] = rnd[0]*Kind[n];
			Sy[n] = rnd[1]*Kind[n];
			Sz[n] = rnd[2]*Kind[n];	
		}	
	break;

	case 1: //homogeneous
		for (int n=0; n<NOS; n++)
		{	
			Sx[n] = chDir[0]*Kind[n];
			Sy[n] = chDir[1]*Kind[n];
			Sz[n] = chDir[2]*Kind[n];	
		}	
	break;

	case 2: // skyrmion Q=1
		if(chSize!=0)
		{
			CreatSkyrmion(px, py, pz, sx, sy, sz, chSize,0,0);	
		}
	break;

	case 3: // skyrmion Q=2
		if(chSize!=0)
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
			F=-2*(-atan2(px[n],py[n]))+PI*0.5;//<-- achiral skyrmion |Q|=2
			//F=-3*(-atan2(px[n],py[n]))+PI*0.5;//<-- achiral skyrmion |Q|=3
			Sx[n] = sin(T)*cos(F)*Kind[n];
			Sy[n] = sin(T)*sin(F)*Kind[n];
			Sz[n] = cos(T)*Kind[n];	
			}	
		}
	break;

	case 4: // skyrmion Q=3
		if(chSize!=0)
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
		if(chSize!=0)
		
		{
			CreatBobber(px, py, pz, sx, sy, sz, chSize, 0, 0, 1);	
		}
	break;

	case 6: // bobber_bottom
		if(chSize!=0)
		
		{
			CreatBobber(px, py, pz, sx, sy, sz, chSize, 0, 0, -1);				
		}
	break;

	case 7: // bobber_lattice
		if(chSize!=0 && !(chDir[0]==0&&chDir[1]==0) )
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
		if(chSize!=0 && !(chDir[0]==0&&chDir[1]==0) )
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
		if(chSize!=0 && !(chDir[0]==0&&chDir[1]==0) )
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
		if(chSize!=0)
		{
			float tmp;
			for (int n=0; n<NOS; n++)
			{	
				r = sqrt(px[n]*px[n]+py[n]*py[n]+pz[n]*pz[n]);
				if (r==0){
					T = 0;
				}else{
					T = pz[n]/r; // angle with respect to the main axis of toroid [0,0,1]
				}
				T = acos(T);
				t = r/chSize;
				t = 1.0 + 4.22/(t*t);
				tmp = PI*(1.0-1.0/sqrt(t));
				t = sin(tmp)*sin(T);
				t = acos(1.0-2.0*t*t);
				F = atan2(py[n],px[n]);
				f = F + atan2( 1.0/(tan(tmp)),cos(T) );
				Sx[n] = sin(t)*cos(f)*Kind[n];
				Sy[n] = sin(t)*sin(f)*Kind[n];
				Sz[n] = cos(t)*Kind[n];
			}
		}
	// CreatBobber(px, py, pz, sx, sy, sz, chSize, -20, 0, 1);
	// CreatSkyrmion(px, py, pz, sx, sy, sz, chSize,20,0);
	break;

	case 11: // spiral
		if(fabs(chDir[0])+fabs(chDir[1])+fabs(chDir[2])!=0)
		{		
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
		}
	break;

	case 12: // skyrmion Lattice
		if(chSize!=0 && !(chDir[0]==0&&chDir[1]==0) )
		{

			float PerSkirmLat=chSize;
			float aSkirmLat=2*PerSkirmLat/(sqrt(3)); 
			float Sk_Radius=0.8*0.5*aSkirmLat; 
			float TranLatX;
			float TranLatY;
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
					TranLatX=aSkirmLat*kx+aSkirmLat*0.5*ky;//Tranlation of skirmion Lattice (along X) without rotate
					TranLatY=PerSkirmLat*ky;//Tranlation of skirmion Lattice (along Y) without rotate
					RotTranLatX=CoChDir*TranLatX-SiChDir*TranLatY;//Tranlation of skirmion Lattice (along Y) with rotate in XY
					RotTranLatY=SiChDir*TranLatX+CoChDir*TranLatY;//Tranlation of skirmion Lattice (along Y) with rotate in XY
					CreatSkyrmion(px, py, pz, sx, sy, sz, Sk_Radius, RotTranLatX, RotTranLatY);
				}
			}	
		}
	break;

	case 13: // globula
		if(chSize!=0)
		{
			float T = 0.f;
			float F = 0.f;
			float r1,r2;
			float rx=0.0f;
			float ry=0.0f;
			float rz=0.0f;

			for (int n=0; n<NOS; n++)
			{
				rx=px[n];
				ry=py[n];
				rz=pz[n];
				r1 = sqrt(rx*rx+ry*ry);
				r2 = sqrt(rx*rx+ry*ry+rz*rz);
				if (r2<chSize*0.5)
				{
					T=PI*exp(-2*r1/chSize);//<-- defines skyrmion profile you may put periodical function to get target like skyrmions
					F= atan2(ry,rx)+PI*0.5;//<--chiral (bloch) skyrmion |Q|=1
					Sx[n] = sin(T)*cos(F)*Kind[n];
					Sy[n] = sin(T)*sin(F)*Kind[n];
					Sz[n] = cos(T)*Kind[n];			
				}else{
					// Sx[n] = 0;
					// Sy[n] = 0;
					// Sz[n] = 1;						
				}
			}
		}
	break;

	case 14: // 3D lattice of Bloch points
		if(chSize!=0)
		{
			float T1 = 0.f;
			float T2 = 0.f;
			float F1 = 0.f;
			float F2 = 0.f;
			float rx=0.0f;
			float ry=0.0f;
			float rz=0.0f;

			for (int n=0; n<NOS; n++)
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
			}
		}
	break;
	
	case 15: // Normalize all spins
	float inv_abs_S;
	for (int n=0; n<NOS; n++)
			{	
				inv_abs_S=1.0f/sqrt(Sx[n]*Sx[n]+Sy[n]*Sy[n]+Sz[n]*Sz[n]);
				Sx[n] = Sx[n]*inv_abs_S;
				Sy[n] = Sy[n]*inv_abs_S;
				Sz[n] = Sz[n]*inv_abs_S;
			}
	break;
	}	
}