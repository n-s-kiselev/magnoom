// ////////////////////////////////////////////// GEOMETRY ////////////////////////////////////////////
// //translation vectors cubic:
// //float		abc[3][3] = {
// // 				{	1.0f, 0.0f, 0.0f }, // a
// // 				{	0.0f, 1.0f, 0.0f }, // b
// // 				{	0.0f, 0.0f, 1.0f }};// c
				
// //atom positions in basic domain: {x0,y0,z0, x1,y1,z1,...}
// //not more then 100 atoms per basic domain
// //FCC
// //float	Block[4*3] = 
// // 		{	0.0f,	0.0f,	0.0f,	
// // 			0.5f,	0.5f,	0.0f,
// // 			0.5f,	0.0f,	0.5f,
// // 			0.0f,	0.5f,	0.5f	
// // 		};
// //BCC
// /*
// float Block[2*3] = 
//     { 0.0f, 0.0f, 0.0f, 
//       0.5f, 0.5f, 0.5f, 
//     };
//  */




// //Simple Cubic 1 [001]
// // float		abc[3][3] = {
// // 				{	1.0f, 0.0f, 0.0f }, // a
// // 				{	0.0f, 1.0f, 0.0f }, // b
// // 				{	0.0f, 0.0f, 1.0f }};// c
// // float			Block[][3] = { 
// // 				{0.5f, 0.5f, 0.5f},  
// // 				};

// /*float		abc[3][3] = {
// 				{	3.0f, 0.0f, 0.0f }, // a
// 				{	0.0f, 1.0f, 0.0f }, // b
// 				{	0.0f, 0.0f, 10.0f }};// c
// float			Block[][3] = { 
// 				{0.0f, 0.5f, 0.0f}, 
// 				{1.5f, 0.5f, 0.0f}, 
// 				{3.0f, 0.5f, 0.0f}, 
// 				{0.75f, 0.0f, 0.0f},
// 				{2.25f, 0.0f, 0.0f},
// 				{0.75f, 1.0f, 0.0f},
// 				{2.25f, 1.0f, 0.0f}  
// 				};*/

// // float		abc[3][3] = {
// // 				{	1.0f, 0.0f, 0.0f }, // a
// // 				{	0.0f, 1.0f, 0.0f }, // b
// // 				{	0.0f, 0.0f, 1.0f }};// c
// // float			Block[][3] = { 
// // 				{0.5f, 0.5f, 0.5f},  
// // 				};
// // Simple Cubic 2 [011]
// // float		abc[3][3] = {
// // 				{	1.0f, 0.0f, 0.0f }, // a
// // 				{	0.0f, sqrt(2.f), 0.0f }, // b
// // 				{	0.0f, 0.0f, sqrt(2.f) }};// c				
// // float		Block[][3] = {		
// // 				{0,		      0,        0},	
// // 				{0.5,	sqrt(2.f)/2,	0}
// // 				};	

// //Simple Cubic 3 [111]
// // float		abc[3][3] = {
// // 				{	sqrt(2.f), 0.0f, 0.0f }, // a
// // 				{	0.0f, sqrt(3.f)*sqrt(2.f), 0.0f }, // b
// // 				{	0.0f, 0.0f, sqrt(3.f) }};// c				
// // float		Block[][3] = {		
// // 				{0,		      0,	        0},	
// // 				{sqrt(2.f)/2, sqrt(6.f)/2,	0},
// // 				{sqrt(2.f)/2, sqrt(6.f)/6,  sqrt(3.f)/3.f},
// // 				{0, 2*sqrt(6.f)/3,  sqrt(3.f)/3.f},
// // 				{0, 1*sqrt(6.f)/3, -sqrt(3.f)/3.f},
// // 				{sqrt(2.f)/2, 5*sqrt(6.f)/6, -sqrt(3.f)/3.f}
// // 				};

// //FCC 2
// // float		abc[3][3] = {
// // 				{	1.0f, 0.0f, 0.0f }, // a
// // 				{	0.0f, sqrt(2.f), 0.0f }, // b
// // 				{	0.0f, 0.0f, sqrt(2.f) }};// c				
// // float		Block[][3] = {		
// // 				{0,		      0,	          0},	
// // 				{0,	sqrt(2.f)/2,			  0},
// // 				{0.5,	sqrt(2.f)/4,sqrt(2.f)/4},
// // 				{0.5, 3*sqrt(2.f)/4, sqrt(2.f)/4},
// // 				{0,	sqrt(2.f)/2,	sqrt(2.f)/2}
// // 				};	
// //FCC 3
// // float		abc[3][3] = {
// // 				{sqrt(2.f), 0.0f, 0.0f }, // a
// // 				{sqrt(2.f)/2, sqrt(6.f)/2, 0.0f }, // b
// // 				{sqrt(2.f)/2, sqrt(6.f)/6, sqrt(3.f)/sqrt(2.f) }};// c				
// // float		Block[][3] = {		
// // 				{0,		      0,	          0},	
// // 				};


// //B20(1)
// // float		abc[3][3] = {
// // 				{	1.0f, 0.0f, 0.0f }, // a
// // 				{	0.0f, 1.0f, 0.0f }, // b
// // 				{	0.0f, 0.0f, 1.0f }};// c
// // #define uB20		0.138//0.138f//0.133
// // float		Block[][3] = {		
// // 				{   uB20,		   uB20,	    uB20},	
// // 				{0.5f+uB20,	0.5f-uB20,	1.0f-uB20},
// // 				{1.0f-uB20,	0.5f+uB20,	0.5f-uB20},	
// // 				{0.5f-uB20,	1.0f-uB20,	0.5f+uB20}	
// // 				};
// //B20(2)

// float		abc[3][3] = {
// 				{	1.0f, 0.0f, 0.0f }, // a
// 				{	0.0f, 1.0f, 0.0f }, // b
// 				{	0.0f, 0.0f, 1.0f }};// c

// float		Block[][3] = {		
// 				{   0.5,		   0.5,	    0.5}	
// 				};


// // #define uB20		0.138f
// // float		Block[][3] = {		
// // 				{   0.,		   0.,	    0.},	
// // 				{0.5f,	0.5f-2*uB20,	1.0f-2*uB20},
// // 				{1.0f-2*uB20,	0.5f,	0.5f-2*uB20},	
// // 				{0.5f-2*uB20,	1.0f-2*uB20,	0.5f}	
// // 				};

// //B20(Maria)
// // float		abc[3][3] = {
// // 				{	1.0f, 0.0f, 0.0f }, // a
// // 				{	0.0f, 1.0f, 0.0f }, // b
// // 				{	0.0f, 0.0f, 1.0f }};// c
// // //#define uB20		0.138f//MnSi
// // // float		Block[][3] = {		
// // // 				{       0.0f,        0.0f,        0.0f},//r1= u, u, u	
// // // 				{       0.5f, 0.5f-2*uB20,     -2*uB20},//r2= 0.5+u, 0.5-u,-u
// // // 				{    -2*uB20,        0.5f, 0.5f-2*uB20},//r3= -u, 0.5+u, 0.5 -u	
// // // 				{0.5f-2*uB20,     -2*uB20,        0.5f} //r4=0.5-u, -u, 0.5+u
// // // 				};
// // #define uB20		0.135f//FeGe
// // float		Block[][3] = {		
// // 				{       0.0f,        0.0f,        0.0f},//r1= u, u, u	
// // 				{       0.5f, 0.5f-2*uB20,     -2*uB20},//r2= 0.5+u, 0.5-u,-u
// // 				{    -2*uB20,        0.5f, 0.5f-2*uB20},//r3= -u, 0.5+u, 0.5 -u	
// // 				{0.5f-2*uB20,     -2*uB20,        0.5f} //r4=0.5-u, -u, 0.5+u
// // 				};

// //number of translations for the basic domain along a,b, and c verctors respectively 
// //int			uABC[3] = {2,147,2};//Grid dimensionality along translation vectors a, b, c; uABC[i]>0 
// int			uABC[3] = {10,10,10};//Grid dimensionality along translation vectors a, b, c; uABC[i]>0 
// //int			uABC[3] = {6,6,6};//Grid dimensionality along translation vectors a, b, c; uABC[i]>0 
// //int			uABC[3] = {80,80,20};//Grid dimensionality along translation vectors a, b, c; uABC[i]>0 


// int			ShellNumber = 1;
// int			AtomsPerBlock = sizeof(Block)/sizeof(float)/3;
// float*		RadiusOfShell = (float *)calloc(ShellNumber , sizeof(float));  
// int*		NeighborsPerAtom = (int *)calloc(AtomsPerBlock, sizeof(int));
// // total number of neighbour pairs per whole map of neighbours
// int			NOS=AtomsPerBlock*uABC[0]*uABC[1]*uABC[2]; // number of spins
// int			NOS_AL=AtomsPerBlock*uABC[1]*uABC[2]; // number of spins per A layer
// int			NOS_BL=AtomsPerBlock*uABC[0]*uABC[2]; // number of spins per B layer
// int			NOS_CL=AtomsPerBlock*uABC[0]*uABC[1]; // number of spins per C layer

// int 		NOSK = 0;

// double 		iNOS = 1.0/NOS;

// int			NOB=uABC[0]*uABC[1]*uABC[2]; // number of Blocks
// int			NOB_AL=uABC[1]*uABC[2]; // number of spins per A layer
// int			NOB_BL=uABC[0]*uABC[2]; // number of spins per B layer
// int			NOB_CL=uABC[0]*uABC[1]; // number of spins per C layer

void	
GetShells(float abc[][3], float block[][3], int atomsPerBlock, int ShellNum, float * R)
{
	//atomsPerBlock is a number of atoms in the basic domain
	float min_distance, test_distance, curr_Radius, Radius, test_Radius;
	float x0,x1,y0,y1,z0,z1,dx,dy,dz;
	int Jin, Kin, Lin;
	//int I0min, I1min, Jmin, Kmin, Lmin;// for each of shell we are looking for the lits of Imin, Jmin, Kmin, Lmin
	if (abc[0][0]==0 && (abc[0][1]==0 && abc[0][2]==0)) {Jin = 0;} else {Jin = 15;} //|a|=0
	if (abc[1][0]==0 && (abc[1][1]==0 && abc[1][2]==0)) {Kin = 0;} else {Kin = 15;} //|b|=0
	if (abc[2][0]==0 && (abc[2][1]==0 && abc[2][2]==0)) {Lin = 0;} else {Lin = 15;} //|c|=0
		Radius = 0.0f;// first we are looking for the minimal interatomic distance closest to 0
		for (int shell = 0; shell < ShellNum; ++shell) // runs ove shells
		{	// at each 'shell'-step we use previouse shell radius starting from R=0.0f
			curr_Radius = Radius;	 
			min_distance = 100.f; // initial value for the distance between current and the tested shell
			for(int I0=0; I0 < atomsPerBlock; I0++) // runs over atoms in the block 
			{	//position of the atom with index Idx0 in the block:
				x0 = block[I0][0];
				y0 = block[I0][1];
				z0 = block[I0][2];
				for(int J=-Jin;J<=Jin;J++) // translation of block along vector 'a' J times
				{
					for(int K=-Kin;K<=Kin;K++)// translation of block along vector 'b' K times
					{
						for( int L=-Lin;L<=Lin;L++)// translation of block along vector 'c' L times
						{	
							for(int I1=0; I1<atomsPerBlock; I1++) // runs over atoms in one of neghbouring block
							{//tested atom position:
								x1 = block[I1][0] + abc[0][0]*J + abc[1][0]*K + abc[2][0]*L; 
								dx = x0-x1; //printf("dx = %f\n", dx );
								y1 = block[I1][1] + abc[0][1]*J + abc[1][1]*K + abc[2][1]*L; 
								dy = y0-y1; //printf("dy = %f\n", dy);
								z1 = block[I1][2] + abc[0][2]*J + abc[1][2]*K + abc[2][2]*L; 
								dz = z0-z1; //printf("dz = %f\n", dz );
								if ( !( (I0==I1 && 0==J) && (0==K && 0==L) ) ) // if it's not the same atom
								{	
									test_Radius = sqrt(dx*dx + dy*dy + dz*dz);
									test_distance = test_Radius - curr_Radius;// distance between shells
									//printf("test_Radius = %f\n", test_Radius );
									//printf("test_distance = %f\n", test_distance );
									if (test_distance > 1e-6 && test_distance < min_distance)
									{
										// if(0==J*K && 0==J*L && 0==L*K)//metka only for for micromagnetics
										{
											min_distance = test_distance;
											//I0min = I0; I1min = I1; Jmin = J; Kmin = K; Lmin = L;
											R[ shell ] = Radius = test_Radius;
										}

									}
								} 
							}//Idx1	
						}//L
					}//K
				}//J
			}//Idx0
			//printf("R[%d] =%f, IdxMin =%d, Imin =%d, Jmin =%d, Kmin = %d, Lmin =%d \n", shell, R[ shell ], IdxMin,Imin,Jmin,Kmin,Lmin);
		}//sIdx
}

int
GetNeighborsNumber(
	float abc[][3], float block[][3], int atomsPerBlock, 
	int ShellNum, float * R, int * NeighborsNum)
{
	//atomsPerBlock is a number of atoms in the basic domain
	float delta, test_Radius, Radius;
	float x0,x1,y0,y1,z0,z1,dx,dy,dz;
	int NumInShell, TotNum=0;
	int Jin, Kin, Lin;
	if (abc[0][0]==0 && (abc[0][1]==0 && abc[0][2]==0)) {Jin = 0;} else {Jin = 15;}
	if (abc[1][0]==0 && (abc[1][1]==0 && abc[1][2]==0)) {Kin = 0;} else {Kin = 15;}
	if (abc[2][0]==0 && (abc[2][1]==0 && abc[2][2]==0)) {Lin = 0;} else {Lin = 15;}
	for(int I0=0; I0 < atomsPerBlock; I0++) // runs over atoms in the block 
		{	//position of the atom with index Idx0 in the block:
			x0 = block[I0][0];
			y0 = block[I0][1];
			z0 = block[I0][2];
			NeighborsNum[I0] = 0; 
			for (int shell = 0; shell < ShellNum; ++shell) // runs ove shells 
			{
				Radius = R [shell];
				NumInShell=0;
				for(int J=-Jin;J<=Jin;J++) // translation of basic domain along vector 'a' J times
				{
					for(int K=-Kin;K<=Kin;K++)// translation of basic domain along vector 'b' K times
					{ 
						for( int L=-Lin;L<=Lin;L++)// translation of basic domain along vector 'c' L times
						{	
							for(int I1=0; I1<atomsPerBlock; I1++) // runs over atoms in the shifted basic domain
							{//tested atom position:
								x1 = block[I1][0] + abc[0][0]*J + abc[1][0]*K + abc[2][0]*L; 
								dx = x0-x1; //printf("dx = %f\n", dx );
								y1 = block[I1][1] + abc[0][1]*J + abc[1][1]*K + abc[2][1]*L; 
								dy = y0-y1; //printf("dy = %f\n", dy);
								z1 = block[I1][2] + abc[0][2]*J + abc[1][2]*K + abc[2][2]*L; 
								dz = z0-z1; //printf("dz = %f\n", dz );
								test_Radius = sqrt(dx*dx + dy*dy + dz*dz); 
								delta = fabs (Radius - test_Radius);
								if (delta<1e-6) //atom Idx1 is in shell 
								{
									NumInShell += 1;
									NeighborsNum[I0] +=1; 
								}
							}//Idx0
						}//L
					}//K
				}//J
			//printf("NumInShell[%d] =%d, NeighborsNum[%d]=%d \n", shell, NumInShell, Idx0,NeighborsNum[Idx0]);	
			}//sell
		TotNum+=NeighborsNum[I0];
		}//sIdx
		return TotNum;
}


void
CreateMapOfNeighbors( 
	float abc[][3], float block[][3], int atomsPerBlock, int shellNum, float * R,
	int* aiBlock, int* niBlock, int* niGridA, int* niGridB, int* niGridC, int* sIdx)
{
	float delta, test_Radius, Radius;
	float x0,x1,y0,y1,z0,z1,dx,dy,dz;
	int Jin, Kin, Lin, N=-1;
	if (abc[0][0]==0 && (abc[0][1]==0 && abc[0][2]==0)) {Jin = 0;} else {Jin = 5;}
	if (abc[1][0]==0 && (abc[1][1]==0 && abc[1][2]==0)) {Kin = 0;} else {Kin = 5;}
	if (abc[2][0]==0 && (abc[2][1]==0 && abc[2][2]==0)) {Lin = 0;} else {Lin = 5;}
	for(int I0=0; I0 < atomsPerBlock; I0++) // runs over atoms in basic domain 
	{	//position of the atom with index Idx0 in the basic domain:
		x0 = block[I0][0];
		y0 = block[I0][1];
		z0 = block[I0][2];
		for (int shell = 0; shell < shellNum; ++shell) // runs ove shells 
		{
			Radius = R [shell];
			for(int J=-Jin;J<=Jin;J++) // translation of basic domain along vector 'a' J times
			{
				for(int K=-Kin;K<=Kin;K++)// translation of basic domain along vector 'b' K times
				{
					for( int L=-Lin;L<=Lin;L++)// translation of basic domain along vector 'c' L times
					{	
						for(int I1=0; I1<atomsPerBlock; I1++) // runs over atoms in the shifted basic domain
						{//tested atom position:
							x1 = block[I1][0] + abc[0][0]*J + abc[1][0]*K + abc[2][0]*L; 
							dx = x0-x1; 
							y1 = block[I1][1] + abc[0][1]*J + abc[1][1]*K + abc[2][1]*L; 
							dy = y0-y1; 
							z1 = block[I1][2] + abc[0][2]*J + abc[1][2]*K + abc[2][2]*L; 
							dz = z0-z1; 
							test_Radius = sqrt(dx*dx + dy*dy + dz*dz); 
							delta = fabs (Radius - test_Radius);
							if (delta<1e-6) //atom I1 is in shell 
							{
								N++; 
								aiBlock[N] = I0;
								niBlock[N] = I1;
								niGridA[N] = J;
								niGridB[N] = K;
								niGridC[N] = L;
								sIdx[N] = shell;
								//printf("I'[%2d] = %2d", N, *aiBlock[N]);
								// printf("I [%2d] = %2d", N, NIdxBlock[N]);
								// printf("J [%2d] = %2d", N, NIdxGridA[N]);
								// printf("K [%2d] = %2d", N, NIdxGridB[N]);
								// printf("L [%2d] = %2d\n", N, NIdxGridC[N]);
							}
						}//Idx0
					}//L
				}//K
			}//J
		//printf("NumInShell[%d] =%d, NeighborsNum[%d]=%d \n", shell, NumInShell, Idx0,NeighborsNum[Idx0]);	
		}//sell
	}//sIdx
}

void
GetPosition( float abc[][3], float block[][3], int I, int J, int K, int L, float * XYZ)
{
	XYZ[0]=block[I][0]+abc[0][0]*J+abc[1][0]*K+abc[2][0]*L;
	XYZ[1]=block[I][1]+abc[0][1]*J+abc[1][1]*K+abc[2][1]*L;
	XYZ[2]=block[I][2]+abc[0][2]*J+abc[1][2]*K+abc[2][2]*L;
}

void
SetExch1( float abc[][3], float block[][3], int arrSize, 
	float* Jijs, //input-> exchange coupling constant in each shell
	float* Bijs, //input-> bi-quadratic exchange coupling constant in each shell
	float* Dijs, //input-> modul of DMI vecotor in each shell
	int* aiBlock, int* niBlock, int* niGridA, int* niGridB, int* niGridC, int* sIdx, 
	float* J0, float* B0, float* D0, float* Dx, float* Dy, float* Dz)
{
	int I0,I1,J,K,L,S;
	double XYZ0[3], XYZ1[3];

    for(int i=0; i<arrSize; i++)
    {
    	I0= aiBlock[i];
    	I1= niBlock[i];
    	J = niGridA[i];
    	K = niGridB[i];
    	L = niGridC[i];
    	S = sIdx[i]; //shell index
    	//spin i position
		XYZ0[0]=block[I0][0];
		XYZ0[1]=block[I0][1];
		XYZ0[2]=block[I0][2];
		//spin j position
		XYZ1[0]=block[I1][0]+abc[0][0]*J+abc[1][0]*K+abc[2][0]*L;
		XYZ1[1]=block[I1][1]+abc[0][1]*J+abc[1][1]*K+abc[2][1]*L;
		XYZ1[2]=block[I1][2]+abc[0][2]*J+abc[1][2]*K+abc[2][2]*L;
		//r_{ij}
		XYZ1[0]-=XYZ0[0];//r1=r1-r0 /x
		XYZ1[1]-=XYZ0[1];//r1=r1-r0 /y
		XYZ1[2]-=XYZ0[2];//r1=r1-r0 /z	
		//metka
		// double norm[3]={0,0,1};
		// if (S==0) {
		// 	Cross(norm, XYZ1, XYZ0); (void)Unit( XYZ0, XYZ0);//interface induced DMI -- skyrmions of Neel type.
		// } else if (S==1){
		// 	(void)Unit( XYZ1, XYZ0);//Bloch skyrmion
		// }
		(void)Unit( XYZ1, XYZ0);//Bloch skyrmion
		Dx[i] = XYZ0[0];
		Dy[i] = XYZ0[1];
		Dz[i] = XYZ0[2];

		J0[i] = Jijs[S];
		D0[i] = Dijs[S];
		B0[i] = Bijs[S];
    }
}

