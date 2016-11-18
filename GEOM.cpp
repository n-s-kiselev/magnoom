////////////////////////////////////////////// GEOMETRY ////////////////////////////////////////////
//translation vectors cubic:
//float		abc[3][3] = {
// 				{	1.0f, 0.0f, 0.0f }, // a
// 				{	0.0f, 1.0f, 0.0f }, // b
// 				{	0.0f, 0.0f, 1.0f }};// c
				
//atom positions in basic domain: {x0,y0,z0, x1,y1,z1,...}
//not more then 100 atoms per basic domain
//FCC
//float	Block[4*3] = 
// 		{	0.0f,	0.0f,	0.0f,	
// 			0.5f,	0.5f,	0.0f,
// 			0.5f,	0.0f,	0.5f,
// 			0.0f,	0.5f,	0.5f	
// 		};
//BCC
/*
float Block[2*3] = 
    { 0.0f, 0.0f, 0.0f, 
      0.5f, 0.5f, 0.5f, 
    };
 */
//Simple Cubic 1 [001]
float		abc[3][3] = {
				{	1.0f, 0.0f, 0.0f }, // a
				{	0.0f, 1.0f, 0.0f }, // b
				{	0.0f, 0.0f, 1.0f }};// c
float			Block[][3] = { 
				{0.5f, 0.5f, 0.5f},  
				};
// Simple Cubic 2 [011]
// float		abc[3][3] = {
// 				{	1.0f, 0.0f, 0.0f }, // a
// 				{	0.0f, sqrt(2.f), 0.0f }, // b
// 				{	0.0f, 0.0f, sqrt(2.f) }};// c				
// float		Block[][3] = {		
// 				{0,		      0,        0},	
// 				{0.5,	sqrt(2.f)/2,	0}
// 				};	

//Simple Cubic 3 [111]
// float		abc[3][3] = {
// 				{	1.0f, 0.0f, 0.0f }, // a
// 				{	0.0f, sqrt(2.f), 0.0f }, // b
// 				{	0.0f, 0.0f, sqrt(2.f) }};// c				
// float		Block[][3] = {		
// 				{0,		      0,	          0},	
// 				{0,	sqrt(2.f)/2,	sqrt(2.f)/2}
// 				};

//FCC 2
// float		abc[3][3] = {
// 				{	1.0f, 0.0f, 0.0f }, // a
// 				{	0.0f, sqrt(2.f), 0.0f }, // b
// 				{	0.0f, 0.0f, sqrt(2.f) }};// c				
// float		Block[][3] = {		
// 				{0,		      0,	          0},	
// 				{0,	sqrt(2.f)/2,			  0},
// 				{0.5,	sqrt(2.f)/4,sqrt(2.f)/4},
// 				{0.5, 3*sqrt(2.f)/4, sqrt(2.f)/4},
// 				{0,	sqrt(2.f)/2,	sqrt(2.f)/2}
// 				};	
//FCC 3
// float		abc[3][3] = {
// 				{sqrt(2.f), 0.0f, 0.0f }, // a
// 				{sqrt(2.f)/2, sqrt(6.f)/2, 0.0f }, // b
// 				{sqrt(2.f)/2, sqrt(6.f)/6, sqrt(3.f)/sqrt(2.f) }};// c				
// float		Block[][3] = {		
// 				{0,		      0,	          0},	
// 				};


//B20(1)
// #define uB20		0.138f
// float		Block[][3] = {		
// 				{   uB20,		   uB20,	    uB20},	
// 				{0.5f+uB20,	0.5f-uB20,	1.0f-uB20},
// 				{1.0f-uB20,	0.5f+uB20,	0.5f-uB20},	
// 				{0.5f-uB20,	1.0f-uB20,	0.5f+uB20}	
// 				};
//B20(2)
/*
float		abc[3][3] = {
				{	1.0f, 0.0f, 0.0f }, // a
				{	0.0f, 1.0f, 0.0f }, // b
				{	0.0f, 0.0f, 1.0f }};// c
#define uB20		0.138f
float		Block[][3] = {		
				{   0.,		   0.,	    0.},	
				{0.5f,	0.5f-2*uB20,	1.0f-2*uB20},
				{1.0f-2*uB20,	0.5f,	0.5f-2*uB20},	
				{0.5f-2*uB20,	1.0f-2*uB20,	0.5f}	
				};
*/
//B20(Maria)
// float		abc[3][3] = {
// 				{	1.0f, 0.0f, 0.0f }, // a
// 				{	0.0f, 1.0f, 0.0f }, // b
// 				{	0.0f, 0.0f, 1.0f }};// c
// //#define uB20		0.138f//MnSi
// // float		Block[][3] = {		
// // 				{       0.0f,        0.0f,        0.0f},//r1= u, u, u	
// // 				{       0.5f, 0.5f-2*uB20,     -2*uB20},//r2= 0.5+u, 0.5-u,-u
// // 				{    -2*uB20,        0.5f, 0.5f-2*uB20},//r3= -u, 0.5+u, 0.5 -u	
// // 				{0.5f-2*uB20,     -2*uB20,        0.5f} //r4=0.5-u, -u, 0.5+u
// // 				};
// #define uB20		0.135f//FeGe
// float		Block[][3] = {		
// 				{       0.0f,        0.0f,        0.0f},//r1= u, u, u	
// 				{       0.5f, 0.5f-2*uB20,     -2*uB20},//r2= 0.5+u, 0.5-u,-u
// 				{    -2*uB20,        0.5f, 0.5f-2*uB20},//r3= -u, 0.5+u, 0.5 -u	
// 				{0.5f-2*uB20,     -2*uB20,        0.5f} //r4=0.5-u, -u, 0.5+u
// 				};

//number of translations for the basic domain along a,b, and c verctors respectively 
//int			ABC[3] = {2,147,2};//Grid dimensionality along translation vectors a, b, c; ABC[i]>0 
//int			ABC[3] = {71,41,1};//Grid dimensionality along translation vectors a, b, c; ABC[i]>0 
int			ABC[3] = {100,100,100};//Grid dimensionality along translation vectors a, b, c; ABC[i]>0 
int			Boundary[3] = {0, 0, 0};// boundary conditions along a, b, c translation vectors

int			ShellNumber = 1;
int			AtomsPerBlock = sizeof(Block)/sizeof(float)/3;
float*		RadiusOfShell = (float *)calloc(ShellNumber , sizeof(float));  
int*		NeighborsPerAtom = (int *)calloc(AtomsPerBlock, sizeof(int));
// total number of neighbour pairs per whole map of neighbours
int			NOS=AtomsPerBlock*ABC[0]*ABC[1]*ABC[2]; // number of spins
int			NOS_AL=AtomsPerBlock*ABC[1]*ABC[2]; // number of spins per A layer
int			NOS_BL=AtomsPerBlock*ABC[0]*ABC[2]; // number of spins per B layer
int			NOS_CL=AtomsPerBlock*ABC[0]*ABC[1]; // number of spins per C layer

double 		iNOS = 1.0/NOS;

int			NOB=ABC[0]*ABC[1]*ABC[2]; // number of Blocks
int			NOB_AL=ABC[1]*ABC[2]; // number of spins per A layer
int			NOB_BL=ABC[0]*ABC[2]; // number of spins per B layer
int			NOB_CL=ABC[0]*ABC[1]; // number of spins per C layer

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
			min_distance = 100.f; // initial value for distance between current and tested shell
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
										min_distance = test_distance;
										//I0min = I0; I1min = I1; Jmin = J; Kmin = K; Lmin = L;
										R[ shell ] = Radius = test_Radius;
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
	double XYZ0[3], XYZ1[3], norm[3];
	norm[0] = 0.0f;
	norm[1] = 0.0f;
	norm[2] = 1.0f;

    for(int i=0; i<arrSize; i++)
    {
    	I0= aiBlock[i];
    	I1= niBlock[i];
    	J = niGridA[i];
    	K = niGridB[i];
    	L = niGridC[i];
    	S = sIdx[i]; //shell index

		XYZ0[0]=block[I0][0];
		XYZ0[1]=block[I0][1];
		XYZ0[2]=block[I0][2];

		XYZ1[0]=block[I1][0]+abc[0][0]*J+abc[1][0]*K+abc[2][0]*L;
		XYZ1[1]=block[I1][1]+abc[0][1]*J+abc[1][1]*K+abc[2][1]*L;
		XYZ1[2]=block[I1][2]+abc[0][2]*J+abc[1][2]*K+abc[2][2]*L;

		XYZ1[0]-=XYZ0[0];//r1=r1-r0 /x
		XYZ1[1]-=XYZ0[1];//r1=r1-r0 /y
		XYZ1[2]-=XYZ0[2];//r1=r1-r0 /z
		//Cross(norm, XYZ1, XYZ0); (void)Unit( XYZ0, XYZ0);//interface induced DMI -- skyrmions of Neel type.

		(void)Unit( XYZ1, XYZ0);//Bloch skyrmion
		Dx[i] = XYZ0[0];
		Dy[i] = XYZ0[1];
		Dz[i] = XYZ0[2];
		J0[i] = Jijs[S];
		D0[i] = Dijs[S];
		B0[i] = Bijs[S];
    }
}

void
SetExchMariya( float abc[][3], float block[][3], int arrSize, 
	float* Jijs, //input-> exchange coupling constant in each shell
	float* Bijs, //input-> bi-quadratic exchange coupling constant in each shell
	float* Dijs, //input-> modul of DMI vecotor in each shell
	int* aiBlock, int* niBlock, int* niGridA, int* niGridB, int* niGridC, int* sIdx, 
	float* J0, float* B0, float* D0, float* Dx, float* Dy, float* Dz)
{
	Dx[	0	]=	-0.328	;		Dy[	0	]=	-0.081	;		Dz[	0	]=	0.457	;	J0[	0	]=	15.396	;
	Dx[	1	]=	0.457	;		Dy[	1	]=	-0.328	;		Dz[	1	]=	-0.081	;	J0[	1	]=	15.396	;
	Dx[	2	]=	-0.081	;		Dy[	2	]=	0.457	;		Dz[	2	]=	-0.328	;	J0[	2	]=	15.396	;
	Dx[	3	]=	0.328	;		Dy[	3	]=	-0.081	;		Dz[	3	]=	0.457	;	J0[	3	]=	15.396	;
	Dx[	4	]=	0.457	;		Dy[	4	]=	0.328	;		Dz[	4	]=	-0.081	;	J0[	4	]=	15.396	;
	Dx[	5	]=	-0.081	;		Dy[	5	]=	0.457	;		Dz[	5	]=	0.328	;	J0[	5	]=	15.396	;
	Dx[	6	]=	0.146	;		Dy[	6	]=	-0.009	;		Dz[	6	]=	0.089	;	J0[	6	]=	-1.715	;
	Dx[	7	]=	-0.146	;		Dy[	7	]=	-0.009	;		Dz[	7	]=	0.089	;	J0[	7	]=	-1.715	;
	Dx[	8	]=	-0.009	;		Dy[	8	]=	0.089	;		Dz[	8	]=	0.146	;	J0[	8	]=	-1.715	;
	Dx[	9	]=	-0.009	;		Dy[	9	]=	0.089	;		Dz[	9	]=	-0.146	;	J0[	9	]=	-1.715	;
	Dx[	10	]=	0.089	;		Dy[	10	]=	0.146	;		Dz[	10	]=	-0.009	;	J0[	10	]=	-1.715	;
	Dx[	11	]=	0.089	;		Dy[	11	]=	-0.146	;		Dz[	11	]=	-0.009	;	J0[	11	]=	-1.715	;
	Dx[	12	]=	0.212	;		Dy[	12	]=	-0.064	;		Dz[	12	]=	-0.005	;	J0[	12	]=	-2.592	;
	Dx[	13	]=	-0.064	;		Dy[	13	]=	-0.005	;		Dz[	13	]=	0.212	;	J0[	13	]=	-2.592	;
	Dx[	14	]=	-0.064	;		Dy[	14	]=	-0.005	;		Dz[	14	]=	-0.212	;	J0[	14	]=	-2.592	;
	Dx[	15	]=	-0.005	;		Dy[	15	]=	0.212	;		Dz[	15	]=	-0.064	;	J0[	15	]=	-2.592	;
	Dx[	16	]=	-0.212	;		Dy[	16	]=	-0.064	;		Dz[	16	]=	-0.005	;	J0[	16	]=	-2.592	;
	Dx[	17	]=	-0.005	;		Dy[	17	]=	-0.212	;		Dz[	17	]=	-0.064	;	J0[	17	]=	-2.592	;
	Dx[	18	]=	-0.034	;		Dy[	18	]=	0.007	;		Dz[	18	]=	-0.059	;	J0[	18	]=	2.867	;
	Dx[	19	]=	-0.059	;		Dy[	19	]=	-0.034	;		Dz[	19	]=	0.007	;	J0[	19	]=	2.867	;
	Dx[	20	]=	0.007	;		Dy[	20	]=	-0.059	;		Dz[	20	]=	-0.034	;	J0[	20	]=	2.867	;
	Dx[	21	]=	-0.007	;		Dy[	21	]=	0.059	;		Dz[	21	]=	0.034	;	J0[	21	]=	2.867	;
	Dx[	22	]=	0.059	;		Dy[	22	]=	0.034	;		Dz[	22	]=	-0.007	;	J0[	22	]=	2.867	;
	Dx[	23	]=	0.034	;		Dy[	23	]=	-0.007	;		Dz[	23	]=	0.059	;	J0[	23	]=	2.867	;
	Dx[	24	]=	0.457	;		Dy[	24	]=	-0.328	;		Dz[	24	]=	0.081	;	J0[	24	]=	15.396	;
	Dx[	25	]=	-0.328	;		Dy[	25	]=	0.081	;		Dz[	25	]=	-0.457	;	J0[	25	]=	15.396	;
	Dx[	26	]=	0.457	;		Dy[	26	]=	0.328	;		Dz[	26	]=	0.081	;	J0[	26	]=	15.396	;
	Dx[	27	]=	-0.081	;		Dy[	27	]=	-0.457	;		Dz[	27	]=	-0.328	;	J0[	27	]=	15.396	;
	Dx[	28	]=	0.328	;		Dy[	28	]=	0.081	;		Dz[	28	]=	-0.457	;	J0[	28	]=	15.396	;
	Dx[	29	]=	-0.081	;		Dy[	29	]=	-0.457	;		Dz[	29	]=	0.328	;	J0[	29	]=	15.396	;
	Dx[	30	]=	0.146	;		Dy[	30	]=	0.009	;		Dz[	30	]=	-0.089	;	J0[	30	]=	-1.715	;
	Dx[	31	]=	-0.009	;		Dy[	31	]=	-0.089	;		Dz[	31	]=	0.146	;	J0[	31	]=	-1.715	;
	Dx[	32	]=	-0.009	;		Dy[	32	]=	-0.089	;		Dz[	32	]=	-0.146	;	J0[	32	]=	-1.715	;
	Dx[	33	]=	-0.146	;		Dy[	33	]=	0.009	;		Dz[	33	]=	-0.089	;	J0[	33	]=	-1.715	;
	Dx[	34	]=	0.089	;		Dy[	34	]=	0.146	;		Dz[	34	]=	0.009	;	J0[	34	]=	-1.715	;
	Dx[	35	]=	0.089	;		Dy[	35	]=	-0.146	;		Dz[	35	]=	0.009	;	J0[	35	]=	-1.715	;
	Dx[	36	]=	-0.064	;		Dy[	36	]=	0.005	;		Dz[	36	]=	0.212	;	J0[	36	]=	-2.592	;
	Dx[	37	]=	-0.064	;		Dy[	37	]=	0.005	;		Dz[	37	]=	-0.212	;	J0[	37	]=	-2.592	;
	Dx[	38	]=	-0.005	;		Dy[	38	]=	0.212	;		Dz[	38	]=	0.064	;	J0[	38	]=	-2.592	;
	Dx[	39	]=	0.212	;		Dy[	39	]=	0.064	;		Dz[	39	]=	0.005	;	J0[	39	]=	-2.592	;
	Dx[	40	]=	-0.005	;		Dy[	40	]=	-0.212	;		Dz[	40	]=	0.064	;	J0[	40	]=	-2.592	;
	Dx[	41	]=	-0.212	;		Dy[	41	]=	0.064	;		Dz[	41	]=	0.005	;	J0[	41	]=	-2.592	;
	Dx[	42	]=	-0.034	;		Dy[	42	]=	-0.007	;		Dz[	42	]=	0.059	;	J0[	42	]=	2.867	;
	Dx[	43	]=	0.059	;		Dy[	43	]=	-0.034	;		Dz[	43	]=	0.007	;	J0[	43	]=	2.867	;
	Dx[	44	]=	-0.007	;		Dy[	44	]=	-0.059	;		Dz[	44	]=	-0.034	;	J0[	44	]=	2.867	;
	Dx[	45	]=	0.007	;		Dy[	45	]=	0.059	;		Dz[	45	]=	0.034	;	J0[	45	]=	2.867	;
	Dx[	46	]=	-0.059	;		Dy[	46	]=	0.034	;		Dz[	46	]=	-0.007	;	J0[	46	]=	2.867	;
	Dx[	47	]=	0.034	;		Dy[	47	]=	0.007	;		Dz[	47	]=	-0.059	;	J0[	47	]=	2.867	;
	Dx[	48	]=	0.081	;		Dy[	48	]=	0.457	;		Dz[	48	]=	-0.328	;	J0[	48	]=	15.396	;
	Dx[	49	]=	0.081	;		Dy[	49	]=	0.457	;		Dz[	49	]=	0.328	;	J0[	49	]=	15.396	;
	Dx[	50	]=	-0.328	;		Dy[	50	]=	-0.081	;		Dz[	50	]=	-0.457	;	J0[	50	]=	15.396	;
	Dx[	51	]=	-0.457	;		Dy[	51	]=	-0.328	;		Dz[	51	]=	0.081	;	J0[	51	]=	15.396	;
	Dx[	52	]=	-0.457	;		Dy[	52	]=	0.328	;		Dz[	52	]=	0.081	;	J0[	52	]=	15.396	;
	Dx[	53	]=	0.328	;		Dy[	53	]=	-0.081	;		Dz[	53	]=	-0.457	;	J0[	53	]=	15.396	;
	Dx[	54	]=	-0.089	;		Dy[	54	]=	0.146	;		Dz[	54	]=	0.009	;	J0[	54	]=	-1.715	;
	Dx[	55	]=	0.146	;		Dy[	55	]=	-0.009	;		Dz[	55	]=	-0.089	;	J0[	55	]=	-1.715	;
	Dx[	56	]=	-0.089	;		Dy[	56	]=	-0.146	;		Dz[	56	]=	0.009	;	J0[	56	]=	-1.715	;
	Dx[	57	]=	0.009	;		Dy[	57	]=	0.089	;		Dz[	57	]=	0.146	;	J0[	57	]=	-1.715	;
	Dx[	58	]=	0.009	;		Dy[	58	]=	0.089	;		Dz[	58	]=	-0.146	;	J0[	58	]=	-1.715	;
	Dx[	59	]=	-0.146	;		Dy[	59	]=	-0.009	;		Dz[	59	]=	-0.089	;	J0[	59	]=	-1.715	;
	Dx[	60	]=	0.212	;		Dy[	60	]=	-0.064	;		Dz[	60	]=	0.005	;	J0[	60	]=	-2.592	;
	Dx[	61	]=	0.064	;		Dy[	61	]=	-0.005	;		Dz[	61	]=	0.212	;	J0[	61	]=	-2.592	;
	Dx[	62	]=	-0.212	;		Dy[	62	]=	-0.064	;		Dz[	62	]=	0.005	;	J0[	62	]=	-2.592	;
	Dx[	63	]=	0.005	;		Dy[	63	]=	0.212	;		Dz[	63	]=	0.064	;	J0[	63	]=	-2.592	;
	Dx[	64	]=	0.064	;		Dy[	64	]=	-0.005	;		Dz[	64	]=	-0.212	;	J0[	64	]=	-2.592	;
	Dx[	65	]=	0.005	;		Dy[	65	]=	-0.212	;		Dz[	65	]=	0.064	;	J0[	65	]=	-2.592	;
	Dx[	66	]=	-0.034	;		Dy[	66	]=	-0.007	;		Dz[	66	]=	-0.059	;	J0[	66	]=	2.867	;
	Dx[	67	]=	0.059	;		Dy[	67	]=	-0.034	;		Dz[	67	]=	-0.007	;	J0[	67	]=	2.867	;
	Dx[	68	]=	0.007	;		Dy[	68	]=	0.059	;		Dz[	68	]=	-0.034	;	J0[	68	]=	2.867	;
	Dx[	69	]=	-0.007	;		Dy[	69	]=	-0.059	;		Dz[	69	]=	0.034	;	J0[	69	]=	2.867	;
	Dx[	70	]=	-0.059	;		Dy[	70	]=	0.034	;		Dz[	70	]=	0.007	;	J0[	70	]=	2.867	;
	Dx[	71	]=	0.034	;		Dy[	71	]=	0.007	;		Dz[	71	]=	0.059	;	J0[	71	]=	2.867	;
	Dx[	72	]=	-0.328	;		Dy[	72	]=	0.081	;		Dz[	72	]=	0.457	;	J0[	72	]=	15.396	;
	Dx[	73	]=	-0.457	;		Dy[	73	]=	-0.328	;		Dz[	73	]=	-0.081	;	J0[	73	]=	15.396	;
	Dx[	74	]=	0.081	;		Dy[	74	]=	-0.457	;		Dz[	74	]=	-0.328	;	J0[	74	]=	15.396	;
	Dx[	75	]=	0.081	;		Dy[	75	]=	-0.457	;		Dz[	75	]=	0.328	;	J0[	75	]=	15.396	;
	Dx[	76	]=	-0.457	;		Dy[	76	]=	0.328	;		Dz[	76	]=	-0.081	;	J0[	76	]=	15.396	;
	Dx[	77	]=	0.328	;		Dy[	77	]=	0.081	;		Dz[	77	]=	0.457	;	J0[	77	]=	15.396	;
	Dx[	78	]=	-0.089	;		Dy[	78	]=	0.146	;		Dz[	78	]=	-0.009	;	J0[	78	]=	-1.715	;
	Dx[	79	]=	-0.089	;		Dy[	79	]=	-0.146	;		Dz[	79	]=	-0.009	;	J0[	79	]=	-1.715	;
	Dx[	80	]=	0.009	;		Dy[	80	]=	-0.089	;		Dz[	80	]=	0.146	;	J0[	80	]=	-1.715	;
	Dx[	81	]=	0.009	;		Dy[	81	]=	-0.089	;		Dz[	81	]=	-0.146	;	J0[	81	]=	-1.715	;
	Dx[	82	]=	0.146	;		Dy[	82	]=	0.009	;		Dz[	82	]=	0.089	;	J0[	82	]=	-1.715	;
	Dx[	83	]=	-0.146	;		Dy[	83	]=	0.009	;		Dz[	83	]=	0.089	;	J0[	83	]=	-1.715	;
	Dx[	84	]=	0.005	;		Dy[	84	]=	0.212	;		Dz[	84	]=	-0.064	;	J0[	84	]=	-2.592	;
	Dx[	85	]=	0.005	;		Dy[	85	]=	-0.212	;		Dz[	85	]=	-0.064	;	J0[	85	]=	-2.592	;
	Dx[	86	]=	0.212	;		Dy[	86	]=	0.064	;		Dz[	86	]=	-0.005	;	J0[	86	]=	-2.592	;
	Dx[	87	]=	0.064	;		Dy[	87	]=	0.005	;		Dz[	87	]=	0.212	;	J0[	87	]=	-2.592	;
	Dx[	88	]=	-0.212	;		Dy[	88	]=	0.064	;		Dz[	88	]=	-0.005	;	J0[	88	]=	-2.592	;
	Dx[	89	]=	0.064	;		Dy[	89	]=	0.005	;		Dz[	89	]=	-0.212	;	J0[	89	]=	-2.592	;
	Dx[	90	]=	-0.034	;		Dy[	90	]=	0.007	;		Dz[	90	]=	0.059	;	J0[	90	]=	2.867	;
	Dx[	91	]=	-0.059	;		Dy[	91	]=	-0.034	;		Dz[	91	]=	-0.007	;	J0[	91	]=	2.867	;
	Dx[	92	]=	-0.007	;		Dy[	92	]=	0.059	;		Dz[	92	]=	-0.034	;	J0[	92	]=	2.867	;
	Dx[	93	]=	0.007	;		Dy[	93	]=	-0.059	;		Dz[	93	]=	0.034	;	J0[	93	]=	2.867	;
	Dx[	94	]=	0.059	;		Dy[	94	]=	0.034	;		Dz[	94	]=	0.007	;	J0[	94	]=	2.867	;
	Dx[	95	]=	0.034	;		Dy[	95	]=	-0.007	;		Dz[	95	]=	-0.059	;	J0[	95	]=	2.867	;
	// float tmp[3];
	// for (int i=0;i<NeighborPairs;i++)//normalize DMI vectors
	// {
	// 	tmp[0]=Dx[i];
	// 	tmp[1]=Dy[i];
	// 	tmp[2]=Dz[i];
	// 	(void)Unit( tmp, tmp);
	// 	Dx[i]=tmp[0];
	// 	Dy[i]=tmp[1];
	// 	Dz[i]=tmp[2];	
	// }
}

void
SetExchMarkus( float abc[][3], float block[][3], int arrSize, 
	int* aiBlock, int* niBlock, int* niGridA, int* niGridB, int* niGridC, int* sIdx, 
	float* J0, float* B0, float* D0, float* Dx, float* Dy, float* Dz)
{
	double tmp[3];
	J0[	0	]=	0.000	;		
	Dx[	0	]=	0.000	;		
	Dy[	0	]=	0.000	;		
	Dz[	0	]=	0.000	;
	for (int i=0;i<NeighborPairs;i++)//normalize DMI vectors
	{
		tmp[0]=Dx[i];
		tmp[1]=Dy[i];
		tmp[2]=Dz[i];
		(void)Unit( tmp, tmp);
		Dx[i]=tmp[0];
		Dy[i]=tmp[1];
		Dz[i]=tmp[2];	
	}
}

void
SetExchChAch( float abc[][3], float block[][3], int arrSize, 
	int* aiBlock, int* niBlock, int* niGridA, int* niGridB, int* niGridC, int* sIdx, 
	float* J0, float* B0, float* D0, float* Dx, float* Dy, float* Dz)
{

J0[	0	]=	1.000	;		Dx[	0	]=	-0.01	;		Dy[	0	]=	-0.015	;		Dz[	0	]=	0.000	;
J0[	1	]=	1.000	;		Dx[	1	]=	 0.01	;		Dy[	1	]=	-0.015	;		Dz[	1	]=	0.000	;
J0[	2	]=	1.000	;		Dx[	2	]=	-0.01	;		Dy[	2	]=	 0.015	;		Dz[	2	]=	0.000	;
J0[	3	]=	1.000	;		Dx[	3	]=	 0.01	;		Dy[	3	]=	 0.015	;		Dz[	3	]=	0.000	;
// float tmp[3];
// for (int i=0;i<NeighborPairs;i++)//normalize DMI vectors
// {
// 	tmp[0]=Dx[i];
// 	tmp[1]=Dy[i];
// 	tmp[2]=Dz[i];
// 	(void)Unit( tmp, tmp);
// 	Dx[i]=tmp[0];
// 	Dy[i]=tmp[1];
// 	Dz[i]=tmp[2];	
// }
}
