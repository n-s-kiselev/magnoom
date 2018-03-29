//*****************************************************************************
//
//  Project Name : Magnoom
//  Author       : Nikolai S. Kiselev
//  Created      : April 2016
//  Modified     : October 2016
//  
//  Compilation  : 
//    Windows7   :  ???
//
// The standard gcc available on OS X through XCode and Clang doesn't support OpenMP. To install the Homebrew version of gcc with OpenMP support you need to install it with
// install brew: ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
// brew doctor
// sudo chown -R $(whoami):admin /usr/local
//
// brew install gcc --without-multilib
// or as pointed out by @Mark Setchell

// brew reinstall gcc --without-multilib
// This will install it to the /usr/local/bin directory. Homebrew will install it as gcc-<version>so as not to clobber the gcc bundled with XCode. The current gcc version available from Homebrew will install as gcc-4.9. You can compile programs with OpenMP support using it via

// gcc-4.9 -fopenmp hello.c
// Alternatively you could put an alias in your .bashrcfile as

// alias gcc='gcc-4.9'
// and then compile using

// gcc -fopenmp hello.c
//
//
//    Mac OS X   :  
//    Before you can use a dynamic library as a dependent library, 
//    the library and its header files must be installed on your computer. 
//    The standard locations for header files are ~/include, /usr/local/include and /usr/include. 
//    The standard locations for dynamic libraries are ~/lib, /usr/local/lib, and /usr/lib.
//    % g++ main.cpp -o magnoom -O3 -Wall -fno-strict-aliasing -I ./include/ -L ./lib -lAntTweakBar -framework GLUT -pthread -framework OpenGL -Wno-deprecated-declarations
//    % echo "" | gcc -xc - -v -E
//    % sudo cp ~/Documents/Nick/AntTweakBar/lib/libAntTweakBar.dylib /usr/local/lib
//    % sudo cp ~/Documents/Nick/AntTweakBar/include/AntTweakBar.h /usr/local/include/
//    % sudo cp -r ~/Documents/Nick/AntTweakBar /usr/local/
//    % g++ main.cpp -o magnoom -O3 -Wall -fno-strict-aliasing -lAntTweakBar -framework GLUT -pthread -framework OpenGL -Wno-deprecated-declarations
//    Note, in OS X: https://lukecyca.com/2008/glui-235-framework-for-mac-os-x.html
//    DYLD_LIBRARY_PATH=DYLD_LIBRARY_PATH=/usr/local/opt/gcc/lib/gcc/6/    
//
//    Raspberry Pi3:
//    pi@raspberrypi: sudo apt-get install libgl1-mesa-dev libgles2-mesa-dev libglew-dev:armhf libglewmx-dev:armhf libglib2.0-dev libglu1-mesa-dev
//    pi@raspberrypi: g++ main.cpp -o magnoom -pthread -O3 -Wall -fno-strict-aliasing -lAntTweakBar -lpthread  -lglut -lGLU -lGLEW -lGL
//
//    Ubuntu   :  
//    $ g++ main.cpp -o magnoom -pthread -O3 -Wall -fno-strict-aliasing -lAntTweakBar -lpthread  -lglut -lGLU -lGLEW -lGL 
//    $ sudo export PATH="/usr/bin:$PATH" {/usr/local/bin/ld: this linker was not configured to use sysroots}
//    CentOS   : 
//    $ g++ JS2v7.cpp -o JS2v7 -pthread -O3 -Wall -fno-strict-aliasing -I ./include -L ./lib -lAntTweakBar -lpthread  -lglut -lGLU -lGLEW -lGL
//			LD_LIBRARY_PATH=~/Desktop/JSpinx4/lib/
//			LD_LIBRARY_PATH=~/Desktop/JSpinx4/include/
//*****************************************************************************
//
//    Arch Linux:
//    you would need instal: sudo pacman -S mesa glu freeglut 
//    the rest works the same as on Ubuntu 
//

//	Windows 32bit + MSVC
// 	$ "C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\bin\vcvars32.bat"
//	$ cl main.cpp /O2

//	Windows 64bit + MSVC
// 	$ "C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\bin\amd64\vcvars64.bat"
//	$ cl main.cpp /O2



#ifdef __APPLE__
    #include "TargetConditionals.h"
    #ifdef TARGET_OS_MAC
        #include <unistd.h>
        #include <GLUT/glut.h>
        #include <OpenGL/OpenGL.h>
    #endif
#else
	#if defined(_WIN32)
		#pragma comment(lib, "glew32.lib")
		#pragma comment(lib, "freeglut.lib")
	#endif
	#include <GL/glew.h>
	#include <GL/freeglut.h>
#endif

#include <AntTweakBar.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h> //I like printf( )!
#include <time.h>
#include "stereomath.h"
// struct timespec tw = {0,100000000};// delay for 0 s and 10^8 ns =0.1s

#if defined(_WIN32)
	#include <windows.h>
	#include <process.h>
	#define uint32_t  unsigned long int
	#define usleep(t)  Sleep( (DWORD)((t)*0.001f) )
	#define snprintf _snprintf			
	#define SEM_VALUE_MAX 32767
	#define SEM_FAILED NULL
	#define pthread_t HANDLE
	#define semaphore_ref HANDLE
	#define pthread_mutex_t CRITICAL_SECTION
	#define sem_open(name, flag, mode, value) CreateSemaphore(NULL,value,SEM_VALUE_MAX,name)
	#define sem_post(sem) ReleaseSemaphore(sem,1,NULL)
	#define sem_wait(sem) WaitForSingleObject(sem,(DWORD)10000)
	#define sem_trywait(sem) WaitForSingleObject(sem,0)
	#define sem_close(sem) CloseHandle(sem)
	#define pthread_create(th_ref, attr_ref, name, arg_ref) ((  ( *(th_ref) = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)(*(name)), arg_ref, 0, NULL) ) == NULL ? -1 : 0  ))
	#define pthread_mutex_init(mutex_ptr,num) InitializeCriticalSection(mutex_ptr)
	#define pthread_mutex_lock(mutex_ptr) EnterCriticalSection(mutex_ptr)
	#define pthread_mutex_unlock(mutex_ptr) LeaveCriticalSection(mutex_ptr)
#else
	#include <unistd.h>
	#include <pthread.h>
	#include <semaphore.h>		
	#define semaphore_ref sem_t*
#endif

pthread_mutex_t  culc_mutex, show_mutex;

enum engine_mutex_flags{DO_IT,WAIT};
enum data_mutex_flags{WAIT_DATA,TAKE_DATA};

int		ENGINE_MUTEX=WAIT;
int		DATA_TRANSFER_MUTEX=WAIT_DATA;

#define THREADS_NUMBER 3
semaphore_ref	sem_in[THREADS_NUMBER];
semaphore_ref	sem_out[THREADS_NUMBER];

//Alternative version of semaphores
/*
semaphore_ref	sem_A1;
semaphore_ref   sem_A2;

void SyncAllThreads()
	{
	if ( sem_trywait(sem_A1)==0 )
		{
		sem_wait(sem_A2);
		} else
		{
		for (int i=0;i<(THREADS_NUMBER-1);++i)
			{sem_post(sem_A1);}
		for (int i=0;i<(THREADS_NUMBER-1);++i)
			{sem_post(sem_A2);}
		}
	};
*/

#include "global.h"		/*All global variables and constants*/
#include "MATH.cpp"		/*All mathematical fuctions*/
#include "GEOM.cpp"		/*All functions salculating size and neighbors*/
#include "ENGINE.cpp"	/*CALC THREAD:LLG solver*/
#include "OPGL.cpp"		/*VISUAL THREAD: All Visualization Functions*/
#include "INITSTATE.cpp"/*Set of functions for initial states*/
#include <fcntl.h> 		/* For O_CREAT and other O_**** constants */


void ReallocateMemoryForSpins(int NOS){
	Sx = (double *)calloc(NOS, sizeof(double));
	Sy = (double *)calloc(NOS, sizeof(double));
	Sz = (double *)calloc(NOS, sizeof(double));	// <-- for 10^6 spins allocated memory for Sz,Sy,Sz = 12 Mega Byte
	Kind = (int *)calloc(NOS, sizeof(int));
}

void ReallocateMemoryForAllOther(int NOS){
	tSx = (double *)calloc(NOS, sizeof(double));
	tSy = (double *)calloc(NOS, sizeof(double));
	tSz = (double *)calloc(NOS, sizeof(double));	// <-- + 12 Mega Byte

	t2Sx = (double *)calloc(NOS, sizeof(double));
	t2Sy = (double *)calloc(NOS, sizeof(double));
	t2Sz = (double *)calloc(NOS, sizeof(double));	// <-- + 12 Mega Byte

	bSx = (double *)calloc(NOS, sizeof(double));
	bSy = (double *)calloc(NOS, sizeof(double));
	bSz = (double *)calloc(NOS, sizeof(double));	// <-- + 12 Mega Byte

	RNx = (float *)calloc(NOS, sizeof(float));
	RNy = (float *)calloc(NOS, sizeof(float));
	RNz = (float *)calloc(NOS, sizeof(float));	// <-- + 12 Mega Byte

	Px = (float *)calloc(NOS, sizeof(float));
	Py = (float *)calloc(NOS, sizeof(float));
	Pz = (float *)calloc(NOS, sizeof(float));	// <-- + 12 Mega Byte

	BPx = (float *)calloc(NOB, sizeof(float));
	BPy = (float *)calloc(NOB, sizeof(float));
	BPz = (float *)calloc(NOB, sizeof(float));	// <-- + 12 Mega Byte

	Heffx = (double *)calloc(NOS, sizeof(double));
	Heffy = (double *)calloc(NOS, sizeof(double));
	Heffz = (double *)calloc(NOS, sizeof(double));// <-- + 12 Mega Byte

	Etot = (double *)calloc(NOS, sizeof(double));// <-- + 4 Mega Byte
	Etot0= (double *)calloc(NOS, sizeof(double));// <-- + 4 Mega Byte
	// For 100x100x100x10^6 spins total allocated memory is about 6*12 = 72 Mega Byte
	// and possibly may reach up to 100 Mb in total with all other variables.
}

void RestartCalcThreads(pthread_t * thread_id, int * thread_args){
		for (int i=0; i<THREADS_NUMBER; i++){
		thread_args[i] = i;		
		if ( pthread_create(&thread_id[i], NULL, CALC_THREAD, (void *)&thread_args[i]) ) {
			fprintf(stderr, "Error in creating CALC_THREAD thread\n"); 
		}		
	}
}

/*************************************************************************/
/*                        Program Main Thread                            */
/*************************************************************************/

int 
main (int argc, char **argv)
{		
	for (int i=0; i<THREADS_NUMBER; i++){
		char name[10]; 
		snprintf(name,10,"inDoor%d\n",i);
		//printf("%s\n",name);
		if ( (sem_in[i] = sem_open(name, O_CREAT, 0644, 0)) == SEM_FAILED ) perror("sem_open");
		snprintf(name,10,"outDoor%d\n",i);
		//printf("%s\n",name); 
		if ( (sem_out[i] = sem_open(name, O_CREAT, 0644, 0)) == SEM_FAILED ) perror("sem_open");
		//int value;		
		//sem_getvalue(sem_in[i], &value); //Function not implemented on Mac OS X!!!
		Max_torque[i]=0;
	} 

	////////////////////////////////////////////////
	srand ( time(NULL) );//init random number seed//
	////////////////////////////////////////////////

	////////////////////////////////////////////////
	readConfigFile();
	////////////////////////////////////////////////



	GetShells(abc, Block, AtomsPerBlock, ShellNumber, RadiusOfShell);
	for(int i=0;i<ShellNumber;i++) printf("R[%d]=%f\n",i,RadiusOfShell[i] );
	NeighborPairs = GetNeighborsNumber(abc, Block, AtomsPerBlock, ShellNumber, RadiusOfShell, NeighborsPerAtom);
	// Allocate arrays for neighbours map:
	AIdxBlock = (int *)calloc(NeighborPairs, sizeof(int));// index of the atom within the block
	NIdxBlock = (int *)calloc(NeighborPairs, sizeof(int));// index of the neighbour within the block
	NIdxGridA = (int *)calloc(NeighborPairs, sizeof(int));// index of the neighbour within the block
	NIdxGridB = (int *)calloc(NeighborPairs, sizeof(int));// index of the relative position of the block of the neighbour in the greed along tr. vect. a
	NIdxGridC = (int *)calloc(NeighborPairs, sizeof(int));// index of the relative position of the block of the neighbour in the greed along tr. vect. b
	SIdx      = (int *)calloc(NeighborPairs, sizeof(int));// index of the shell corresponding to this pair

	CreateMapOfNeighbors( abc, Block, AtomsPerBlock, ShellNumber, RadiusOfShell, AIdxBlock, NIdxBlock, NIdxGridA, NIdxGridB, NIdxGridC, SIdx);
	Jexc = (float *)calloc(NeighborPairs, sizeof(float));
	Bexc = (float *)calloc(NeighborPairs, sizeof(float));
	Dexc = (float *)calloc(NeighborPairs, sizeof(float));
	VDMx = (float *)calloc(NeighborPairs, sizeof(float));
	VDMy = (float *)calloc(NeighborPairs, sizeof(float));
	VDMz = (float *)calloc(NeighborPairs, sizeof(float));
	SetExch1( abc, Block, NeighborPairs, Jij, Bij, Dij, AIdxBlock, NIdxBlock, NIdxGridA, NIdxGridB, NIdxGridC, SIdx, 
	Jexc, Bexc, Dexc, VDMx, VDMy, VDMz);
	// SetExchMarkus( abc, Block, NeighborPairs, AIdxBlock, NIdxBlock, NIdxGridA, NIdxGridB, NIdxGridC, SIdx, Jexc, Bexc, Dexc, VDMx, VDMy, VDMz);
	// SetExchChAch( abc, Block, NeighborPairs, AIdxBlock, NIdxBlock, NIdxGridA, NIdxGridB, NIdxGridC, SIdx, Jexc, Bexc, Dexc, VDMx, VDMy, VDMz);
	// SetExchMariya( abc, Block, NeighborPairs, Jij, Bij, Dij, AIdxBlock, NIdxBlock, NIdxGridA, NIdxGridB, NIdxGridC, SIdx, Jexc, Bexc, Dexc, VDMx, VDMy, VDMz);
	for(int i=0;i<AtomsPerBlock;i++) printf("Neighbors Per Atom[%d]=%d\n",i,NeighborsPerAtom[i] );
	printf("AtomsPerBlock =%2d\n", AtomsPerBlock);
	printf("  ShellNumber =%2d\n", ShellNumber);
	printf("NeighborPairs =%2d\n", NeighborPairs);
	float AXYZ[3] = {0.0, 0.0, 0.0};//atom position
	float NXYZ[3] = {0.0, 0.0, 0.0};//neighbor position

	// Neighbours map to console:	
	printf("  N\t| Ia\t| In\t| Jn\t| Kn\t| Ln\t| Sn\t| Ax\t| Ay\t| Az\t| Nx\t| Ny\t| Nz\t|Dist\t| Jij\t| Dij\t| D_x\t| D_y\t| D_z\t|\n");
	for(int i=0;i<NeighborPairs;i++) {
		printf("%3d  \t|",            i);
		printf("%3d  \t|", AIdxBlock[i]); //I'
		printf("%3d  \t|", NIdxBlock[i]); //I
		printf("%3d  \t|", NIdxGridA[i]); //J
		printf("%3d  \t|", NIdxGridB[i]); //K
		printf("%3d  \t|", NIdxGridC[i]); //L
		printf("%3d  \t|",      SIdx[i]); //S
		GetPosition( abc, Block, AIdxBlock[i], 0, 0, 0, AXYZ);
		printf("%2.3f\t|%2.3f\t|%2.3f\t|",AXYZ[0],AXYZ[1],AXYZ[2]);
		GetPosition( abc, Block, NIdxBlock[i], NIdxGridA[i], NIdxGridB[i], NIdxGridC[i], NXYZ);
		printf("%2.3f\t|%2.3f\t|%2.3f\t|",NXYZ[0],NXYZ[1],NXYZ[2]);
		AXYZ[0]-=NXYZ[0];
		AXYZ[1]-=NXYZ[1];
		AXYZ[2]-=NXYZ[2];
		printf("%2.3f\t|",  sqrt(AXYZ[0]*AXYZ[0]+AXYZ[1]*AXYZ[1]+AXYZ[2]*AXYZ[2])); //Distance
		printf("%2.3f\t|",  Jexc[i]); //Jij 
		printf("%2.3f\t|",  Dexc[i]); //Jij 
		printf("%2.3f\t|",  VDMx[i]); //DMI x
		printf("%2.3f\t|",  VDMy[i]); //DMI y
		printf("%2.3f\t|\n",VDMz[i]); //DMI z
	}

	// Open output file:
	outFile = fopen ("table.csv","w");
	if (outFile!=NULL) {fputs ("iter,time,Mx,My,Mz,E_tot,\n",outFile);}

	// Memory allocation:
	RHue = (float *)calloc(360, sizeof(float));
	GHue = (float *)calloc(360, sizeof(float));
	BHue = (float *)calloc(360, sizeof(float));

	ReallocateMemoryForSpins(NOS);
	ReallocateMemoryForAllOther(NOS);

	pthread_mutex_init(&culc_mutex,0);
	pthread_mutex_init(&show_mutex,0);

	pthread_t thread_id[THREADS_NUMBER];
	int thread_args[THREADS_NUMBER];
	RestartCalcThreads(thread_id, thread_args);

	GetBox(abc, uABC, Box);
	UpdateSpinPositions(abc, uABC, Block, AtomsPerBlock, Box, Px, Py, Pz);
	UpdateKind(Kind, Px, Py, Pz, NOS, NOSK);
	InitSpinComponents( Px, Py, Pz, Sx, Sy, Sz, 11);
	for (int i=0;i<NOS;i++) { bSx[i]=Sx[i]; bSy[i]=Sy[i]; bSz[i]=Sz[i];}

    // Set OpenGL context initial state.
	setupOpenGL();
    //  Allocate memory for vetices, normals, colors and indicies array used in drawing subrutines
	ReallocateArrayDrawing();
	// Fill array for prototype (arrow or cane) array 
	UpdatePrototypeVerNorInd(vertexProto, normalProto, indicesProto, arrowFaces, WhichVectorMode,0);
	// Fill big array for indecies for all arrows, cans, cones or boxes 
	UpdateIndices(indicesProto , IdNumProto, indices, IdNum, VCNumProto); 
	UpdateVerticesNormalsColors(vertexProto, normalProto, VCNumProto, vertices, normals, colors, VCNum, Px, Py, Pz, bSx, bSy, bSz, WhichVectorMode);
	CreateNewVBO();
	UpdateVBO(&vboIdV, &vboIdN, &vboIdC, &iboIdI, vertices, normals, colors, indices);

    VHf[0]=sin(PI*VHtheta/180)*cos(PI*VHphi/180);
	VHf[1]=sin(PI*VHtheta/180)*sin(PI*VHphi/180);
	VHf[2]=cos(PI*VHtheta/180);

	ReallocateArrayDrawing_H();
    UpdatePrototypeVerNorInd(vertexProto_H, normalProto_H, indices_H, arrowFaces_H, ARROW1,1);
    UpdateVerticesNormalsColors_H(vertexProto_H, normalProto_H, VCNum_H, vertices_H, normals_H, colors_H, Box[0][0]*0.6, Box[1][1]*0.6, Box[2][2]*0.6, VHf[0], VHf[1], VHf[2]);
    CreateNewVBO_H();
    UpdateVBO_H(&vboIdV_H, &vboIdN_H, &vboIdC_H, &iboIdI_H, vertices_H, normals_H, colors_H, indices_H);

	ReallocateArrayDrawing_BOX();
	UpdateVerticesNormalsColors_BOX(vertices_BOX, normals_BOX, colors_BOX, indices_BOX, Box);
	CreateNewVBO_BOX();
	UpdateVBO_BOX(&vboIdV_BOX, &vboIdN_BOX, &vboIdC_BOX, &iboIdI_BOX, vertices_BOX, normals_BOX, colors_BOX, indices_BOX);

	ReallocateArrayDrawing_BASIS();
	UpdateVerticesNormalsColors_BASIS(vertices_BASIS, normals_BASIS, colors_BASIS, indices_BASIS, Box);
	CreateNewVBO_BASIS();
	UpdateVBO_BASIS(&vboIdV_BASIS, &vboIdN_BASIS, &vboIdC_BASIS, &iboIdI_BASIS, vertices_BASIS, normals_BASIS, colors_BASIS, indices_BASIS);

	ReallocateArrayDrawing_PBC();
	UpdateVerticesNormalsColors_PBC(0, vertices_PBC_A, normals_PBC_A, colors_PBC_A, indices_PBC_A, Box);
	UpdateVerticesNormalsColors_PBC(1, vertices_PBC_B, normals_PBC_B, colors_PBC_B, indices_PBC_B, Box);
	UpdateVerticesNormalsColors_PBC(2, vertices_PBC_C, normals_PBC_C, colors_PBC_C, indices_PBC_C, Box);  
	CreateNewVBO_PBC();
    UpdateVBO_PBC(&vboIdV_PBC_A, &vboIdN_PBC_A, &vboIdC_PBC_A, &iboIdI_PBC_A, vertices_PBC_A, normals_PBC_A, colors_PBC_A, indices_PBC_A);
    UpdateVBO_PBC(&vboIdV_PBC_B, &vboIdN_PBC_B, &vboIdC_PBC_B, &iboIdI_PBC_B, vertices_PBC_B, normals_PBC_B, colors_PBC_B, indices_PBC_B);
    UpdateVBO_PBC(&vboIdV_PBC_C, &vboIdN_PBC_C, &vboIdC_PBC_C, &iboIdI_PBC_C, vertices_PBC_C, normals_PBC_C, colors_PBC_C, indices_PBC_C);

  	// Start GLUT event processing loop
	glutMainLoop();

	// Deallocate memory before closing the program:
	free(AIdxBlock);
	free(NIdxBlock);
	free(NIdxGridA);
	free(NIdxGridB);
	free(NIdxGridC);
	free(SIdx);

	free(Jexc);  			free(Bexc);  			free(Dexc);
	free(VDMx);  			free(VDMy);  			free(VDMz);

	free(Sx);    			free(Sy);    			free(Sz);
	free(tSx);   			free(tSy);   			free(tSz);
	free(bSx);   			free(bSy);   			free(bSz);
	free(t2Sx);  			free(t2Sy);  			free(t2Sz);
	free(Heffx); 			free(Heffy); 			free(Heffz);
	free(RNx);   			free(RNy);   			free(RNz);
	free(Px);    			free(Py);    			free(Pz);
	free(BPx);   			free(BPy);   			free(BPz);
	free(RHue);  			free(GHue);  			free(BHue);
	free(vertices);			free(vertices_H); 		free(vertices_BOX); 	free(vertices_PBC_A);
	free(normals);			free(normals_H); 		free(normals_BOX); 		free(normals_PBC_A);
	free(colors);			free(colors_H); 		free(colors_BOX); 		free(colors_PBC_A);
	free(indices);			free(indices_H); 		free(indices_BOX); 		free(indices_PBC_A);
	free(vertexProto);		free(vertexProto_H);
	free(normalProto);		free(normalProto_H);
	free(indicesProto);		free(indicesProto_H);

	free(RadiusOfShell);
	free(NeighborsPerAtom);

	fclose (outFile);
/*
	for (int i=0; i<THREADS_NUMBER; i++){
		if (sem_close(sem_in[i]) == -1) {
		    perror("sem_close");
		    exit(EXIT_FAILURE);
		}
		if (sem_close(sem_out[i]) == -1) {
		    perror("sem_close");
		    exit(EXIT_FAILURE);
		}

	}
*/
	return 0;    
}
