//*****************************************************************************
//
//  Project Name : Magnoom
//  Author       : Nikolai S. Kiselev
//  Created      : April 2016
//  Modified     : October 2016
//  
//  Build with the repository's nob.c build system (see README.md).
#include <glad/glad.h>
#include <AntTweakBar.h>
#include "vendor/glfw2/TwGLFW2.h"

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h> //I like printf( )!
#include <time.h>

#if defined(_WIN32)
	#include <windows.h>
	#include <process.h>
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
	#define pthread_join(thread,retval) (WaitForSingleObject((thread), INFINITE), CloseHandle(thread), 0)
#else
	#include <unistd.h>
	#include <pthread.h>
	#include <semaphore.h>		
	#define semaphore_ref sem_t*
#endif

#define ABS(x) ((x)<0?-(x):(x))

enum engine_mutex_flags{DO_IT,WAIT};
enum data_mutex_flags{WAIT_DATA,TAKE_DATA};

#define THREADS_NUMBER 3

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "vendor/stb/stb_image_write.h"

/*****************************************************************************/
/* Constants formerly in magnoom.h                                          */
/*****************************************************************************/
#define PI      3.14159265359
#define PI2     1.570796326795
#define iPI     0.318309886     // 1/Pi
#define TPI     6.28318530718   // 2*Pi
#define iTPI    0.1591549430919 // 1/(2*Pi)
#define D2R     0.01745329251   // Pi/180 degrees to radians
#define R2D     57.295779513    // 180/Pi radians to degrees
#define min_float ((float)(1e-37F))
#define max_float ((float)(3.402823466e+38F))
#define setNAN(x) *((int *)x) = 0x7FC00000
#define isNAN(x)  ( *((int *)x) == 0x7FC00000  ? TRUE : FALSE  )
#define GoldenRatio 1.61803398875 //(1+sqrt[5])/2
#define intMAX  24000000
// used in OpenGL
#define ESCAPE  0x1b // the escape key:
#define SPACE   0x20 // the space  key:
#define PLUS    0x2b // "+" key
#define MINUS   0x2d // "-" key
// active mouse buttons (or them together):
#define LEFT   4
#define MIDDLE 2
#define RIGHT  1
// number of slots for camera positions to save in memory
#define NumCamPosSave 5

/*****************************************************************************/
/* Enum types needed by magnoom_ctx fields (formerly in magnoom.h/visualization.c) */
/*****************************************************************************/
typedef enum    {A_AXIS, B_AXIS, C_AXIS, FILTER} enSliceMode;
typedef enum    {SIN_FIELD, GAUSSIAN_FIELD, SINC_FIELD, CIRCULAR_FIELD} enACField;
enum            IntegrationScheme{HEUN,SIB,RK23,RK45,RELAX};
enum            Average_mode{ALONG_A,ALONG_B, ALONG_C, ALONG_0};

typedef enum    {ORTHO, PERSP} enProjections;
typedef enum    {RND, HOMO, SKYRM1, SKYRM2, SKYRM3, BOBBER_T, BOBBER_B, BOBBER_L, BOBBER_L_T, BOBBER_L_B,
                 HOPFION1, SPIRAL, SKYRMION_L, GLOBULA, MultyQ, NORM} enIniState;
typedef enum    {DEFAULT_G, CILINDER_G, SPHERE_G} enGeom;
typedef enum    {WHITE, BLACK, RED, GREEN, BLUE, MANUAL} enColors;
typedef enum    {ARROW1, CONE1, CANE, uPOINT, BOX1} enVectorMode;

/*****************************************************************************/
/* magnoom_ctx: single struct holding all program state (formerly globals    */
/* scattered across magnoom.h, magnoom.c, and visualization.c).              */
/*****************************************************************************/
typedef struct magnoom_ctx {
	/* Slicing parameters */
	int             A_layer_min;
	int             B_layer_min;
	int             C_layer_min;
	int             A_layer_max;
	int             B_layer_max;
	int             C_layer_max;
	enSliceMode     WhichSliceMode;

	/* Cartesian components of spins */
	double*         Sx;
	double*         Sy;
	double*         Sz;
	double*         tSx;
	double*         tSy;
	double*         tSz;
	double*         t2Sx;
	double*         t2Sy;
	double*         t2Sz;
	double*         t3Sx;
	double*         t3Sy;
	double*         t3Sz;
	double**        Image_x;
	double**        Image_y;
	double**        Image_z;
	double**        dImage_x;
	double**        dImage_y;
	double**        dImage_z;
	double*         bSx;
	double*         bSy;
	double*         bSz;
	float*          IsoLineX;
	float*          IsoLineY;
	float*          IsoLineZ;
	float*          RNx;
	float*          RNy;
	float*          RNz;
	float*          Px;
	float*          Py;
	float*          Pz;
	float*          BPx;
	float*          BPy;
	float*          BPz;
	int*            Kind;
	double*         Heffx;
	double*         Heffy;
	double*         Heffz;
	bool*           Proj;
	double          ALPHA;
	int             Boundary[3];

	/* Energy */
	double*         Etot;
	double*         Etot0;
	double          totalEnergy;
	double          perSpEnergy;
	double          totalEnergyFerro;
	double          perSpEnergyMinusFerro;
	double          Mtot[3];
	double          mtot[3];
	double          Max_torque[THREADS_NUMBER];
	double          MAX_TORQUE;
	double          coef[THREADS_NUMBER][6];
	double          BigDataBank[10][1500];
	int             recordsCounter;
	int             WhereAmI;
	int             NumOfPoints;
	int             TempPoints;
	int             TempInt;

	/* Color scheme */
	int             HueMapRGB[6];
	int             HueMapRYGB[6];
	int             HueMap[6];
	float*          RHue;
	float*          GHue;
	float*          BHue;

	/* Neighbor map */
	int             NeighborPairs;
	int*            AIdxBlock;
	int*            NIdxBlock;
	int*            NIdxGridA;
	int*            NIdxGridB;
	int*            NIdxGridC;
	int*            SIdx;
	float           Box[3][3];
	float*          Jexc;
	float*          Bexc;
	float*          Dexc;
	float*          VDMx;
	float*          VDMy;
	float*          VDMz;

	/* Heisenberg / biquadratic / DMI exchange (indexed by shell) */
	float           Jij[7];
	float           Bij[6];
	float           Dij[6];

	/* Magnetocrystalline anisotropy */
	float           VKu1[3];
	float           VKu2[3];
	float           Ku1;
	float           Ku2;
	float           Kc;

	/* DC/AC applied field */
	float           VHf[3];
	float           VHtheta;
	float           VHphi;
	float           Hf;
	float           Bdc[3];
	float           VHac[3];
	float           Hac;
	double          Bac[3];
	double          Omega_dc;
	double          Period_dc;
	float           GPulseWidth;
	float           t_offset;
	float           HacTime;
	enACField       WhichACField;

	int             WhichIntegrationScheme;
	int             WhichAverageMode;
	int             save_slice;

	/* Spin-torque / current */
	float           VCu[3];
	float           Cu;

	/* LLG dynamics parameters */
	int             Precession;
	float           damping;
	float           t_step;
	float           Temperature;

	/* Zhang-Li */
	float           Xi;
	float           Curr_u;

	/* GUI control */
	int             Play;
	int             SpecialEvent;
	int             DataTransfer;
	unsigned int    ITERATION;
	unsigned int    Max_Numb_Iteration;
	int             Record;
	int             AC_FIELD_ON;
	int             AC_MODE_REC;

	/* FPS & IPS */
	int             currentTime;
	int             previousTime;
	int             frameCount;
	int             timeInterval;
	int             previousIteration;
	int             currentIteration;
	float           FPS;
	float           IPS;
	FILE*           outFile;
	unsigned int    rec_iteration;
	int             rec_num_mode;
	int             current_rec_num_mode;
	int             num_images;
	char            BuferString[800];
	double          outputEtotal;
	double          outputMtotal[3];
	int             SleepTime;
	char            shortBufer[200];
	char            inputfilename[64];
	char            outputfilename[64];

	/* Geometry */
	float           abc[3][3];
	float           Block[1][3];
	int             uABC[3];
	int             ShellNumber;
	int             AtomsPerBlock;
	float*          RadiusOfShell;
	int*            NeighborsPerAtom;
	int             NOS;
	int             NOS_AL;
	int             NOS_BL;
	int             NOS_CL;
	int             NOSK;
	double          iNOS;
	int             NOB;
	int             NOB_AL;
	int             NOB_BL;
	int             NOB_CL;

	/* Color map control */
	int             ColorShift;
	int             InvertHue;
	int             InvertValue;

	/* magnoom.c: cross-thread concurrency primitives */
	pthread_mutex_t culc_mutex;
	pthread_mutex_t show_mutex;
	int             ENGINE_MUTEX;
	int             DATA_TRANSFER_MUTEX;
	volatile bool   EngineShutdown;
	volatile bool   EngineShutdownRequested;
	semaphore_ref   sem_in[THREADS_NUMBER];
	semaphore_ref   sem_out[THREADS_NUMBER];

	/* visualization.c: window / GLFW / display state */
	GLFWwindow*     MainWindow;
	const char*     WINDOWTITLE;
	int             window_width;
	int             window_height;
	float           asp_rat;
	float           asp_rat_inv;
	double          MouseScaleX;
	double          MouseScaleY;
	int             MouseWheelPosition;

	/* mouse/button state */
	int             ActiveButton;
	int             Xmouse;
	int             Ymouse;

	/* camera / rotation / translation / projection */
	float           CameraEye[3];
	float           CameraC[3];
	float           CameraUp[3];
	float           axisX[3];
	float           axisY[3];
	float           axisZ[3];
	float           PerspSet[4];
	enProjections   WhichProjection;
	float           RotAxis[3];
	float           q_Rotation[4];
	float           Rot[3];
	float           dRot[3];
	float           RotSpeed;
	float           TransXYZ[3];
	float           dTransXYZ[3];
	float           TransSpeed;
	int             CurrentCameraPositionBank;
	float           CameraPosition[NumCamPosSave][7];
	float           axis[3];
	float           angle;

	/* OpenGL material / lighting */
	GLfloat         ambient[4];
	GLfloat         diffuse[4];
	GLfloat         specular[4];
	GLfloat         shininess;
	float           g_LightMultiplier;
	float           g_LightDirection[3];
	int             Light_On;

	/* initial-state generation parameters */
	enIniState      WhichInitialState;
	enGeom          WhichGeometry;
	float           chSizeG;
	float           chSize;
	float           chDir[3];
	float           RotateAllSpins;

	/* color scheme (visualization) */
	int             my_background_color[3];
	enColors        WhichBackgroundColor;
	int             temp_color[3];
	GLfloat         BackgroundColors[6][3];

	/* vector-drawing mode / appearance */
	enVectorMode    WhichVectorMode;
	float           Scale;
	float           Pivot;
	float           WireWidth;
	float           Scale_H;
	int             arrowFaces;
	int             arrowFaces_H;
	int             N_Multisample;

	/* axes/box toggles + legacy display-list handles */
	int             AxesOn;
	int             BoxOn;
	GLuint          AxesList;
	GLuint          BoxList;
	GLuint          BoundaryListA;
	GLuint          BoundaryListB;
	GLuint          BoundaryListC;
	GLfloat         AXES_WIDTH;
	GLfloat         AXES_LENGTH;

	/* AntTweakBar bar handles */
	TwBar*          help_bar;
	TwBar*          view_bar;
	TwBar*          control_bar;
	TwBar*          initial_bar;
	TwBar*          ac_field_bar;
	TwBar*          info_bar;
	TwBar*          my_window;

	/* slice/filter thresholds & flags */
	int             N_filter;
	int             SpinFilter1;
	int             SpinFilter2;
	int             SpinFilter3;
	int             theta_max1;
	float           Sz_min1;
	int             theta_min1;
	float           Sz_max1;
	int             theta_max2;
	float           Sz_min2;
	int             theta_min2;
	float           Sz_max2;
	int             theta_max3;
	float           Sz_min3;
	int             theta_min3;
	float           Sz_max3;
	int             phi_max1;
	int             phi_min1;
	int             phi_max2;
	int             phi_min2;
	int             phi_max3;
	int             phi_min3;
	int             GreedFilter;
	int             GreedFilterInvert;
	int             GreedFilterMaxA;
	int             GreedFilterMinA;
	int             GreedFilterMaxB;
	int             GreedFilterMinB;
	int             GreedFilterMaxC;
	int             GreedFilterMinC;

	/* drawing arrays: vertices/normals/colors/indices and variants */
	GLfloat*        vertexProto;
	GLfloat*        normalProto;
	GLuint*         indicesProto;
	GLfloat*        vertexProto_H;
	GLfloat*        normalProto_H;
	GLuint*         indicesProto_H;
	GLfloat*        vertices;
	GLfloat*        normals;
	GLfloat*        colors;
	GLuint*         indices;
	GLfloat*        vertices_H;
	GLfloat*        normals_H;
	GLfloat*        colors_H;
	GLuint*         indices_H;
	GLfloat*        vertices_BOX;
	GLfloat*        normals_BOX;
	GLfloat*        colors_BOX;
	GLuint*         indices_BOX;
	GLfloat*        vertices_BASIS;
	GLfloat*        normals_BASIS;
	GLfloat*        colors_BASIS;
	GLuint*         indices_BASIS;
	GLfloat*        vertices_AC_phase;
	GLfloat*        colors_AC_phase;
	GLuint*         indices_AC_phase;
	GLfloat*        vertices_PBC_A;
	GLfloat*        normals_PBC_A;
	GLfloat*        colors_PBC_A;
	GLuint*         indices_PBC_A;
	GLfloat*        vertices_PBC_B;
	GLfloat*        normals_PBC_B;
	GLfloat*        colors_PBC_B;
	GLuint*         indices_PBC_B;
	GLfloat*        vertices_PBC_C;
	GLfloat*        normals_PBC_C;
	GLfloat*        colors_PBC_C;
	GLuint*         indices_PBC_C;

	/* VBO/IBO ids */
	GLuint          vboIdV;
	GLuint          vboIdN;
	GLuint          vboIdC;
	GLuint          iboIdI;
	GLuint          vboIdV_H;
	GLuint          vboIdN_H;
	GLuint          vboIdC_H;
	GLuint          iboIdI_H;
	GLuint          vboIdV_BOX;
	GLuint          vboIdN_BOX;
	GLuint          vboIdC_BOX;
	GLuint          iboIdI_BOX;
	GLuint          vboIdV_BASIS;
	GLuint          vboIdN_BASIS;
	GLuint          vboIdC_BASIS;
	GLuint          iboIdI_BASIS;
	GLuint          vboIdV_AC_phase;
	GLuint          vboIdC_AC_phase;
	GLuint          iboIdI_AC_phase;
	GLuint          vboIdV_PBC_A, vboIdV_PBC_B, vboIdV_PBC_C;
	GLuint          vboIdN_PBC_A, vboIdN_PBC_B, vboIdN_PBC_C;
	GLuint          vboIdC_PBC_A, vboIdC_PBC_B, vboIdC_PBC_C;
	GLuint          iboIdI_PBC_A, iboIdI_PBC_B, iboIdI_PBC_C;

	/* geometry-sizing counters */
	int             ElNumProto;
	int             IdNumProto;
	int             VCNumProto;
	int             ElNum;
	int             IdNum;
	int             VCNum;
	int             IdNum_H;
	int             VCNum_H;
	int             IdNum_BOX;
	int             VCNum_BOX;
	int             IdNum_BASIS;
	int             VCNum_BASIS;
	int             IdNum_AC_phase;
	int             VCNum_AC_phase;
	int             IdNum_PBC;
	int             VCNum_PBC;
} magnoom_ctx;

magnoom_ctx mag_ctx;

/* Per-thread argument for CALC_THREAD (solvers.c): carries both the thread's */
/* index and ctx, since pthread's start-routine signature has room for only  */
/* one void* argument. */
typedef struct { int index; magnoom_ctx *ctx; } calc_thread_arg;

/*****************************************************************************/
/* Transitional macros: every migrated name NOT reused as a function         */
/* parameter/local elsewhere resolves transparently to its mag_ctx field, so */
/* not-yet-converted code keeps compiling and behaving identically while the */
/* refactor proceeds file-by-file. Removed once every reference uses ctx->   */
/* explicitly (see docs/plans/magnoom-ctx-refactor.md, section 9.8).         */
/* Excluded on purpose (real name shadowing elsewhere in the codebase):      */
/* abc, uABC, NOS, Jij, Bij, Dij, VDMx, VDMy, VDMz, Etot, Mtot, Xi, Curr_u,   */
/* coef, ALPHA, axisX, axisY, axisZ, axis, window_width, window_height, Px,  */
/* Py, Pz, Sx, Sy, Sz, Kind, vertices, normals, colors, indices, Rot,        */
/* WhichSliceMode.                                                          */
/*****************************************************************************/

/*****************************************************************************/
/* magnoom_ctx_init: sets every field to its current compile-time default,   */
/* folding in the former InitializeGlobalState()'s derived-value computation.*/
/* Returns false (matching InitializeGlobalState's contract) if the neighbor-*/
/* map allocations fail.                                                    */
/*****************************************************************************/
bool magnoom_ctx_init(magnoom_ctx *ctx)
{
	/* NOTE: ctx is always ctx (there is exactly one instance). Fields   */
	/* that still resolve through a transitional macro (see the macro block  */
	/* above) MUST be written as bare names here, not ctx->name -- the macro */
	/* expansion after "->" would otherwise produce invalid syntax. Only the */
	/* names excluded from macros (real shadows elsewhere: abc, uABC, NOS,   */
	/* Jij, Bij, Dij, axisX/Y/Z, axis, window_width, window_height,          */
	/* WhichSliceMode) are written through ctx-> below.                      */

	/* Slicing parameters */
	ctx->A_layer_min = 1;
	ctx->B_layer_min = 1;
	ctx->C_layer_min = 1;
	ctx->A_layer_max = 1;
	ctx->B_layer_max = 1;
	ctx->C_layer_max = 1;
	ctx->WhichSliceMode = C_AXIS;

	/* Boundary conditions */
	ctx->Boundary[0] = 1; ctx->Boundary[1] = 1; ctx->Boundary[2] = 0;

	/* Color scheme */
	ctx->HueMapRGB[0]=0;  ctx->HueMapRGB[1]=60;  ctx->HueMapRGB[2]=120;
	ctx->HueMapRGB[3]=180; ctx->HueMapRGB[4]=240; ctx->HueMapRGB[5]=300;
	ctx->HueMapRYGB[0]=0;  ctx->HueMapRYGB[1]=90;  ctx->HueMapRYGB[2]=180;
	ctx->HueMapRYGB[3]=225; ctx->HueMapRYGB[4]=270; ctx->HueMapRYGB[5]=315;
	ctx->HueMap[0]=0;  ctx->HueMap[1]=60;  ctx->HueMap[2]=120;
	ctx->HueMap[3]=180; ctx->HueMap[4]=240; ctx->HueMap[5]=300;

	/* Heisenberg / biquadratic / DMI exchange */
	ctx->Jij[0]=1.0f; ctx->Jij[1]=0.0f; ctx->Jij[2]=0.0f; ctx->Jij[3]=0.0f;
	ctx->Jij[4]=0.0f; ctx->Jij[5]=0.0f; ctx->Jij[6]=0.0f;
	ctx->Bij[0]=0.0f; ctx->Bij[1]=0.0f; ctx->Bij[2]=0.0f;
	ctx->Bij[3]=0.0f; ctx->Bij[4]=0.0f; ctx->Bij[5]=0.0f;
	ctx->Dij[0]=0.251327412f; ctx->Dij[1]=0.0f; ctx->Dij[2]=0.0f;
	ctx->Dij[3]=0.0f; ctx->Dij[4]=0.0f; ctx->Dij[5]=0.0f;

	/* Magnetocrystalline anisotropy */
	ctx->VKu1[0]=0.0f; ctx->VKu1[1]=0.0f; ctx->VKu1[2]=1.0f;
	ctx->VKu2[0]=0.0f; ctx->VKu2[1]=0.0f; ctx->VKu2[2]=1.0f;

	/* DC/AC applied field */
	ctx->VHf[0]=0.0f; ctx->VHf[1]=0.0f; ctx->VHf[2]=1.0f;
	ctx->VHac[0]=0.0f; ctx->VHac[1]=0.0f; ctx->VHac[2]=1.0f;
	ctx->Omega_dc = 0.005;
	ctx->GPulseWidth = 20.0f;
	ctx->WhichACField = SIN_FIELD;

	ctx->WhichIntegrationScheme = SIB;
	ctx->WhichAverageMode = ALONG_0;

	/* Spin-torque / current */
	ctx->VCu[0]=0.0f; ctx->VCu[1]=0.0f; ctx->VCu[2]=1.0f;

	/* LLG dynamics parameters */
	ctx->Precession = 1;
	ctx->damping = 0.01f;
	ctx->t_step = 0.1f;

	/* GUI control */
	ctx->SpecialEvent = 1;
	ctx->DataTransfer = 1;
	ctx->Max_Numb_Iteration = 100000;

	/* FPS & IPS / recording */
	ctx->rec_iteration = 100;
	ctx->rec_num_mode = 10;
	ctx->num_images = 32;
	ctx->SleepTime = 1000;
	strcpy(ctx->inputfilename, "input.csv");
	strcpy(ctx->outputfilename, "output.csv");

	/* Geometry */
	ctx->abc[0][0]=1.0f; ctx->abc[0][1]=0.0f; ctx->abc[0][2]=0.0f;
	ctx->abc[1][0]=0.0f; ctx->abc[1][1]=1.0f; ctx->abc[1][2]=0.0f;
	ctx->abc[2][0]=0.0f; ctx->abc[2][1]=0.0f; ctx->abc[2][2]=1.0f;
	ctx->Block[0][0]=0.5f; ctx->Block[0][1]=0.5f; ctx->Block[0][2]=0.5f;
	ctx->uABC[0]=10; ctx->uABC[1]=10; ctx->uABC[2]=10;
	ctx->ShellNumber = 1;
	ctx->AtomsPerBlock = sizeof(ctx->Block)/sizeof(float)/3;

	/* magnoom.c concurrency primitives */
	ctx->ENGINE_MUTEX = WAIT;
	ctx->DATA_TRANSFER_MUTEX = WAIT_DATA;
	ctx->EngineShutdown = false;
	ctx->EngineShutdownRequested = false;

	/* visualization.c: window / GLFW / display state */
	ctx->WINDOWTITLE = "Magnoom v1.0";
	ctx->window_width = 1400;
	ctx->window_height = 800;
	ctx->MouseScaleX = 1.0;
	ctx->MouseScaleY = 1.0;

	/* camera / rotation / translation / projection */
	ctx->CameraEye[0]=0.0f; ctx->CameraEye[1]=0.0f; ctx->CameraEye[2]=150.0f;
	ctx->CameraC[0]=0.0f; ctx->CameraC[1]=0.0f; ctx->CameraC[2]=0.0f;
	ctx->CameraUp[0]=0.0f; ctx->CameraUp[1]=1.0f; ctx->CameraUp[2]=0.0f;
	ctx->axisX[0]=1.0f; ctx->axisX[1]=0.0f; ctx->axisX[2]=0.0f;
	ctx->axisY[0]=0.0f; ctx->axisY[1]=1.0f; ctx->axisY[2]=0.0f;
	ctx->axisZ[0]=0.0f; ctx->axisZ[1]=0.0f; ctx->axisZ[2]=1.0f;
	ctx->PerspSet[0]=60.0f; ctx->PerspSet[1]=0.0f; ctx->PerspSet[2]=1.0f; ctx->PerspSet[3]=5000.0f;
	ctx->WhichProjection = PERSP;
	ctx->q_Rotation[0]=0.0f; ctx->q_Rotation[1]=0.0f; ctx->q_Rotation[2]=0.0f; ctx->q_Rotation[3]=1.0f;
	ctx->RotSpeed = 1.0f;
	ctx->TransSpeed = 3.0f;
	ctx->Scale = 1.5f;
	ctx->Pivot = 0.55f;
	ctx->WireWidth = 0.2f;
	ctx->Scale_H = 10.0f;
	ctx->axis[0]=0.7f; ctx->axis[1]=0.7f; ctx->axis[2]=0.0f;
	ctx->angle = 0.8f;

	/* OpenGL material / lighting */
	ctx->ambient[0]=0.33f; ctx->ambient[1]=0.22f; ctx->ambient[2]=0.03f; ctx->ambient[3]=0.0f;
	ctx->diffuse[0]=0.78f; ctx->diffuse[1]=0.57f; ctx->diffuse[2]=0.11f; ctx->diffuse[3]=1.0f;
	ctx->specular[0]=0.1f; ctx->specular[1]=0.1f; ctx->specular[2]=0.08f; ctx->specular[3]=1.0f;
	ctx->shininess = 200.0f;
	ctx->g_LightMultiplier = 1.0f;
	ctx->g_LightDirection[0]=0.5f; ctx->g_LightDirection[1]=0.5f; ctx->g_LightDirection[2]=-0.5f;
	ctx->Light_On = 1;

	/* initial-state generation parameters */
	ctx->WhichInitialState = RND;
	ctx->WhichGeometry = DEFAULT_G;
	ctx->chSizeG = 50.0f;
	ctx->chSize = 12.0f;
	ctx->chDir[0]=0.0f; ctx->chDir[1]=1.0f; ctx->chDir[2]=0.0f;

	/* color scheme (visualization) */
	ctx->my_background_color[0]=55; ctx->my_background_color[1]=55; ctx->my_background_color[2]=155;
	ctx->WhichBackgroundColor = MANUAL;
	ctx->temp_color[0]=55; ctx->temp_color[1]=55; ctx->temp_color[2]=155;
	ctx->BackgroundColors[0][0]=1;    ctx->BackgroundColors[0][1]=1;    ctx->BackgroundColors[0][2]=1;
	ctx->BackgroundColors[1][0]=0.1f; ctx->BackgroundColors[1][1]=0.1f; ctx->BackgroundColors[1][2]=0.1f;
	ctx->BackgroundColors[2][0]=1.0f; ctx->BackgroundColors[2][1]=0.8f; ctx->BackgroundColors[2][2]=0.8f;
	ctx->BackgroundColors[3][0]=0.8f; ctx->BackgroundColors[3][1]=1.0f; ctx->BackgroundColors[3][2]=0.8f;
	ctx->BackgroundColors[4][0]=0.8f; ctx->BackgroundColors[4][1]=0.8f; ctx->BackgroundColors[4][2]=1.0f;
	ctx->BackgroundColors[5][0]=0.3f; ctx->BackgroundColors[5][1]=0.3f; ctx->BackgroundColors[5][2]=0.3f;

	/* vector-drawing mode / appearance */
	ctx->WhichVectorMode = BOX1;
	ctx->arrowFaces = 6;
	ctx->arrowFaces_H = 30;
	ctx->N_Multisample = 16;

	/* axes/box toggles */
	ctx->AxesOn = 1;
	ctx->BoxOn = 1;
	ctx->AXES_WIDTH = 2.0f;
	ctx->AXES_LENGTH = 2.0f;

	/* slice/filter thresholds & flags */
	ctx->SpinFilter1 = 1;
	ctx->theta_max1 = 95; ctx->theta_min1 = 85;
	ctx->theta_max2 = 180; ctx->theta_min2 = 160;
	ctx->theta_max3 = 10; ctx->theta_min3 = 0;
	ctx->phi_max1 = 360; ctx->phi_max2 = 360; ctx->phi_max3 = 360;

	/*************************************************************/
	/* Computed defaults (formerly InitializeGlobalState's body). */
	/*************************************************************/
	ctx->Bdc[0] = ctx->Hf*ctx->VHf[0];
	ctx->Bdc[1] = ctx->Hf*ctx->VHf[1];
	ctx->Bdc[2] = ctx->Hf*ctx->VHf[2];
	ctx->Hac = ctx->Hf*(float)sin(PI/180.0);
	ctx->Period_dc = TPI/ctx->Omega_dc;

	ctx->RadiusOfShell = (float *)calloc((size_t)ctx->ShellNumber, sizeof(float));
	ctx->NeighborsPerAtom = (int *)calloc((size_t)ctx->AtomsPerBlock, sizeof(int));
	if (ctx->RadiusOfShell == NULL || ctx->NeighborsPerAtom == NULL) return false;

	ctx->NOS = ctx->AtomsPerBlock*ctx->uABC[0]*ctx->uABC[1]*ctx->uABC[2];
	ctx->NOS_AL = ctx->AtomsPerBlock*ctx->uABC[1]*ctx->uABC[2];
	ctx->NOS_BL = ctx->AtomsPerBlock*ctx->uABC[0]*ctx->uABC[2];
	ctx->NOS_CL = ctx->AtomsPerBlock*ctx->uABC[0]*ctx->uABC[1];
	ctx->iNOS = 1.0/ctx->NOS;
	ctx->NOB = ctx->uABC[0]*ctx->uABC[1]*ctx->uABC[2];
	ctx->NOB_AL = ctx->uABC[1]*ctx->uABC[2];
	ctx->NOB_BL = ctx->uABC[0]*ctx->uABC[2];
	ctx->NOB_CL = ctx->uABC[0]*ctx->uABC[1];

	ctx->asp_rat = (float)ctx->window_width/(float)ctx->window_height;
	ctx->asp_rat_inv = (float)ctx->window_height/(float)ctx->window_width;
	ctx->PerspSet[1] = ctx->asp_rat;
	ctx->Sz_min1 = (float)cos(ctx->theta_max1*PI/180.0);
	ctx->Sz_max1 = (float)cos(ctx->theta_min1*PI/180.0);
	ctx->Sz_min2 = (float)cos(ctx->theta_max2*PI/180.0);
	ctx->Sz_max2 = (float)cos(ctx->theta_min2*PI/180.0);
	ctx->Sz_min3 = (float)cos(ctx->theta_max3*PI/180.0);
	ctx->Sz_max3 = (float)cos(ctx->theta_min3*PI/180.0);
	ctx->GreedFilterMaxA = ctx->uABC[0]-1;
	ctx->GreedFilterMaxB = ctx->uABC[1]-1;
	ctx->GreedFilterMaxC = ctx->uABC[2]-1;
	return true;
}

/*****************************************************************************/
/* Functions formerly embedded in magnoom.h (color/HSV conversion and        */
/* OVF/VTK/PNG/BIN file I/O). Left with their original parameter lists       */
/* (including Sx/Sy/Sz/Px/Py/Pz/box/WhichSliceMode parameters that already   */
/* correctly shadow the corresponding globals) for now; every OTHER global   */
/* they touch is a non-shadowed name and keeps resolving through the         */
/* transitional macros above, exactly like every other not-yet-converted     */
/* function (see docs/plans/magnoom-ctx-refactor.md, section 9.8).           */
/*****************************************************************************/

float color_function(float H ){
    float h, result=0;
    if (H>360){
        h=H-360.f;
    }else{
        h=H;
    }
    if (60<=h && h<180){
        result=1.0f;
    } else if (240<=h && h<360){
        result=0.0f;
    }
    if (0<=h && h<60){
        result=h/60.0f;
    }else if (180<=h && h<240){
        result=4.0f-h/60.0f;
    }
    return result;
}


void HSVtoRGB(magnoom_ctx *ctx, float Vec[3], float rgb[3], int inV, int inH ){
    // int H = inH*359+(1-2*inH)*((int) atan2int (Vec[1]+0.01,Vec[0])); //it's fast atan function [int deg], see math_utils.c
    float S=sqrt(Vec[0]*Vec[0]+Vec[1]*Vec[1]+Vec[2]*Vec[2]);
    float F = atan2(Vec[1]/S,Vec[0]/S);
    float H = inH*360+(1-2*inH)*(F > 0 ? F : (TPI + F))*R2D;
    float maxV, minV, dV;
    if ((1-2*inV)*Vec[2]/S<0){
        // maxV = 1 - fabs(Vec[2]/S);
        maxV = 1 - fabs(Vec[2]/S);
        minV=0;
    }else{
        maxV = 1;
        minV = fabs(Vec[2]/S);
    }
    //metka
    maxV = 0.5+S*(maxV-0.5);
    minV = 0.5-S*(0.5-minV);


    dV = maxV-minV;

    // RGB[0] = RHue[H]*dV+minV;
    // RGB[1] = GHue[H]*dV+minV;
    // RGB[2] = BHue[H]*dV+minV;
    //metka: for director field color code//
    //dV = 1 - fabs(Vec[2]);              //
    //minV =0.5*(1-dV);                   //
    ////////////////////////////////////////


    rgb[0] = color_function(H+120+ctx->ColorShift)*dV+minV;//rad
    rgb[1] = color_function(H+000+ctx->ColorShift)*dV+minV;//green
    rgb[2] = color_function(H-120+ctx->ColorShift)*dV+minV;//blue

    // rgb[0] = colorGreen(H+120+ColorShift)*dV*S+minV+dV*(1-S)*0.5;//rad
    // rgb[1] = colorGreen(H    +ColorShift)*dV*S+minV+dV*(1-S)*0.5;//green
    // rgb[2] = colorGreen(H-120+ColorShift)*dV*S+minV+dV*(1-S)*0.5;//blue
}




void Save_OVF_b8(magnoom_ctx *ctx, double* Sx, double* Sy, double* Sz, char ovf_filename[64]){
    float temp0 = 0;
    float temp1 = 0;
    float temp2 = 0;
    float temp3 = 0;
    float a_lattice = 1.0e-9;
    //char ovf_filename[64] = "mode.ovf";
    FILE * pFile = fopen (ovf_filename,"wb");
    if(pFile!=NULL) {
         fputs ("# OOMMF OVF 2.0\n",pFile);
         fputs ("# Segment count: 1\n",pFile);
         fputs ("# Begin: Segment\n",pFile);
         fputs ("# Begin: Header\n",pFile);
         fputs ("# Title: m\n",pFile);
         fputs ("# meshtype: rectangular\n",pFile);
         fputs ("# meshunit: m\n",pFile);
         fputs ("# xmin: 0\n",pFile);
         fputs ("# ymin: 0\n",pFile);
         fputs ("# zmin: 0\n",pFile);

         temp0  = ctx->abc[0][0]*ctx->abc[0][0];
         temp0 += ctx->abc[0][1]*ctx->abc[0][1];
         temp0 += ctx->abc[0][2]*ctx->abc[0][2];
         temp1 = sqrt(temp0);
         snprintf(ctx->shortBufer,80,"# xmax: %.6g\n",ctx->uABC[0]*temp1*a_lattice);
         fputs (ctx->shortBufer,pFile);

         temp0  = ctx->abc[1][0]*ctx->abc[1][0];
         temp0 += ctx->abc[1][1]*ctx->abc[1][1];
         temp0 += ctx->abc[1][2]*ctx->abc[1][2];
         temp2 = sqrt(temp0);
         snprintf(ctx->shortBufer,80,"# ymax: %.6g\n",ctx->uABC[1]*temp2*a_lattice);
         fputs (ctx->shortBufer,pFile);

         temp0  = ctx->abc[2][0]*ctx->abc[2][0];
         temp0 += ctx->abc[2][1]*ctx->abc[2][1];
         temp0 += ctx->abc[2][2]*ctx->abc[2][2];
         temp3 = sqrt(temp0);
         snprintf(ctx->shortBufer,80,"# zmax: %.6g\n",ctx->uABC[2]*temp3*a_lattice);
         fputs (ctx->shortBufer,pFile);
         fputs ("# valuedim: 3\n",pFile);
         fputs ("# valuelabels: m_x m_y m_z\n",pFile);
         fputs ("# valueunits: 1 1 1\n",pFile);
         fputs ("# Desc: Total simulation time:  0  s\n",pFile);

         snprintf(ctx->shortBufer,80,"# xbase: %.6g\n",temp1*0.5*a_lattice);
         fputs (ctx->shortBufer,pFile);

         snprintf(ctx->shortBufer,80,"# ybase: %.6g\n",temp2*0.5*a_lattice);
         fputs (ctx->shortBufer,pFile);

         snprintf(ctx->shortBufer,80,"# zbase: %.6g\n",temp3*0.5*a_lattice);
         fputs (ctx->shortBufer,pFile);

         snprintf(ctx->shortBufer,80,"# xnodes: %d\n",ctx->uABC[0]);
         fputs (ctx->shortBufer,pFile);
         snprintf(ctx->shortBufer,80,"# ynodes: %d\n",ctx->uABC[1]);
         fputs (ctx->shortBufer,pFile);
         snprintf(ctx->shortBufer,80,"# znodes: %d\n",ctx->uABC[2]);
         fputs (ctx->shortBufer,pFile);

         snprintf(ctx->shortBufer,80,"# xstepsize:  %.6g\n",temp1*a_lattice);
         fputs (ctx->shortBufer,pFile);

         snprintf(ctx->shortBufer,80,"# ystepsize: %.6g\n",temp2*a_lattice);
         fputs (ctx->shortBufer,pFile);

         snprintf(ctx->shortBufer,80,"# zstepsize: %.6g\n",temp3*a_lattice);
         fputs (ctx->shortBufer,pFile);

         fputs ("# End: Header\n",pFile);
         fputs ("# Begin: Data Binary 8\n",pFile);
         double Temp1[]= {123456789012345.0};
         fwrite (Temp1, sizeof(double), 1, pFile);
         for (int cn = 0; cn<ctx->uABC[2]; cn++)
            {
             for (int bn = 0; bn<ctx->uABC[1]; bn++)
                {
                 for (int an = 0; an<ctx->uABC[0]; an++)
                    {
                     int n = an+bn*ctx->uABC[0]+cn*ctx->uABC[0]*ctx->uABC[1];// index of the block
                     n = n*ctx->AtomsPerBlock;//index of the first spin in the block
                     for (int atom=0; atom<ctx->AtomsPerBlock; atom++)
                        {
                         int N = n + atom;
                         double Temp[]= {Sx[N], Sy[N], Sz[N]};
                         fwrite (Temp , sizeof(double), 3, pFile);
                     }
                 }
             }
         }
        fputs ("# End: Data Binary 4\n",pFile);
        fputs ("# End: Segment\n",pFile);
        fclose (pFile);
    }
    printf("Recording to the file %s is done!\n", ovf_filename);
}


void SaveBin(magnoom_ctx *ctx, double* Sx, double* Sy, double* Sz, char bin_filename[64]){
    unsigned short int num = 65535;
    struct tfshortint {
        unsigned short int t;
        unsigned short int f;
    };
    FILE * pFile = fopen (bin_filename,"wb");
    for(int k = 0; k < ctx->uABC[2]; k++){
        for(int j = 0; j < ctx->uABC[1]; j++){
            for(int i = 0; i < ctx->uABC[0]; i++){
            int n = i+j*ctx->uABC[0]+k*ctx->uABC[0]*ctx->uABC[1];

            double nx = Sx[n], ny = Sy[n], nz = Sz[n];

            double T, F;
            T = acos(nz)/PI;
            F = atan2(ny,nx)/PI;
            if(F <= 0) F += 2.0;
            F /= 2;

            unsigned short int p=0, q=0;

            q = T*num;
            p = F*num;

            struct tfshortint my_par = {q, p};
            fwrite(&my_par, sizeof(struct tfshortint), 1, pFile);

      }
    }
  }
  fclose (pFile);
  printf("Recording to the file %s is done!\n", bin_filename);
}


void SavePng(magnoom_ctx *ctx, double* Sx, double* Sy, double* Sz, char png_filename[64], enSliceMode WhichSliceMode, int x1, int y1, int z1){
  if (WhichSliceMode!=FILTER)
  {
    int scale = 10;
    int N1 = ctx->uABC[0], N2 = ctx->uABC[1];
    switch (WhichSliceMode){
      case A_AXIS:
          N1 = ctx->uABC[1]; N2 = ctx->uABC[2];
          break;
      case B_AXIS:
          N1 = ctx->uABC[0]; N2 = ctx->uABC[2];
          break;
      case C_AXIS:
          N1 = ctx->uABC[0]; N2 = ctx->uABC[1];
          break;
      default: N1 = ctx->uABC[0]; N2 = ctx->uABC[1];
          break;
    }

    int width = scale*N1;
    int height = scale*N2;
    unsigned char* image = (unsigned char*)calloc((size_t)width*height*3, sizeof(unsigned char));
    if (image == NULL) {
      printf("Unable to allocate memory for PNG image %s\n", png_filename);
      return;
    }

    float rgb[3],vec[3];
    for(int i = 0; i < N1; i++){
      for(int j = 0; j < N2; j++){
        int n;
        switch (WhichSliceMode){
          case A_AXIS:
              n = x1 + i*ctx->uABC[0] + j*ctx->uABC[0]*ctx->uABC[1];
              break;
          case B_AXIS:
              n = N1-i-1 + y1*ctx->uABC[0] + j*ctx->uABC[0]*ctx->uABC[1];
              break;
          case C_AXIS:
              n = i + j*ctx->uABC[0] + z1*ctx->uABC[0]*ctx->uABC[1];
              break;
          default: n = i + j*ctx->uABC[0] + 0*ctx->uABC[0]*ctx->uABC[1];
              break;
        }
        rgb[0] = 0; rgb[1] = 0; rgb[2] = 0;
        vec[0] = Sx[n]; vec[1] = Sy[n]; vec[2] = Sz[n];
        // vec[0] = 0; vec[1] = 0; vec[2] =1;
        HSVtoRGB(ctx, vec, rgb, ctx->InvertValue, ctx->InvertHue);///1,0
        rgb[0] *= 255;
        rgb[1] *= 255;
        rgb[2] *= 255;

        int ii = scale*i;
        for(int k1 = 0;k1<scale;k1++)
            for(int k2 = 0;k2<scale; k2++) {
                int x = ii+k1;
                int y = scale*((N2-1)-j)+k2;
                size_t pixel = ((size_t)y*width+x)*3;
                image[pixel+0] = (unsigned char)rgb[0];
                image[pixel+1] = (unsigned char)rgb[1];
                image[pixel+2] = (unsigned char)rgb[2];
            }
      }
    }
    if (stbi_write_png(png_filename, width, height, 3, image, width*3))
      printf("Image has been saved to %s\n",png_filename);
    else
      printf("Unable to save PNG image to %s\n",png_filename);
    free(image);
  }
}

void Save_VTS_b4(magnoom_ctx *ctx, double* Sx, double* Sy, double* Sz, float * Px, float * Py, float * Pz, float box[][3], char vts_filename[64]){
    float a_lattice = 1.0e-9;
    FILE * pFile = fopen (vts_filename,"wb");
    if(pFile!=NULL) {
        fputs ("<?xml version=\"1.0\"?>\n",pFile);
        fputs ("<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n",pFile);
        snprintf(ctx->shortBufer,80,"\t<StructuredGrid WholeExtent=\"0 %d 0 %d 0 %d\">\n",ctx->uABC[0]-1, ctx->uABC[1]-1, ctx->uABC[2]-1);
        fputs (ctx->shortBufer,pFile);
        snprintf(ctx->shortBufer,80,"\t\t<Piece Extent=\"0 %d 0 %d 0 %d \">\n",ctx->uABC[0]-1, ctx->uABC[1]-1, ctx->uABC[2]-1);
        fputs (ctx->shortBufer,pFile);
        fputs ("\t\t\t<PointData Vectors=\"m\">\n",pFile);
        fputs ("\t\t\t\t<DataArray type=\"Float32\" Name=\"m\" NumberOfComponents=\"3\" format=\"binary\">\n",pFile);
        fputs ("\t\t\t\t\t",pFile);
        float Tr[3] = {
                        (box[0][0]+box[1][0]+box[2][0])/2.f,
                        (box[0][1]+box[1][1]+box[2][1])/2.f,
                        (box[0][2]+box[1][2]+box[2][2])/2.f
                      };
        for (int cn = 0; cn<ctx->uABC[2]; cn++){
            for (int bn = 0; bn<ctx->uABC[1]; bn++){
                for (int an = 0; an<ctx->uABC[0]; an++){
                    int n = an+bn*ctx->uABC[0]+cn*ctx->uABC[0]*ctx->uABC[1];// index of the block
                    n = n*ctx->AtomsPerBlock;//index of the first spin in the block
                    for (int atom=0; atom<ctx->AtomsPerBlock; atom++){
                        int N = n + atom;
                        float Temp[]= {(Px[N]+Tr[0])*a_lattice, (Py[N]+Tr[1])*a_lattice, (Pz[N]+Tr[2])*a_lattice};
                        fwrite ((const char*)Temp , sizeof(int), 3, pFile);
                    }
                }
            }
        }
        fputs ("\n",pFile);
        fputs ("\t\t\t\t</DataArray>\n",pFile);
        fputs ("\t\t\t</PointData>\n",pFile);
        fputs ("\t\t\t<Points>\n",pFile);
        fputs ("\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n",pFile);
        fputs ("\t\t\t\t\t",pFile);
        for (int cn = 0; cn<ctx->uABC[2]; cn++){
            for (int bn = 0; bn<ctx->uABC[1]; bn++){
                for (int an = 0; an<ctx->uABC[0]; an++){
                    int n = an+bn*ctx->uABC[0]+cn*ctx->uABC[0]*ctx->uABC[1];// index of the block
                    n = n*ctx->AtomsPerBlock;//index of the first spin in the block
                    for (int atom=0; atom<ctx->AtomsPerBlock; atom++){
                        int N = n + atom;
                        float Temp[]= {(float)Sx[N], (float)Sy[N], (float)Sz[N]};
                        fwrite ((const char*)Temp , sizeof(int), 3, pFile);
                    }
                }
            }
        }
        fputs ("\n",pFile);
        fputs ("\t\t\t\t</DataArray>\n",pFile);
        fputs ("\t\t\t</Points>\n",pFile);
        fputs ("\t\t</Piece>\n",pFile);
        fputs ("\t</StructuredGrid>\n",pFile);
        fputs ("</VTKFile>\n",pFile);
        fclose (pFile);
    }
    printf("Recording to the file %s is done!\n", vts_filename);
}

void Save_VTS_ascii(magnoom_ctx *ctx, double* Sx, double* Sy, double* Sz, float * Px, float * Py, float * Pz, float box[][3], char vts_filename[64]){
    float a_lattice = 1.0;//.0e-9;
    FILE * pFile = fopen (vts_filename,"wb");
    if(pFile!=NULL) {
        fputs ("<?xml version=\"1.0\"?>\n",pFile);
        fputs ("<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n",pFile);
        snprintf(ctx->shortBufer,80,"\t<StructuredGrid WholeExtent=\"0 %d 0 %d 0 %d\">\n",ctx->uABC[0]-1, ctx->uABC[1]-1, ctx->uABC[2]-1);
        fputs (ctx->shortBufer,pFile);
        snprintf(ctx->shortBufer,80,"\t\t<Piece Extent=\"0 %d 0 %d 0 %d \">\n",ctx->uABC[0]-1, ctx->uABC[1]-1, ctx->uABC[2]-1);
        fputs (ctx->shortBufer,pFile);
        fputs ("\t\t\t<PointData Vectors=\"m\">\n",pFile);
        fputs ("\t\t\t\t<DataArray type=\"Float32\" Name=\"m\" NumberOfComponents=\"3\" format=\"ascii\">\n",pFile);
        fputs ("\t\t\t\t\t",pFile);
        float Tr[3] = {
                        (box[0][0]+box[1][0]+box[2][0])/2.f,
                        (box[0][1]+box[1][1]+box[2][1])/2.f,
                        (box[0][2]+box[1][2]+box[2][2])/2.f
                      };
        for (int cn = 0; cn<ctx->uABC[2]; cn++){
            for (int bn = 0; bn<ctx->uABC[1]; bn++){
                for (int an = 0; an<ctx->uABC[0]; an++){
                    int n = an+bn*ctx->uABC[0]+cn*ctx->uABC[0]*ctx->uABC[1];// index of the block
                    n = n*ctx->AtomsPerBlock;//index of the first spin in the block
                    for (int atom=0; atom<ctx->AtomsPerBlock; atom++){
                        int N = n + atom;
                        snprintf(ctx->shortBufer,80,"%.6g %.6g %.6g ",(Px[N]+Tr[0])*a_lattice, (Py[N]+Tr[1])*a_lattice, (Pz[N]+Tr[2])*a_lattice);
                        fputs (ctx->shortBufer,pFile);
                    }
                }
            }
        }
        fputs ("\n",pFile);
        fputs ("\t\t\t\t</DataArray>\n",pFile);
        fputs ("\t\t\t</PointData>\n",pFile);
        fputs ("\t\t\t<Points>\n",pFile);
        fputs ("\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n",pFile);
        fputs ("\t\t\t\t\t",pFile);
        for (int cn = 0; cn<ctx->uABC[2]; cn++){
            for (int bn = 0; bn<ctx->uABC[1]; bn++){
                for (int an = 0; an<ctx->uABC[0]; an++){
                    int n = an+bn*ctx->uABC[0]+cn*ctx->uABC[0]*ctx->uABC[1];// index of the block
                    n = n*ctx->AtomsPerBlock;//index of the first spin in the block
                    for (int atom=0; atom<ctx->AtomsPerBlock; atom++){
                        int N = n + atom;
                        snprintf(ctx->shortBufer,80,"%.6g %.6g %.6g ",Sx[N], Sy[N], Sz[N]);
                        fputs (ctx->shortBufer,pFile);
                    }
                }
            }
        }
        fputs ("\n",pFile);
        fputs ("\t\t\t\t</DataArray>\n",pFile);
        fputs ("\t\t\t</Points>\n",pFile);
        fputs ("\t\t</Piece>\n",pFile);
        fputs ("\t</StructuredGrid>\n",pFile);
        fputs ("</VTKFile>\n",pFile);
        fclose (pFile);
    }
    printf("Recording to the file %s is done!\n", vts_filename);
}

/// Swap endianness of a value
float end_swap_float(const float val, const int size)
    {
    float ret = 0;
    for(int i = 0; i < size; ++i)
        ((char*)&ret)[size-1-i] = ((char*)&val)[i];
    return ret;
    }

double end_swap_double(const double val, const int size)
    {
    double ret = 0;
    for(int i = 0; i < size; ++i)
        ((char*)&ret)[size-1-i] = ((char*)&val)[i];
    return ret;
    }

void Save_VTK(magnoom_ctx *ctx, double* Sx, double* Sy, double* Sz, const int mode, char vtk_filename[64])
{
    float temp0 = 0;
    float temp1 = 0;
    float temp2 = 0;
    float temp3 = 0;
    float a_lattice = 1.0e-9;
    FILE * pFile = fopen (vtk_filename,"wb");
    if(pFile!=NULL) {
        fputs ("# vtk DataFile Version 2.0\n",pFile);
        fputs ("Field data file\n",pFile);
        snprintf(ctx->shortBufer,80,"%s\n",(mode==0 ? "BINARY" : "ASCII" ));
        fputs (ctx->shortBufer,pFile);
        fputs ("\n",pFile);
        fputs ("DATASET STRUCTURED_POINTS\n",pFile);
        snprintf(ctx->shortBufer,80,"DIMENSIONS %d %d %d \n",ctx->uABC[0], ctx->uABC[1], ctx->uABC[2]);
        fputs (ctx->shortBufer,pFile);
        fputs ("ORIGIN 0 0 0\n",pFile);

        temp0  = ctx->abc[0][0]*ctx->abc[0][0];
        temp0 += ctx->abc[0][1]*ctx->abc[0][1];
        temp0 += ctx->abc[0][2]*ctx->abc[0][2];
        temp1 = sqrt(temp0);

        temp0  = ctx->abc[1][0]*ctx->abc[1][0];
        temp0 += ctx->abc[1][1]*ctx->abc[1][1];
        temp0 += ctx->abc[1][2]*ctx->abc[1][2];
        temp2 = sqrt(temp0);

        temp0  = ctx->abc[2][0]*ctx->abc[2][0];
        temp0 += ctx->abc[2][1]*ctx->abc[2][1];
        temp0 += ctx->abc[2][2]*ctx->abc[2][2];
        temp3 = sqrt(temp0);

        snprintf(ctx->shortBufer,80,"SPACING %.6g %.6g %.6g \n",temp1, temp2, temp3);
        fputs (ctx->shortBufer,pFile);

        snprintf(ctx->shortBufer,80,"POINT_DATA %d \n",ctx->uABC[0]*ctx->uABC[1]*ctx->uABC[2]);
        fputs (ctx->shortBufer,pFile);
        fputs ("\n",pFile);
        fputs ("SCALARS m float 3\n",pFile);
        fputs ("LOOKUP_TABLE default\n",pFile);
        for (int cn = 0; cn<ctx->uABC[2]; cn++){
            for (int bn = 0; bn<ctx->uABC[1]; bn++){
                for (int an = 0; an<ctx->uABC[0]; an++){
                    int n = an+bn*ctx->uABC[0]+cn*ctx->uABC[0]*ctx->uABC[1];// index of the block
                    n = n*ctx->AtomsPerBlock;//index of the first spin in the block
                    temp0 = 0.f;
                    temp1 = 0.f;
                    temp2 = 0.f;
                    for (int atom=0; atom<ctx->AtomsPerBlock; atom++){
                        int N = n + atom;
                        temp0 += Sx[N];
                        temp1 += Sy[N];
                        temp2 += Sz[N];
                    }
                    temp0 = temp0/ctx->AtomsPerBlock;
                    temp1 = temp1/ctx->AtomsPerBlock;
                    temp2 = temp2/ctx->AtomsPerBlock;
                    if (mode==0){
                        temp0 = end_swap_float(temp0, sizeof(temp0));
                        temp1 = end_swap_float(temp1, sizeof(temp1));
                        temp2 = end_swap_float(temp2, sizeof(temp2));
                        fwrite ((const char*)&temp0, sizeof(float), 1, pFile);
                        fwrite ((const char*)&temp1, sizeof(float), 1, pFile);
                        fwrite ((const char*)&temp2, sizeof(float), 1, pFile);
                    }else{
                        snprintf(ctx->shortBufer,80,"%.6g %.6g %.6g ",temp0, temp1, temp2);
                        fputs (ctx->shortBufer,pFile);
                    }
                }
            }
        }
        fputs ("\n",pFile);
        fputs ("SCALARS F1 float 1\n",pFile);
        fputs ("LOOKUP_TABLE default\n",pFile);
        for (int cn = 0; cn<ctx->uABC[2]; cn++){
            for (int bn = 0; bn<ctx->uABC[1]; bn++){
                for (int an = 0; an<ctx->uABC[0]; an++){
                    int n = an+bn*ctx->uABC[0]+cn*ctx->uABC[0]*ctx->uABC[1];// index of the block
                    n = n*ctx->AtomsPerBlock;//index of the first spin in the block
                    temp0 = 0.f;
                    temp1 = 0.f;
                    temp2 = 0.f;
                    for (int atom=0; atom<ctx->AtomsPerBlock; atom++){
                        int N = n + atom;
                        temp0 += Sx[N];
                        temp1 += Sy[N];
                        temp2 += Sz[N];
                    }
                        temp0 = temp0/ctx->AtomsPerBlock;
                        temp1 = temp1/ctx->AtomsPerBlock;
                        temp2 = temp2/ctx->AtomsPerBlock;
                        temp3 = sqrt(temp0*temp0+temp1*temp1+temp2*temp2);
                        temp0 = temp0/temp3;
                        temp1 = temp1/temp3;
                        temp2 = temp2/temp3;
                        temp0 = atan2(temp1,temp0) * 180 / PI;
                    if (mode==0){
                        temp0 = end_swap_float(temp0, sizeof(temp0));
                        fwrite ((const char*)&temp0, sizeof(float), 1, pFile);
                    }else{
                        snprintf(ctx->shortBufer,80,"%.6g ",temp0);
                        fputs (ctx->shortBufer,pFile);
                    }
                }
            }
        }
        fputs ("\n",pFile);
        fputs ("SCALARS F2 float 1\n",pFile);
        fputs ("LOOKUP_TABLE default\n",pFile);
        for (int cn = 0; cn<ctx->uABC[2]; cn++){
            for (int bn = 0; bn<ctx->uABC[1]; bn++){
                for (int an = 0; an<ctx->uABC[0]; an++){
                    int n = an+bn*ctx->uABC[0]+cn*ctx->uABC[0]*ctx->uABC[1];// index of the block
                    n = n*ctx->AtomsPerBlock;//index of the first spin in the block
                    temp0 = 0.f;
                    temp1 = 0.f;
                    temp2 = 0.f;
                    for (int atom=0; atom<ctx->AtomsPerBlock; atom++){
                        int N = n + atom;
                        temp0 += Sx[N];
                        temp1 += Sy[N];
                        temp2 += Sz[N];
                    }
                        temp0 = temp0/ctx->AtomsPerBlock;
                        temp1 = temp1/ctx->AtomsPerBlock;
                        temp2 = temp2/ctx->AtomsPerBlock;
                        temp3 = sqrt(temp0*temp0+temp1*temp1+temp2*temp2);
                        temp0 = temp0/temp3;
                        temp1 = temp1/temp3;
                        temp2 = temp2/temp3;
                        temp0 = atan2(temp1,-temp0) * 180 / PI;
                    if (mode==0){
                        temp0 = end_swap_float(temp0, sizeof(temp0));
                        fwrite ((const char*)&temp0, sizeof(float), 1, pFile);
                    }else{
                        snprintf(ctx->shortBufer,80,"%.6g ",temp0);
                        fputs (ctx->shortBufer,pFile);
                    }
                }
            }
        }
        fclose (pFile);
        printf("Recording to the file %s is done!\n", vtk_filename);
    }
}

void Save_VTK_6(magnoom_ctx *ctx, double* Sx, double* Sy, double* Sz,
                double* dSx, double* dSy, double* dSz, const int mode, char vtk_filename[64])
{
    float temp0 = 0;
    float temp1 = 0;
    float temp2 = 0;
    float temp3 = 0;
    float a_lattice = 1.0e-9;
    FILE * pFile = fopen (vtk_filename,"wb");
    if(pFile!=NULL) {
        fputs ("# vtk DataFile Version 2.0\n",pFile);
        fputs ("Field data file\n",pFile);
        snprintf(ctx->shortBufer,80,"%s\n",(mode==0 ? "BINARY" : "ASCII" ));
        fputs (ctx->shortBufer,pFile);
        fputs ("\n",pFile);
        fputs ("DATASET STRUCTURED_POINTS\n",pFile);
        snprintf(ctx->shortBufer,80,"DIMENSIONS %d %d %d \n",ctx->uABC[0], ctx->uABC[1], ctx->uABC[2]);
        fputs (ctx->shortBufer,pFile);
        fputs ("ORIGIN 0 0 0\n",pFile);

        temp0  = ctx->abc[0][0]*ctx->abc[0][0];
        temp0 += ctx->abc[0][1]*ctx->abc[0][1];
        temp0 += ctx->abc[0][2]*ctx->abc[0][2];
        temp1 = sqrt(temp0);

        temp0  = ctx->abc[1][0]*ctx->abc[1][0];
        temp0 += ctx->abc[1][1]*ctx->abc[1][1];
        temp0 += ctx->abc[1][2]*ctx->abc[1][2];
        temp2 = sqrt(temp0);

        temp0  = ctx->abc[2][0]*ctx->abc[2][0];
        temp0 += ctx->abc[2][1]*ctx->abc[2][1];
        temp0 += ctx->abc[2][2]*ctx->abc[2][2];
        temp3 = sqrt(temp0);

        snprintf(ctx->shortBufer,80,"SPACING %.6g %.6g %.6g \n",temp1*a_lattice, temp2*a_lattice, temp3*a_lattice);
        fputs (ctx->shortBufer,pFile);

        snprintf(ctx->shortBufer,80,"POINT_DATA %d \n",ctx->uABC[0]*ctx->uABC[1]*ctx->uABC[2]);
        fputs (ctx->shortBufer,pFile);
        fputs ("\n",pFile);
        fputs ("SCALARS m float 3\n",pFile);
        fputs ("LOOKUP_TABLE default\n",pFile);
        for (int cn = 0; cn<ctx->uABC[2]; cn++){
            for (int bn = 0; bn<ctx->uABC[1]; bn++){
                for (int an = 0; an<ctx->uABC[0]; an++){
                    int n = an+bn*ctx->uABC[0]+cn*ctx->uABC[0]*ctx->uABC[1];// index of the block
                    n = n*ctx->AtomsPerBlock;//index of the first spin in the block
                    temp0 = 0.f;
                    temp1 = 0.f;
                    temp2 = 0.f;
                    for (int atom=0; atom<ctx->AtomsPerBlock; atom++){
                        int N = n + atom;
                        temp0 += Sx[N];
                        temp1 += Sy[N];
                        temp2 += Sz[N];
                    }
                    temp0 = temp0/ctx->AtomsPerBlock;
                    temp1 = temp1/ctx->AtomsPerBlock;
                    temp2 = temp2/ctx->AtomsPerBlock;
                    if (mode==0){
                        temp0 = end_swap_float(temp0, sizeof(temp0));
                        temp1 = end_swap_float(temp1, sizeof(temp1));
                        temp2 = end_swap_float(temp2, sizeof(temp2));
                        fwrite ((const char*)&temp0, sizeof(float), 1, pFile);
                        fwrite ((const char*)&temp1, sizeof(float), 1, pFile);
                        fwrite ((const char*)&temp2, sizeof(float), 1, pFile);
                    }else{
                        snprintf(ctx->shortBufer,80,"%.6g %.6g %.6g ",temp0, temp1, temp2);
                        fputs (ctx->shortBufer,pFile);
                    }
                }
            }
        }

        fputs ("\n",pFile);
        fputs ("SCALARS dm float 3\n",pFile);
        fputs ("LOOKUP_TABLE default\n",pFile);
        for (int cn = 0; cn<ctx->uABC[2]; cn++){
            for (int bn = 0; bn<ctx->uABC[1]; bn++){
                for (int an = 0; an<ctx->uABC[0]; an++){
                    int n = an+bn*ctx->uABC[0]+cn*ctx->uABC[0]*ctx->uABC[1];// index of the block
                    n = n*ctx->AtomsPerBlock;//index of the first spin in the block
                    temp0 = 0.f;
                    temp1 = 0.f;
                    temp2 = 0.f;
                    for (int atom=0; atom<ctx->AtomsPerBlock; atom++){
                        int N = n + atom;
                        temp0 += dSx[N];
                        temp1 += dSy[N];
                        temp2 += dSz[N];
                    }
                    temp0 = temp0/ctx->AtomsPerBlock;
                    temp1 = temp1/ctx->AtomsPerBlock;
                    temp2 = temp2/ctx->AtomsPerBlock;
                    if (mode==0){
                        temp0 = end_swap_float(temp0, sizeof(temp0));
                        temp1 = end_swap_float(temp1, sizeof(temp1));
                        temp2 = end_swap_float(temp2, sizeof(temp2));
                        fwrite ((const char*)&temp0, sizeof(float), 1, pFile);
                        fwrite ((const char*)&temp1, sizeof(float), 1, pFile);
                        fwrite ((const char*)&temp2, sizeof(float), 1, pFile);
                    }else{
                        snprintf(ctx->shortBufer,80,"%.6g %.6g %.6g ",temp0, temp1, temp2);
                        fputs (ctx->shortBufer,pFile);
                    }
                }
            }
        }

        fputs ("\n",pFile);
        fputs ("SCALARS F1 float 1\n",pFile);
        fputs ("LOOKUP_TABLE default\n",pFile);
        for (int cn = 0; cn<ctx->uABC[2]; cn++){
            for (int bn = 0; bn<ctx->uABC[1]; bn++){
                for (int an = 0; an<ctx->uABC[0]; an++){
                    int n = an+bn*ctx->uABC[0]+cn*ctx->uABC[0]*ctx->uABC[1];// index of the block
                    n = n*ctx->AtomsPerBlock;//index of the first spin in the block
                    temp0 = 0.f;
                    temp1 = 0.f;
                    temp2 = 0.f;
                    for (int atom=0; atom<ctx->AtomsPerBlock; atom++){
                        int N = n + atom;
                        temp0 += dSx[N];
                        temp1 += dSy[N];
                        temp2 += dSz[N];
                    }
                        temp0 = temp0/ctx->AtomsPerBlock;
                        temp1 = temp1/ctx->AtomsPerBlock;
                        temp2 = temp2/ctx->AtomsPerBlock;
                        temp3 = sqrt(temp0*temp0+temp1*temp1+temp2*temp2);
                        temp0 = temp0/temp3;
                        temp1 = temp1/temp3;
                        temp2 = temp2/temp3;
                        temp0 = atan2(temp1,temp0) * 180 / PI;
                    if (mode==0){
                        temp0 = end_swap_float(temp0, sizeof(temp0));
                        fwrite ((const char*)&temp0, sizeof(float), 1, pFile);
                    }else{
                        snprintf(ctx->shortBufer,80,"%.6g ",temp0);
                        fputs (ctx->shortBufer,pFile);
                    }
                }
            }
        }
        fputs ("\n",pFile);
        fputs ("SCALARS F2 float 1\n",pFile);
        fputs ("LOOKUP_TABLE default\n",pFile);
        for (int cn = 0; cn<ctx->uABC[2]; cn++){
            for (int bn = 0; bn<ctx->uABC[1]; bn++){
                for (int an = 0; an<ctx->uABC[0]; an++){
                    int n = an+bn*ctx->uABC[0]+cn*ctx->uABC[0]*ctx->uABC[1];// index of the block
                    n = n*ctx->AtomsPerBlock;//index of the first spin in the block
                    temp0 = 0.f;
                    temp1 = 0.f;
                    temp2 = 0.f;
                    for (int atom=0; atom<ctx->AtomsPerBlock; atom++){
                        int N = n + atom;
                        temp0 += Sx[N];
                        temp1 += Sy[N];
                        temp2 += Sz[N];
                    }
                        temp0 = temp0/ctx->AtomsPerBlock;
                        temp1 = temp1/ctx->AtomsPerBlock;
                        temp2 = temp2/ctx->AtomsPerBlock;
                        temp3 = sqrt(temp0*temp0+temp1*temp1+temp2*temp2);
                        temp0 = temp0/temp3;
                        temp1 = temp1/temp3;
                        temp2 = temp2/temp3;
                        temp0 = atan2(temp1,-temp0) * 180 / PI;
                    if (mode==0){
                        temp0 = end_swap_float(temp0, sizeof(temp0));
                        fwrite ((const char*)&temp0, sizeof(float), 1, pFile);
                    }else{
                        snprintf(ctx->shortBufer,80,"%.6g ",temp0);
                        fputs (ctx->shortBufer,pFile);
                    }
                }
            }
        }
        fclose (pFile);
        printf("Recording to the file %s is done!\n", vtk_filename);
    }
}

int ReadHeaderLine(FILE * fp, char * line)
{
    char c;//single character
    int pos=0;
    do{ c = (char)fgetc(fp);//get current char and move pointer to the next position
        if (c != EOF && c != '\n') { line[pos++] = c;}//if it's not the end of the file
    }while(c != EOF && c != '\n');//if it's not the end of the file or end of the line
    line[pos] = 0;//complite the readed line
    if ((pos==0 || line[0]!='#') && c != EOF){
        return ReadHeaderLine(fp, line);// recursive call for ReadHeaderLine if the current line is empty
    }
    return pos-1;// the last symbol is the line end symbol
}

void ReadDataLine(FILE * fp, char * line)
{
    char c;//single character
    int pos=0;
    do{ c = (char)fgetc(fp);
        if (c != EOF && c != '\n') { line[pos++] = c;}
    }while(c != EOF && c != '\n');
    line[pos] = 0;
}

int ReadVTKLines(FILE * fp, char * line)
{
    char c;//single character
    int pos=0;
    do{ c = (char)fgetc(fp);//get current char and move pointer to the next position
        if (c != EOF && c != '\n') { line[pos++] = c;}//if it's not the end of the file
    }while(c != EOF && c != '\n');//if it's not the end of the file or end of the line
    line[pos] = 0;//complite the readed line
    if ((pos==0 || line[0]==' ') && c != EOF){
        return ReadVTKLines(fp, line);// recursive call for ReadHeaderLine if the current line is empty
    }
    return pos-1;// the last symbol is the line end symbol
}

void Read_VTK(magnoom_ctx *ctx, double* Sx, double* Sy, double* Sz, char vtk_filename[64]){
char  line[256];//whole line of header should be not longer then 256 characters
int   lineLength=0;
int   valuedim=3;
int   xnodes;
int   ynodes;
int   znodes;
char  keyW1 [256];//key word 1
char  keyW2 [256];//key word 2
char  keyW3 [256];//key word 3
char  keyW4 [256];//key word 4
int   binType = 0;
float temp4_x, temp4_y, temp4_z;
double temp8_x, temp8_y, temp8_z;
FILE * FilePointer = fopen(vtk_filename, "rb");
if(FilePointer!=NULL) {
        lineLength=ReadHeaderLine(FilePointer, line);//read and check the first nonempty line which starts with '#'
        if (lineLength==-1) {// if there are no one line which starts with '#'
            printf("%s has unknown format! \n", vtk_filename);
        }else{
            sscanf(line, "# %s %s %s %s", keyW1, keyW2, keyW3, keyW4);
            if(strncmp(keyW1, "vtk",3)!=0 ||
               strncmp(keyW2, "DataFile",  8)!=0 ||
               strncmp(keyW3, "Version",  7)!=0 ||
               strncmp(keyW4, "2.0",  3)!=0){
                //if the first line isn't "vtk DataFile Version 2.0"

                printf("%s has wrong header or wrong file format! (not vtk 2.0 file format):\n", ctx->inputfilename);
                printf("%s %s %s %s\n", keyW1, keyW2, keyW3, keyW4);
                lineLength=-1;
            }
        }
        //READING HEADER
        if (lineLength!=-1){
            do{
                lineLength = ReadVTKLines(FilePointer, line);
                sscanf(line, "%s %s %s %s", keyW1, keyW2, keyW3, keyW4 );
                //printf("%s %s %s\n", keyW1, keyW2, keyW3);
                if (strncmp(keyW1, "DIMENSIONS",10)==0) {
                    sscanf(keyW2, "%d", &xnodes ); printf("xnodes=%d\n", xnodes);
                    sscanf(keyW3, "%d", &ynodes ); printf("ynodes=%d\n", ynodes);
                    sscanf(keyW4, "%d", &znodes ); printf("znodes=%d\n", znodes);
                }else if (strncmp(keyW1, "SCALARS",7)==0) {
                    printf("DataName=%s\n", keyW2);
                    if (strncmp(keyW3, "float",5)==0 && strncmp(keyW4, "3",1)==0) binType = 4;
                    if (strncmp(keyW3, "double",6)==0 && strncmp(keyW4, "3",1)==0) binType = 8;
                }
            }while(!(strncmp(keyW1, "LOOKUP_TABLE",12)==0 && strncmp(keyW2, "default",7)==0) && lineLength != -1 );
        }
        //READING DATA
        if (valuedim!=0 && xnodes!=0 && ynodes!=0 && znodes!=0){
            sscanf(line, "#%*s %s %s %s", keyW1, keyW2, keyW3 );
            int n;
            if (binType == 4){
                for (int k=0; k<znodes; k++){
                    for (int j=0; j<ynodes; j++){
                        for (int i=0; i<xnodes; i++){
                            if (k<ctx->uABC[2] && j<ctx->uABC[1] && i<ctx->uABC[0]){
                                n = i + j*xnodes + k*xnodes*ynodes; //index of the block!
                                //printf("n=%d\n", n);
                                if (binType==4){
                                    if(!fread(&temp4_x,binType,1,FilePointer)) break;
                                    if(!fread(&temp4_y,binType,1,FilePointer)) break;
                                    if(!fread(&temp4_z,binType,1,FilePointer)) break;
                                    for (int t=0; t<ctx->AtomsPerBlock; t++){
                                        int I=n*ctx->AtomsPerBlock+t;
                                        temp4_x = end_swap_float(temp4_x, sizeof(temp4_x));
                                        temp4_y = end_swap_float(temp4_y, sizeof(temp4_y));
                                        temp4_z = end_swap_float(temp4_z, sizeof(temp4_z));
                                        Sx[I]=ctx->bSx[I]=(double)temp4_x;
                                        Sy[I]=ctx->bSy[I]=(double)temp4_y;
                                        Sz[I]=ctx->bSz[I]=(double)temp4_z;
                                    }
                                }else if (binType==8){
                                    if(!fread(&temp8_x,binType,1,FilePointer)) break;
                                    if(!fread(&temp8_y,binType,1,FilePointer)) break;
                                    if(!fread(&temp8_z,binType,1,FilePointer)) break;
                                    for (int t=0; t<ctx->AtomsPerBlock; t++){
                                        int I=n*ctx->AtomsPerBlock+t;
                                        temp8_x = end_swap_double(temp8_x, sizeof(temp8_x));
                                        temp8_y = end_swap_double(temp8_y, sizeof(temp8_y));
                                        temp8_z = end_swap_double(temp8_z, sizeof(temp8_z));
                                        Sx[I]=ctx->bSx[I]=temp8_x;
                                        Sy[I]=ctx->bSy[I]=temp8_y;
                                        Sz[I]=ctx->bSz[I]=temp8_z;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }else{
            printf("%s cannot read data in vtk file!\n", ctx->inputfilename);
        }
        fclose(FilePointer);
        printf("Done!\n");// when everything is done
    }else{printf("Cannot open file: %s \n", ctx->inputfilename);}
}

#include "linmath.h"		/*All global variables and constants*/

void ReallocateMemoryForImages(magnoom_ctx *ctx, int NumImages, int NOS){
    ctx->Image_x = (double **) calloc(NumImages, sizeof(double *));
    for(int i=0;i<NumImages;i++) ctx->Image_x[i] = (double *) calloc(NOS, sizeof(double));
    ctx->Image_y = (double **) calloc(NumImages, sizeof(double *));
    for(int i=0;i<NumImages;i++) ctx->Image_y[i] = (double *) calloc(NOS, sizeof(double));
    ctx->Image_z = (double **) calloc(NumImages, sizeof(double *));
    for(int i=0;i<NumImages;i++) ctx->Image_z[i] = (double *) calloc(NOS, sizeof(double));

    ctx->dImage_x = (double **) calloc(NumImages, sizeof(double *));
    for(int i=0;i<NumImages;i++) ctx->dImage_x[i] = (double *) calloc(NOS, sizeof(double));
    ctx->dImage_y = (double **) calloc(NumImages, sizeof(double *));
    for(int i=0;i<NumImages;i++) ctx->dImage_y[i] = (double *) calloc(NOS, sizeof(double));
    ctx->dImage_z = (double **) calloc(NumImages, sizeof(double *));
    for(int i=0;i<NumImages;i++) ctx->dImage_z[i] = (double *) calloc(NOS, sizeof(double));
}

#include "math_utils.c"		/*All mathematical fuctions*/
#include "lattice_geometry.c"	/*All functions salculating size and neighbors*/
#include "solvers.c"		/*CALC THREAD: LLG solvers*/
#include "visualization.c"	/*VISUAL THREAD: All Visualization Functions*/
#include "initial_states.c"	/*Set of functions for initial states*/

#include <fcntl.h> 		/* For O_CREAT and other O_**** constants */

void ReallocateMemoryForSpins(magnoom_ctx *ctx, int NOS){
	ctx->Sx = (double *)calloc(NOS, sizeof(double));
	ctx->Sy = (double *)calloc(NOS, sizeof(double));
	ctx->Sz = (double *)calloc(NOS, sizeof(double));	// <-- for 10^6 spins allocated memory for ctx->Sz,ctx->Sy,ctx->Sz = 12 Mega Byte
	ctx->Kind = (int *)calloc(NOS, sizeof(int));
}

void ReallocateMemoryForAllOther(magnoom_ctx *ctx, int NOS){
	ctx->tSx = (double *)calloc(NOS, sizeof(double));
	ctx->tSy = (double *)calloc(NOS, sizeof(double));
	ctx->tSz = (double *)calloc(NOS, sizeof(double));	// <-- + 12 Mega Byte

	ctx->t2Sx = (double *)calloc(NOS, sizeof(double));
	ctx->t2Sy = (double *)calloc(NOS, sizeof(double));
	ctx->t2Sz = (double *)calloc(NOS, sizeof(double));	// <-- + 12 Mega Byte

	ctx->t3Sx = (double *)calloc(NOS, sizeof(double));
	ctx->t3Sy = (double *)calloc(NOS, sizeof(double));
	ctx->t3Sz = (double *)calloc(NOS, sizeof(double));	// <-- + 12 Mega Byte

	ctx->bSx = (double *)calloc(NOS, sizeof(double));
	ctx->bSy = (double *)calloc(NOS, sizeof(double));
	ctx->bSz = (double *)calloc(NOS, sizeof(double));	// <-- + 12 Mega Byte

	ctx->RNx = (float *)calloc(NOS, sizeof(float));
	ctx->RNy = (float *)calloc(NOS, sizeof(float));
	ctx->RNz = (float *)calloc(NOS, sizeof(float));	// <-- + 12 Mega Byte

	ctx->Px = (float *)calloc(NOS, sizeof(float));
	ctx->Py = (float *)calloc(NOS, sizeof(float));
	ctx->Pz = (float *)calloc(NOS, sizeof(float));	// <-- + 12 Mega Byte

	ctx->BPx = (float *)calloc(ctx->NOB, sizeof(float));
	ctx->BPy = (float *)calloc(ctx->NOB, sizeof(float));
	ctx->BPz = (float *)calloc(ctx->NOB, sizeof(float));	// <-- + 12 Mega Byte

	ctx->Heffx = (double *)calloc(NOS, sizeof(double));
	ctx->Heffy = (double *)calloc(NOS, sizeof(double));
	ctx->Heffz = (double *)calloc(NOS, sizeof(double));// <-- + 12 Mega Byte

	ctx->Etot = (double *)calloc(NOS, sizeof(double));// <-- + 4 Mega Byte
	ctx->Etot0= (double *)calloc(NOS, sizeof(double));// <-- + 4 Mega Byte
	// For 100x100x100x10^6 spins total allocated memory is about 6*12 = 72 Mega Byte
	// and possibly may reach up to 100 Mb in total with all other variables.

	ctx->Proj = (bool *)calloc(NOS, sizeof(bool));

}

void RestartCalcThreads(magnoom_ctx *ctx, pthread_t * thread_id, calc_thread_arg * thread_args){
		for (int i=0; i<THREADS_NUMBER; i++){
		thread_args[i].index = i;
		thread_args[i].ctx = ctx;
		if ( pthread_create(&thread_id[i], NULL, CALC_THREAD, (void *)&thread_args[i]) ) {
			fprintf(stderr, "Error in creating CALC_THREAD thread\n");
		}
	}
}

/*************************************************************************/
/*                        Program Main Thread                            */
/*************************************************************************/

int 
main (int argc, char **argv){
	if (!magnoom_ctx_init(&mag_ctx)) {
		fprintf(stderr, "Unable to initialize Magnoom data.\n");
		return 1;
	}
		
	for (int i=0; i<THREADS_NUMBER; i++){
		char name[10]; 
		snprintf(name,10,"inDoor%d\n",i);
		//printf("%s\n",name);
		if ( (mag_ctx.sem_in[i] = sem_open(name, O_CREAT, 0644, 0)) == SEM_FAILED ) perror("sem_open");
		snprintf(name,10,"outDoor%d\n",i);
		//printf("%s\n",name); 
		if ( (mag_ctx.sem_out[i] = sem_open(name, O_CREAT, 0644, 0)) == SEM_FAILED ) perror("sem_open");
		//int value;		
		//sem_getvalue(sem_in[i], &value); //Function not implemented on Mac OS X!!!
		mag_ctx.Max_torque[i]=0;
	} 


	////////////////////////////////////////////////
	srand ( time(NULL) );//init random number seed//
	////////////////////////////////////////////////

	////////////////////////////////////////////////
	readConfigFile(&mag_ctx);
	////////////////////////////////////////////////



	GetShells(mag_ctx.abc, mag_ctx.Block, mag_ctx.AtomsPerBlock, mag_ctx.ShellNumber, mag_ctx.RadiusOfShell);
	for(int i=0;i<mag_ctx.ShellNumber;i++) printf("R[%d]=%f\n",i,mag_ctx.RadiusOfShell[i] );
	mag_ctx.NeighborPairs = GetNeighborsNumber(mag_ctx.abc, mag_ctx.Block, mag_ctx.AtomsPerBlock, mag_ctx.ShellNumber, mag_ctx.RadiusOfShell, mag_ctx.NeighborsPerAtom);
	// Allocate arrays for neighbours map:
	mag_ctx.AIdxBlock = (int *)calloc(mag_ctx.NeighborPairs, sizeof(int));// index of the atom within the block
	mag_ctx.NIdxBlock = (int *)calloc(mag_ctx.NeighborPairs, sizeof(int));// index of the neighbour within the block
	mag_ctx.NIdxGridA = (int *)calloc(mag_ctx.NeighborPairs, sizeof(int));// index of the neighbour within the block
	mag_ctx.NIdxGridB = (int *)calloc(mag_ctx.NeighborPairs, sizeof(int));// index of the relative position of the block of the neighbour in the greed along tr. vect. a
	mag_ctx.NIdxGridC = (int *)calloc(mag_ctx.NeighborPairs, sizeof(int));// index of the relative position of the block of the neighbour in the greed along tr. vect. b
	mag_ctx.SIdx      = (int *)calloc(mag_ctx.NeighborPairs, sizeof(int));// index of the shell corresponding to this pair

	CreateMapOfNeighbors( mag_ctx.abc, mag_ctx.Block, mag_ctx.AtomsPerBlock, mag_ctx.ShellNumber, mag_ctx.RadiusOfShell, mag_ctx.AIdxBlock, mag_ctx.NIdxBlock, mag_ctx.NIdxGridA, mag_ctx.NIdxGridB, mag_ctx.NIdxGridC, mag_ctx.SIdx);
	mag_ctx.Jexc = (float *)calloc(mag_ctx.NeighborPairs, sizeof(float));
	mag_ctx.Bexc = (float *)calloc(mag_ctx.NeighborPairs, sizeof(float));
	mag_ctx.Dexc = (float *)calloc(mag_ctx.NeighborPairs, sizeof(float));
	mag_ctx.VDMx = (float *)calloc(mag_ctx.NeighborPairs, sizeof(float));
	mag_ctx.VDMy = (float *)calloc(mag_ctx.NeighborPairs, sizeof(float));
	mag_ctx.VDMz = (float *)calloc(mag_ctx.NeighborPairs, sizeof(float));
	SetExch1( mag_ctx.abc, mag_ctx.Block, mag_ctx.NeighborPairs, mag_ctx.Jij, mag_ctx.Bij, mag_ctx.Dij, mag_ctx.AIdxBlock, mag_ctx.NIdxBlock, mag_ctx.NIdxGridA, mag_ctx.NIdxGridB, mag_ctx.NIdxGridC, mag_ctx.SIdx, 
	mag_ctx.Jexc, mag_ctx.Bexc, mag_ctx.Dexc, mag_ctx.VDMx, mag_ctx.VDMy, mag_ctx.VDMz);
	// SetExchMarkus( mag_ctx.abc, Block, NeighborPairs, AIdxBlock, NIdxBlock, NIdxGridA, NIdxGridB, NIdxGridC, SIdx, Jexc, Bexc, Dexc, mag_ctx.VDMx, mag_ctx.VDMy, mag_ctx.VDMz);
	// SetExchChAch( mag_ctx.abc, Block, NeighborPairs, AIdxBlock, NIdxBlock, NIdxGridA, NIdxGridB, NIdxGridC, SIdx, Jexc, Bexc, Dexc, mag_ctx.VDMx, mag_ctx.VDMy, mag_ctx.VDMz);
	// SetExchMariya( mag_ctx.abc, Block, NeighborPairs, mag_ctx.Jij, mag_ctx.Bij, mag_ctx.Dij, AIdxBlock, NIdxBlock, NIdxGridA, NIdxGridB, NIdxGridC, SIdx, Jexc, Bexc, Dexc, mag_ctx.VDMx, mag_ctx.VDMy, mag_ctx.VDMz);
	for(int i=0;i<mag_ctx.AtomsPerBlock;i++) printf("Neighbors Per Atom[%d]=%d\n",i,mag_ctx.NeighborsPerAtom[i] );
	printf("AtomsPerBlock =%2d\n", mag_ctx.AtomsPerBlock);
	printf("  ShellNumber =%2d\n", mag_ctx.ShellNumber);
	printf("NeighborPairs =%2d\n", mag_ctx.NeighborPairs);
	float AXYZ[3] = {0.0, 0.0, 0.0};//atom position
	float NXYZ[3] = {0.0, 0.0, 0.0};//neighbor position

	// Neighbours map to console:	
	printf("  N\t| Ia\t| In\t| Jn\t| Kn\t| Ln\t| Sn\t| Ax\t| Ay\t| Az\t| Nx\t| Ny\t| Nz\t|Dist\t| Jij\t| Dij\t| D_x\t| D_y\t| D_z\t|\n");
	for(int i=0;i<mag_ctx.NeighborPairs;i++) {
		printf("%3d  \t|",            i);
		printf("%3d  \t|", mag_ctx.AIdxBlock[i]); //I'
		printf("%3d  \t|", mag_ctx.NIdxBlock[i]); //I
		printf("%3d  \t|", mag_ctx.NIdxGridA[i]); //J
		printf("%3d  \t|", mag_ctx.NIdxGridB[i]); //K
		printf("%3d  \t|", mag_ctx.NIdxGridC[i]); //L
		printf("%3d  \t|",      mag_ctx.SIdx[i]); //S
		GetPosition( mag_ctx.abc, mag_ctx.Block, mag_ctx.AIdxBlock[i], 0, 0, 0, AXYZ);
		printf("%2.3f\t|%2.3f\t|%2.3f\t|",AXYZ[0],AXYZ[1],AXYZ[2]);
		GetPosition( mag_ctx.abc, mag_ctx.Block, mag_ctx.NIdxBlock[i], mag_ctx.NIdxGridA[i], mag_ctx.NIdxGridB[i], mag_ctx.NIdxGridC[i], NXYZ);
		printf("%2.3f\t|%2.3f\t|%2.3f\t|",NXYZ[0],NXYZ[1],NXYZ[2]);
		AXYZ[0]-=NXYZ[0];
		AXYZ[1]-=NXYZ[1];
		AXYZ[2]-=NXYZ[2];
		printf("%2.3f\t|",  sqrt(AXYZ[0]*AXYZ[0]+AXYZ[1]*AXYZ[1]+AXYZ[2]*AXYZ[2])); //Distance
		printf("%2.3f\t|",  mag_ctx.Jexc[i]); //mag_ctx.Jij 
		printf("%2.3f\t|",  mag_ctx.Dexc[i]); //mag_ctx.Jij 
		printf("%2.3f\t|",  mag_ctx.VDMx[i]); //DMI x
		printf("%2.3f\t|",  mag_ctx.VDMy[i]); //DMI y
		printf("%2.3f\t|\n",mag_ctx.VDMz[i]); //DMI z
	}

	// Open output file:
	mag_ctx.outFile = fopen ("table.csv","w");
	if (mag_ctx.outFile!=NULL) {fputs ("iter,time,Mx,My,Mz,E_tot,\n",mag_ctx.outFile);}

	// Memory allocation:
	mag_ctx.RHue = (float *)calloc(360, sizeof(float));
	mag_ctx.GHue = (float *)calloc(360, sizeof(float));
	mag_ctx.BHue = (float *)calloc(360, sizeof(float));

	ReallocateMemoryForSpins(&mag_ctx, mag_ctx.NOS);
	ReallocateMemoryForAllOther(&mag_ctx, mag_ctx.NOS);
	ReallocateMemoryForImages(&mag_ctx, mag_ctx.num_images, mag_ctx.NOS);


	pthread_mutex_init(&mag_ctx.culc_mutex,0);
	pthread_mutex_init(&mag_ctx.show_mutex,0);
	pthread_t thread_id[THREADS_NUMBER];
	calc_thread_arg thread_args[THREADS_NUMBER];
	RestartCalcThreads(&mag_ctx, thread_id, thread_args);

	GetBox(&mag_ctx, mag_ctx.abc, mag_ctx.uABC, mag_ctx.Box);
	UpdateSpinPositions(&mag_ctx, mag_ctx.abc, mag_ctx.uABC, mag_ctx.Block, mag_ctx.AtomsPerBlock, mag_ctx.Box, mag_ctx.Px, mag_ctx.Py, mag_ctx.Pz);
	UpdateKind(&mag_ctx, mag_ctx.Kind, mag_ctx.Px, mag_ctx.Py, mag_ctx.Pz, mag_ctx.NOS, mag_ctx.NOSK);
	InitSpinComponents( &mag_ctx, mag_ctx.Px, mag_ctx.Py, mag_ctx.Pz, mag_ctx.Sx, mag_ctx.Sy, mag_ctx.Sz, 0);
	for (int i=0;i<mag_ctx.NOS;i++) { mag_ctx.bSx[i]=mag_ctx.Sx[i]; mag_ctx.bSy[i]=mag_ctx.Sy[i]; mag_ctx.bSz[i]=mag_ctx.Sz[i];}

    // Set OpenGL context initial state.
	setupOpenGL(&mag_ctx);
    //  Allocate memory for vetices, mag_ctx.normals, mag_ctx.colors and indicies array used in drawing subrutines
	ReallocateArrayDrawing(&mag_ctx);
	// Fill array for prototype (arrow or cane) array 
	UpdatePrototypeVerNorInd(&mag_ctx, mag_ctx.vertexProto, mag_ctx.normalProto, mag_ctx.indicesProto, mag_ctx.arrowFaces, mag_ctx.WhichVectorMode,0);
	// Fill big array for indecies for all arrows, cans, cones or boxes 
	UpdateIndices(&mag_ctx, mag_ctx.indicesProto , mag_ctx.IdNumProto, mag_ctx.indices, mag_ctx.IdNum, mag_ctx.VCNumProto); 
	UpdateVerticesNormalsColors(&mag_ctx, mag_ctx.vertexProto, mag_ctx.normalProto, mag_ctx.VCNumProto, mag_ctx.vertices, mag_ctx.normals, mag_ctx.colors, mag_ctx.VCNum, mag_ctx.Px, mag_ctx.Py, mag_ctx.Pz, mag_ctx.bSx, mag_ctx.bSy, mag_ctx.bSz, mag_ctx.WhichVectorMode);
	CreateNewVBO(&mag_ctx);
	UpdateVBO(&mag_ctx, &mag_ctx.vboIdV, &mag_ctx.vboIdN, &mag_ctx.vboIdC, &mag_ctx.iboIdI, mag_ctx.vertices, mag_ctx.normals, mag_ctx.colors, mag_ctx.indices);

    mag_ctx.VHf[0]=sin(PI*mag_ctx.VHtheta/180)*cos(PI*mag_ctx.VHphi/180);
	mag_ctx.VHf[1]=sin(PI*mag_ctx.VHtheta/180)*sin(PI*mag_ctx.VHphi/180);
	mag_ctx.VHf[2]=cos(PI*mag_ctx.VHtheta/180);

	ReallocateArrayDrawing_H(&mag_ctx);
    UpdatePrototypeVerNorInd(&mag_ctx, mag_ctx.vertexProto_H, mag_ctx.normalProto_H, mag_ctx.indices_H, mag_ctx.arrowFaces_H, ARROW1,1);
    UpdateVerticesNormalsColors_H(&mag_ctx, mag_ctx.vertexProto_H, mag_ctx.normalProto_H, mag_ctx.VCNum_H, mag_ctx.vertices_H, mag_ctx.normals_H, mag_ctx.colors_H, mag_ctx.Box[0][0]*0.6, mag_ctx.Box[1][1]*0.6, mag_ctx.Box[2][2]*0.6, mag_ctx.VHf[0], mag_ctx.VHf[1], mag_ctx.VHf[2]);
    CreateNewVBO_H(&mag_ctx);
    UpdateVBO_H(&mag_ctx, &mag_ctx.vboIdV_H, &mag_ctx.vboIdN_H, &mag_ctx.vboIdC_H, &mag_ctx.iboIdI_H, mag_ctx.vertices_H, mag_ctx.normals_H, mag_ctx.colors_H, mag_ctx.indices_H);

	ReallocateArrayDrawing_BOX(&mag_ctx);
	UpdateVerticesNormalsColors_BOX(&mag_ctx, mag_ctx.vertices_BOX, mag_ctx.normals_BOX, mag_ctx.colors_BOX, mag_ctx.indices_BOX, mag_ctx.Box);
	CreateNewVBO_BOX(&mag_ctx);
	UpdateVBO_BOX(&mag_ctx, &mag_ctx.vboIdV_BOX, &mag_ctx.vboIdN_BOX, &mag_ctx.vboIdC_BOX, &mag_ctx.iboIdI_BOX, mag_ctx.vertices_BOX, mag_ctx.normals_BOX, mag_ctx.colors_BOX, mag_ctx.indices_BOX);

	ReallocateArrayDrawing_BASIS(&mag_ctx);
	UpdateVerticesNormalsColors_BASIS(&mag_ctx, mag_ctx.vertices_BASIS, mag_ctx.normals_BASIS, mag_ctx.colors_BASIS, mag_ctx.indices_BASIS, mag_ctx.Box);
	CreateNewVBO_BASIS(&mag_ctx);
	UpdateVBO_BASIS(&mag_ctx, &mag_ctx.vboIdV_BASIS, &mag_ctx.vboIdN_BASIS, &mag_ctx.vboIdC_BASIS, &mag_ctx.iboIdI_BASIS, mag_ctx.vertices_BASIS, mag_ctx.normals_BASIS, mag_ctx.colors_BASIS, mag_ctx.indices_BASIS);

	ReallocateArrayDrawing_PBC(&mag_ctx);
	UpdateVerticesNormalsColors_PBC(&mag_ctx, 0, mag_ctx.vertices_PBC_A, mag_ctx.normals_PBC_A, mag_ctx.colors_PBC_A, mag_ctx.indices_PBC_A, mag_ctx.Box);
	UpdateVerticesNormalsColors_PBC(&mag_ctx, 1, mag_ctx.vertices_PBC_B, mag_ctx.normals_PBC_B, mag_ctx.colors_PBC_B, mag_ctx.indices_PBC_B, mag_ctx.Box);
	UpdateVerticesNormalsColors_PBC(&mag_ctx, 2, mag_ctx.vertices_PBC_C, mag_ctx.normals_PBC_C, mag_ctx.colors_PBC_C, mag_ctx.indices_PBC_C, mag_ctx.Box);  
	CreateNewVBO_PBC(&mag_ctx);
    UpdateVBO_PBC(&mag_ctx, &mag_ctx.vboIdV_PBC_A, &mag_ctx.vboIdN_PBC_A, &mag_ctx.vboIdC_PBC_A, &mag_ctx.iboIdI_PBC_A, mag_ctx.vertices_PBC_A, mag_ctx.normals_PBC_A, mag_ctx.colors_PBC_A, mag_ctx.indices_PBC_A);
    UpdateVBO_PBC(&mag_ctx, &mag_ctx.vboIdV_PBC_B, &mag_ctx.vboIdN_PBC_B, &mag_ctx.vboIdC_PBC_B, &mag_ctx.iboIdI_PBC_B, mag_ctx.vertices_PBC_B, mag_ctx.normals_PBC_B, mag_ctx.colors_PBC_B, mag_ctx.indices_PBC_B);
    UpdateVBO_PBC(&mag_ctx, &mag_ctx.vboIdV_PBC_C, &mag_ctx.vboIdN_PBC_C, &mag_ctx.vboIdC_PBC_C, &mag_ctx.iboIdI_PBC_C, mag_ctx.vertices_PBC_C, mag_ctx.normals_PBC_C, mag_ctx.colors_PBC_C, mag_ctx.indices_PBC_C);

	// Explicit GLFW loop, matching the AntTweakBar Legacy examples.
	while (!glfwWindowShouldClose(mag_ctx.MainWindow)) {
		idle(&mag_ctx);
		Display(&mag_ctx);
		glfwSwapBuffers(MainWindow);
		glfwPollEvents();
	}

	mag_ctx.EngineShutdownRequested = true;
	mag_ctx.ENGINE_MUTEX = DO_IT;
	for (int i = 0; i < THREADS_NUMBER; ++i) pthread_join(thread_id[i], NULL);
	TwTerminate();
	glfwTerminate();

	// Deallocate memory before closing the program:
	free(mag_ctx.AIdxBlock);
	free(mag_ctx.NIdxBlock);
	free(mag_ctx.NIdxGridA);
	free(mag_ctx.NIdxGridB);
	free(mag_ctx.NIdxGridC);
	free(mag_ctx.SIdx);

	free(mag_ctx.Jexc);  			free(mag_ctx.Bexc);  			free(mag_ctx.Dexc);
	free(mag_ctx.VDMx);  			free(mag_ctx.VDMy);  			free(mag_ctx.VDMz);

	free(mag_ctx.Sx);    			free(mag_ctx.Sy);    			free(mag_ctx.Sz);
	free(mag_ctx.tSx);   			free(mag_ctx.tSy);   			free(mag_ctx.tSz);
	free(mag_ctx.bSx);   			free(mag_ctx.bSy);   			free(mag_ctx.bSz);
	free(mag_ctx.t2Sx);  			free(mag_ctx.t2Sy);  			free(mag_ctx.t2Sz);
	free(mag_ctx.t3Sx);  			free(mag_ctx.t3Sy);  			free(mag_ctx.t3Sz);
	free(mag_ctx.Heffx); 			free(mag_ctx.Heffy); 			free(mag_ctx.Heffz);
	free(mag_ctx.RNx);   			free(mag_ctx.RNy);   			free(mag_ctx.RNz);
	free(mag_ctx.Px);    			free(mag_ctx.Py);    			free(mag_ctx.Pz);
	free(mag_ctx.BPx);   			free(mag_ctx.BPy);   			free(mag_ctx.BPz);
	free(mag_ctx.RHue);  			free(mag_ctx.GHue);  			free(mag_ctx.BHue);
	free(mag_ctx.vertices);			free(mag_ctx.vertices_H); 		free(mag_ctx.vertices_BOX); 	free(mag_ctx.vertices_PBC_A);
	free(mag_ctx.normals);			free(mag_ctx.normals_H); 		free(mag_ctx.normals_BOX); 		free(mag_ctx.normals_PBC_A);
	free(mag_ctx.colors);			free(mag_ctx.colors_H); 		free(mag_ctx.colors_BOX); 		free(mag_ctx.colors_PBC_A);
	free(mag_ctx.indices);			free(mag_ctx.indices_H); 		free(mag_ctx.indices_BOX); 		free(mag_ctx.indices_PBC_A);
	free(mag_ctx.vertexProto);		free(mag_ctx.vertexProto_H);
	free(mag_ctx.normalProto);		free(mag_ctx.normalProto_H);
	free(mag_ctx.indicesProto);		free(mag_ctx.indicesProto_H);

	free(mag_ctx.RadiusOfShell);
	free(mag_ctx.NeighborsPerAtom);

	free(mag_ctx.Proj);

	fclose (mag_ctx.outFile);

	return 0;    
}
