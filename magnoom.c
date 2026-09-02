//*****************************************************************************
//
//  Project Name : Magnoom
//  Author       : Nikolai S. Kiselev
//  Created      : April 2016
//  Modified     : October 2016
//  
//  Build with the repository's nob.c build system (see README.md).
#if !defined(_WIN32) && !defined(_POSIX_C_SOURCE)
#define _POSIX_C_SOURCE 200809L
#endif

#include <glad/glad.h>
#include <AntTweakBar.h>
// This app never calls a GLU function, so skip glfw.h's unconditional
// GL/glu.h include (not guaranteed to be installed).
#define GLFW_NO_GLU
#include <GL/glfw.h>

#include <math.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <ctype.h>
#include <stdarg.h>
#include <stddef.h>

#if defined(_WIN32)
	#include <windows.h>
	#include <process.h>
	#define snprintf _snprintf			
	#define SEM_VALUE_MAX 32767
	#define SEM_FAILED NULL
	#define pthread_t HANDLE
	#define semaphore_ref HANDLE
	#define pthread_mutex_t CRITICAL_SECTION
	#define sem_open(name, flag, mode, value) CreateSemaphore(NULL,value,SEM_VALUE_MAX,name)
	#define sem_post(sem) ReleaseSemaphore(sem,1,NULL)
	#define sem_wait(sem) WaitForSingleObject(sem,INFINITE)
	#define sem_trywait(sem) WaitForSingleObject(sem,0)
	#define sem_close(sem) CloseHandle(sem)
	#define pthread_create(th_ref, attr_ref, name, arg_ref) ((  ( *(th_ref) = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)(*(name)), arg_ref, 0, NULL) ) == NULL ? -1 : 0  ))
	#define pthread_mutex_init(mutex_ptr,num) InitializeCriticalSection(mutex_ptr)
	#define pthread_mutex_lock(mutex_ptr) EnterCriticalSection(mutex_ptr)
	#define pthread_mutex_unlock(mutex_ptr) LeaveCriticalSection(mutex_ptr)
	#define pthread_mutex_destroy(mutex_ptr) DeleteCriticalSection(mutex_ptr)
	#define pthread_join(thread,retval) (WaitForSingleObject((thread), INFINITE), CloseHandle(thread), 0)
#else
	#include <unistd.h>
	#include <pthread.h>
	#include <semaphore.h>		
	#define semaphore_ref sem_t*
#endif

#if defined(__APPLE__)
	#include <mach-o/dyld.h>
#endif

static const char SOFTWARE_NAME[] = "Magnoom";
static const char SOFTWARE_VERSION[] = "1.06";

#define ABS(x) ((x)<0?-(x):(x))

enum engine_mutex_flags{DO_IT,WAIT,STOP_REQUESTED};
enum data_mutex_flags{WAIT_DATA,TAKE_DATA};

#define THREADS_NUMBER 3
#define MAX_ATOMS_PER_BLOCK 5
#define MAX_SHELL_NUM 16
#define MAGNOOM_PATH_CAPACITY 4096
#define MAGNOOM_MODAL_MESSAGE_CAPACITY 512

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "vendor/stb/stb_image_write.h"

/*****************************************************************************/
/* Constants formerly in magnoom.h                                          */
/*****************************************************************************/
#define PI      3.14159265359
#define TPI     (2.0 * PI)
#define iTPI    (1.0 / TPI)
#define D2R     (PI / 180.0)
#define R2D     (180.0 / PI)
#define HIDDEN_VECTOR_SCALE 24000000
// active mouse buttons (or them together):
#define LEFT   4
#define MIDDLE 2
#define RIGHT  1
// number of slots for camera positions to save in memory
#define CAMERA_POSITION_SLOTS 5

/*****************************************************************************/
/* Enum types needed by magnoom_ctx fields (formerly in magnoom.h/visualization.c) */
/*****************************************************************************/
typedef enum    {A_AXIS, B_AXIS, C_AXIS, FILTER} enSliceMode;
typedef enum    {BEXT_AC_SIN, BEXT_AC_GAUSSIAN, BEXT_AC_SINC, BEXT_AC_CIRCULAR} enBextACWaveform;
typedef enum    {FILE_FORMAT_UNKNOWN = 0, FILE_FORMAT_CSV, FILE_FORMAT_OVF,
                 FILE_FORMAT_VTK, FILE_FORMAT_BIN, FILE_FORMAT_PNG,
                 FILE_FORMAT_COUNT} FileFormatEnum;
enum            IntegrationScheme{HEUN,SIB,RK2,RK4,RELAX};
enum            Average_mode{ALONG_A,ALONG_B, ALONG_C, ALONG_0};

typedef enum    {ORTHO, PERSP} enProjections;
typedef enum    {RND, HOMO, SKYRM1, SKYRM2, SKYRM3, BOBBER_T, BOBBER_B, BOBBER_L, BOBBER_L_T, BOBBER_L_B,
                 HOPFION1, SPIRAL, SKYRMION_L, GLOBULA, MultyQ, NORM} enIniState;
typedef enum    {DEFAULT_G, CILINDER_G, SPHERE_G} enGeom;
typedef enum    {WHITE, BLACK, RED, GREEN, BLUE, MANUAL} enColors;
typedef enum    {ARROW1, CONE1, CANE, uPOINT, BOX1} enVectorMode;
typedef enum    {ANISOTROPY_GLOBAL, ANISOTROPY_INDIVIDUAL} AnisotropyMode;
typedef enum    {LEGACY_ASCII_VTK, LEGACY_BINARY_VTK} VtkFormatMode;
typedef enum    {ANISOTROPY_RECORD_K2, ANISOTROPY_RECORD_K4, ANISOTROPY_RECORD_QUATERNION} AnisotropyRecordKind;
typedef enum    {ANISOTROPY_COMPONENT_K2, ANISOTROPY_COMPONENT_K4} AnisotropyComponentKind;
enum { ANISOTROPY_K2_COMPONENT_COUNT = 6, ANISOTROPY_K4_COMPONENT_COUNT = 15 };

typedef struct AnisotropyTensor {
	double K2[3][3];
	double K4[3][3][3][3];
} AnisotropyTensor;

/* Defined in anisotropy.c (included near the end of this file); forward
 * declared here so magnoom_ctx_init() can set compiled-in tensor defaults. */
bool k2_set(double K2[3][3], int i, int j, double value);
bool k4_set(double K4[3][3][3][3], int i, int j, int k, int l, double value);

#define MAX_ANISOTROPY_CONFIG_RECORDS 4096
typedef struct AnisotropyConfigRecord {
	AnisotropyRecordKind kind;
	int atom;
	int index[4];
	double value;
	int line;
} AnisotropyConfigRecord;

typedef struct AnisotropyComponentControl {
	struct magnoom_ctx *ctx;
	AnisotropyComponentKind kind;
	int component;
} AnisotropyComponentControl;
typedef enum    {LIGHT_OFF, LIGHT_FIXED, LIGHT_ADAPTIVE} enLightingMode;

/*****************************************************************************/
/* vbo_mesh: descriptor for GPU buffer groups used in fixed-function drawing */
/*****************************************************************************/
typedef struct {
	GLuint          vertex_buffer;
	GLuint          normal_buffer;
	GLuint          color_buffer;
	GLuint          index_buffer;
	size_t          component_count;
	size_t          component_capacity;
	size_t          index_count;
	size_t          index_capacity;
	GLenum          usage;
	int             uses_normals;
} vbo_mesh;

enum {
	VBO_UPLOAD_VERTICES = 1 << 0,
	VBO_UPLOAD_NORMALS  = 1 << 1,
	VBO_UPLOAD_COLORS   = 1 << 2,
	VBO_UPLOAD_INDICES  = 1 << 3,
	VBO_UPLOAD_ATTRIBUTES = VBO_UPLOAD_VERTICES | VBO_UPLOAD_NORMALS | VBO_UPLOAD_COLORS,
	VBO_UPLOAD_ALL = VBO_UPLOAD_ATTRIBUTES | VBO_UPLOAD_INDICES
};

/*****************************************************************************/
/* NeighborMap: flattened template of neighbor-pair bonds in the basic       */
/* block, grouped into coordination shells distinguished by distance and    */
/* local symmetry (see magnoom_ctx_build_neighbor_map in lattice_geometry.c).*/
/* Purely topological -- exchange constants (Jij/Bij/Dij) and per-pair DMI  */
/* directions (VDMX/Y/Z) live directly on magnoom_ctx, indexed by ShellIdx. */
/*****************************************************************************/
typedef struct NeighborMap {
	int             ShellCount;                /* realized shell count, <= MAX_SHELL_NUM */
	float           ShellRadius[MAX_SHELL_NUM]; /* radius of shell s */

	int             PairCount;                  /* total flattened neighbor-pair entries */
	int*            AIdxBlock;                  /* [PairCount] central atom index within the block */
	int*            NIdxBlock;                  /* [PairCount] neighbor atom index within the block */
	int*            NIdxGridA;                  /* [PairCount] neighbor block offset along vector a */
	int*            NIdxGridB;                  /* [PairCount] neighbor block offset along vector b */
	int*            NIdxGridC;                  /* [PairCount] neighbor block offset along vector c */
	int*            ShellIdx;                   /* [PairCount] shell index (0..ShellCount-1) of this pair */
} NeighborMap;

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

	/* Cartesian components of spins, interleaved (S[3*n+0..2] via VEC_X/Y/Z) */
	double*         S;      /* main spin array, length 3*NOS                */
	double*         bS;     /* backup/secondary spin array, length 3*NOS    */
	double*         tS;     /* RK-integrator temporary, length 3*NOS        */
	double*         t2S;    /* RK-integrator temporary, length 3*NOS        */
	double*         t3S;    /* RK-integrator temporary, length 3*NOS        */
	double*         Image;  /* flat snapshot block, IMAGE_COMPONENT-indexed */
	double*         dImage; /* flat derivative-snapshot block, same layout  */
	float*          IsoLineX;
	float*          IsoLineY;
	float*          IsoLineZ;
	float*          NoiseX;
	float*          NoiseY;
	float*          NoiseZ;
	float*          PosX;
	float*          PosY;
	float*          PosZ;
	float*          BlockPosX;
	float*          BlockPosY;
	float*          BlockPosZ;
	int*            Kind;
	double*         HeffX;
	double*         HeffY;
	double*         HeffZ;
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
	double          RecordingBuffer[10][1500];
	int             recordsCounter;
	int             NumOfPoints;

	/* Color scheme */
	int             HueMapRGB[6];
	int             HueMapRYGB[6];
	int             HueMap[6];
	float*          RHue;
	float*          GHue;
	float*          BHue;

	/* Neighbor map */
	NeighborMap     Neighbors;
	float           Box[3][3];
	float*          Jexc;
	float*          Bexc;
	float*          Dexc;
	float*          VDMX;
	float*          VDMY;
	float*          VDMZ;

	/* Heisenberg / biquadratic / DMI exchange (indexed by shell) */
	float           Jij[MAX_SHELL_NUM];
	float           Bij[MAX_SHELL_NUM];
	float           Dij[MAX_SHELL_NUM];

	/* Magnetocrystalline anisotropy */
	AnisotropyTensor anisotropy_local[MAX_ATOMS_PER_BLOCK];
	AnisotropyTensor anisotropy_global[MAX_ATOMS_PER_BLOCK];
	double          anisotropy_quaternion[MAX_ATOMS_PER_BLOCK][4];
	float           anisotropy_rotation_axis[3];
	float           anisotropy_rotation_angle;
	AnisotropyMode  anisotropy_mode;
	AnisotropyConfigRecord anisotropy_config_records[MAX_ANISOTROPY_CONFIG_RECORDS];
	int             anisotropy_config_record_count;
	int             anisotropy_selected_atom;
	AnisotropyComponentControl anisotropy_component_controls[
		ANISOTROPY_K2_COMPONENT_COUNT + ANISOTROPY_K4_COMPONENT_COUNT];
	double          anisotropy_map_theta_step; /* polar-angle sampling step, radians */
	double          anisotropy_map_phi_step;   /* azimuthal-angle sampling step, radians */
	VtkFormatMode   anisotropy_map_vtk_format;

	/* External magnetic field: static (DC) and time-dependent (AC) components */
	float           BextDCDirection[3];
	float           BextDCTheta;
	float           BextDCPhi;
	float           BextDCMagnitude;
	float           BextDC[3];
	float           BextACDirection[3];
	float           BextACAmplitude;
	double          BextAC[3];
	double          BextACOmega;
	double          BextACPeriod;
	float           BextACPulseWidth;
	float           BextACTimeOffset;
	float           BextACScalar;
	enBextACWaveform       BextACWaveform;

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
	unsigned int    ITERATION;
	unsigned int    Max_Numb_Iteration;
	int             Record;
	int             BextACEnabled;
	int             BextACModeRecording;

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
	char            BufferString[800];
	double          outputEtotal;
	double          outputMtotal[3];
	int             SleepTime;
	char            ShortBuffer[200];
	char            input_directory[MAGNOOM_PATH_CAPACITY];
	char            output_directory[MAGNOOM_PATH_CAPACITY];
	char            inputfilename[MAGNOOM_PATH_CAPACITY];
	char            outputfilename[MAGNOOM_PATH_CAPACITY];
	char            record_path[MAGNOOM_PATH_CAPACITY];

	/* Geometry */
	float           abc[3][3];
	float           (*Block)[3];
	int             uABC[3];
	int             ShellNumber;
	int             AtomsPerBlock;
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
	pthread_mutex_t record_mutex;
	int             EngineRunState;
	int             DataTransferState;
	bool            EngineIdle;
	volatile bool   EngineShutdown;
	volatile bool   EngineShutdownRequested;
	semaphore_ref   sem_in[THREADS_NUMBER];
	semaphore_ref   sem_out[THREADS_NUMBER];

	/* visualization.c: window / GLFW / display state */
	bool            WindowShouldClose;
	const char*     WINDOWTITLE;
	int             window_width;
	int             window_height;
	float           asp_rat;
	float           asp_rat_inv;
	double          ContentScaleX;
	double          ContentScaleY;
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
	float           CameraPosition[CAMERA_POSITION_SLOTS][7];
	float           axis[3];
	float           angle;

	/* OpenGL material / lighting */
	GLfloat         ambient[4];
	GLfloat         diffuse[4];
	GLfloat         specular[4];
	GLfloat         shininess;
	float           LightMultiplier;
	float           LightDirection[3];
	enLightingMode  WhichLightingMode;

	/* initial-state generation parameters */
	enIniState      WhichInitialState;
	enGeom          WhichGeometry;
	float           chSizeG;
	float           chSize;
	float           chDir[3];
	float           RotateAllSpins;

	/* color scheme (visualization) */
	enColors        WhichBackgroundColor;
	int             temp_color[3];
	GLfloat         BackgroundColors[6][3];

	/* vector-drawing mode / appearance */
	enVectorMode    WhichVectorMode;
	float           Scale;
	float           Pivot;
	float           WireWidth;
	float           Scale_BextDC;
	int             arrowFaces;
	int             arrowFaces_BextDC;
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
	TwBar*          BextAC_bar;
	TwBar*          anisotropy_bar;
	TwBar*          info_bar;
	TwBar*          modal_bar;
	bool            modal_open_requested;
	bool            modal_close_requested;
	bool            modal_is_warning; /* false = error dialog (red), true = warning dialog (dark magenta) */
	char            modal_message[MAGNOOM_MODAL_MESSAGE_CAPACITY];

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
	GLfloat*        vertexProto_BextDC;
	GLfloat*        normalProto_BextDC;
	GLuint*         indicesProto_BextDC;
	GLfloat*        vertices;
	GLfloat*        normals;
	GLfloat*        colors;
	GLuint*         indices;
	GLfloat*        vertices_BextDC;
	GLfloat*        normals_BextDC;
	GLfloat*        colors_BextDC;
	GLuint*         indices_BextDC;
	GLfloat*        vertices_BOX;
	GLfloat*        normals_BOX;
	GLfloat*        colors_BOX;
	GLuint*         indices_BOX;
	GLfloat*        vertices_BASIS;
	GLfloat*        normals_BASIS;
	GLfloat*        colors_BASIS;
	GLuint*         indices_BASIS;
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

	/* geometry-sizing counters */
	int             ElNumProto;
	int             IdNumProto;
	int             VCNumProto;
	int             ElNum;
	int             IdNum;
	int             VCNum;
	int             IdNum_BextDC;
	int             VCNum_BextDC;
	int             IdNum_BOX;
	int             VCNum_BOX;
	int             IdNum_BASIS;
	int             VCNum_BASIS;
	int             IdNum_PBC;
	int             VCNum_PBC;

	/* vbo_mesh descriptors replacing parallel GPU-ID and count fields */
	vbo_mesh        spin_mesh;
	vbo_mesh        BextDC_mesh;
	vbo_mesh        box_mesh;
	vbo_mesh        basis_mesh;
	vbo_mesh        pbc_mesh[3];
} magnoom_ctx;

/* Per-thread argument for CALC_THREAD (solvers.c): carries both the thread's */
/* index and ctx, since pthread's start-routine signature has room for only  */
/* one void* argument. */
typedef struct { int index; magnoom_ctx *ctx; } calc_thread_arg;

typedef struct magnoom_bin_spin {
	unsigned short int t;
	unsigned short int f;
} magnoom_bin_spin;

static const char *fileFormatNames[FILE_FORMAT_COUNT] = {
	"unknown", "CSV", "OVF", "VTK", "BIN", "PNG"
};

static bool magnoom_copy_path(char *destination, size_t capacity, const char *source)
{
	size_t length;
	if (destination == NULL || capacity == 0 || source == NULL) return false;
	length = strlen(source);
	if (length >= capacity) return false;
	memcpy(destination, source, length + 1);
	return true;
}

static bool magnoom_is_path_separator(char c)
{
#if defined(_WIN32)
	return c == '/' || c == '\\';
#else
	return c == '/';
#endif
}

static const char *magnoom_get_file_extension(const char *filename)
{
	const char *basename, *dot, *cursor;

	if (filename == NULL || filename[0] == '\0') return NULL;
	basename = filename;
	for (cursor = filename; *cursor != '\0'; ++cursor) {
		if (magnoom_is_path_separator(*cursor)) basename = cursor + 1;
	}
	dot = strrchr(basename, '.');
	if (dot == NULL || dot == basename || dot[1] == '\0') return NULL;
	return dot;
}

static FileFormatEnum GetFileFormatFromExtension(const char *filename)
{
	const char *extension = magnoom_get_file_extension(filename);
	if (extension == NULL) return FILE_FORMAT_UNKNOWN;
	if (strcmp(extension, ".csv") == 0) return FILE_FORMAT_CSV;
	if (strcmp(extension, ".ovf") == 0) return FILE_FORMAT_OVF;
	if (strcmp(extension, ".vtk") == 0) return FILE_FORMAT_VTK;
	if (strcmp(extension, ".bin") == 0) return FILE_FORMAT_BIN;
	if (strcmp(extension, ".png") == 0) return FILE_FORMAT_PNG;
	return FILE_FORMAT_UNKNOWN;
}

static bool magnoom_file_format_can_import(FileFormatEnum format)
{
	return format == FILE_FORMAT_CSV || format == FILE_FORMAT_OVF ||
	       format == FILE_FORMAT_VTK || format == FILE_FORMAT_BIN;
}

static bool magnoom_file_format_can_export(FileFormatEnum format)
{
	return format > FILE_FORMAT_UNKNOWN && format < FILE_FORMAT_COUNT;
}

static bool magnoom_read_first_nonempty_line(FILE *file, char *line, size_t capacity)
{
	while (fgets(line, (int)capacity, file) != NULL) {
		const char *cursor = line;
		while (*cursor == ' ' || *cursor == '\t' || *cursor == '\r' || *cursor == '\n') ++cursor;
		if (*cursor != '\0') {
			if (cursor != line) memmove(line, cursor, strlen(cursor) + 1);
			return true;
		}
	}
	return false;
}

static bool magnoom_csv_line_is_spin(const char *line)
{
	const char *cursor = line;
	char *end;

	for (int component = 0; component < 6; ++component) {
		double value;
		while (*cursor == ' ' || *cursor == '\t') ++cursor;
		errno = 0;
		value = strtod(cursor, &end);
		if (end == cursor || errno == ERANGE || !isfinite(value) || fabs(value) > FLT_MAX) return false;
		cursor = end;
		while (*cursor == ' ' || *cursor == '\t') ++cursor;
		if (component < 5) {
			if (*cursor != ',') return false;
			++cursor;
		}
	}
	while (*cursor == ' ' || *cursor == '\t') ++cursor;
	if (*cursor == ',') ++cursor;
	while (*cursor == ' ' || *cursor == '\t' || *cursor == '\r' || *cursor == '\n') ++cursor;
	return *cursor == '\0';
}

static bool magnoom_validate_csv_content(FILE *file, const magnoom_ctx *ctx)
{
	char line[512];
	int rows = 0;

	if (ctx == NULL || ctx->NOS <= 0) return false;
	while (fgets(line, sizeof(line), file) != NULL) {
		size_t length = strlen(line);
		const char *cursor = line;
		if (length >= 120 || (length > 0 && line[length - 1] != '\n' && !feof(file))) return false;
		while (*cursor == ' ' || *cursor == '\t' || *cursor == '\r' || *cursor == '\n') ++cursor;
		if (*cursor == '\0') return false;
		if (!magnoom_csv_line_is_spin(cursor)) return false;
		rows++;
	}
	return rows == ctx->NOS;
}

static bool magnoom_validate_ovf_content(FILE *file, const magnoom_ctx *ctx)
{
	char line[512];
	int xnodes = 0, ynodes = 0, znodes = 0, valuedim = 0;
	bool data_marker = false;

	if (ctx == NULL || !magnoom_read_first_nonempty_line(file, line, sizeof(line)) ||
		strncmp(line, "# OOMMF OVF 2.0", strlen("# OOMMF OVF 2.0")) != 0) return false;
	while (fgets(line, sizeof(line), file) != NULL) {
		if (sscanf(line, "# valuedim: %d", &valuedim) == 1) continue;
		if (sscanf(line, "# xnodes: %d", &xnodes) == 1) continue;
		if (sscanf(line, "# ynodes: %d", &ynodes) == 1) continue;
		if (sscanf(line, "# znodes: %d", &znodes) == 1) continue;
		if (strstr(line, "# Begin: Data") != NULL) {
			data_marker = true;
			break;
		}
	}
	return data_marker && valuedim == 3 &&
	       xnodes == ctx->uABC[0] && ynodes == ctx->uABC[1] && znodes == ctx->uABC[2];
}

static bool magnoom_validate_vtk_content(FILE *file, const magnoom_ctx *ctx)
{
	char line[512], scalar_name[128], scalar_type[128];
	int xnodes = 0, ynodes = 0, znodes = 0, components = 0;
	bool lookup_table = false, payload_complete = false;

	if (ctx == NULL || !magnoom_read_first_nonempty_line(file, line, sizeof(line)) ||
		strncmp(line, "# vtk DataFile Version 2.0", strlen("# vtk DataFile Version 2.0")) != 0) return false;
	if (fgets(line, sizeof(line), file) == NULL) return false; /* title */
	if (fgets(line, sizeof(line), file) == NULL) return false;
	line[strcspn(line, "\r\n")] = '\0';
	if (strcmp(line, "BINARY") != 0) return false;
	if (fgets(line, sizeof(line), file) == NULL) return false;
	line[strcspn(line, "\r\n")] = '\0';
	if (strcmp(line, "DATASET STRUCTURED_POINTS") != 0) return false;
	while (fgets(line, sizeof(line), file) != NULL) {
		if (sscanf(line, "DIMENSIONS %d %d %d", &xnodes, &ynodes, &znodes) == 3) continue;
		if (sscanf(line, "SCALARS %127s %127s %d", scalar_name, scalar_type, &components) == 3) continue;
		if (strncmp(line, "LOOKUP_TABLE default", strlen("LOOKUP_TABLE default")) == 0) {
			long payload_start = ftell(file);
			lookup_table = true;
			if (payload_start >= 0 && xnodes > 0 && ynodes > 0 && znodes > 0) {
				size_t nx = (size_t)xnodes, ny = (size_t)ynodes, nz = (size_t)znodes;
				if (nx <= SIZE_MAX/ny && nx*ny <= SIZE_MAX/nz) {
					size_t cells = nx*ny*nz;
					if (cells <= (size_t)LONG_MAX/(3*sizeof(float))) {
						long payload_size = (long)(cells*3*sizeof(float));
						if (payload_start <= LONG_MAX - payload_size &&
							fseek(file, payload_start + payload_size, SEEK_SET) == 0) {
							int boundary = fgetc(file);
							payload_complete = boundary == EOF || boundary == '\r' || boundary == '\n';
						}
					}
				}
			}
			break;
		}
	}
	return lookup_table && payload_complete && components == 3 &&
	       strcmp(scalar_type, "float") == 0 &&
	       xnodes == ctx->uABC[0] && ynodes == ctx->uABC[1] && znodes == ctx->uABC[2];
}

static bool ValidateFileFormat(const magnoom_ctx *ctx, const char *path,
	FileFormatEnum format, char *error_message, size_t error_capacity)
{
	static const unsigned char png_signature[8] = {137, 80, 78, 71, 13, 10, 26, 10};
	FILE *file;
	bool valid = false;

	if (error_message == NULL || error_capacity == 0) return false;
	error_message[0] = '\0';
	file = fopen(path, "rb");
	if (file == NULL) {
		snprintf(error_message, error_capacity, "Cannot open the input file: %s", strerror(errno));
		error_message[error_capacity - 1] = '\0';
		return false;
	}

	switch (format) {
		case FILE_FORMAT_CSV:
			valid = magnoom_validate_csv_content(file, ctx);
			break;
		case FILE_FORMAT_OVF:
			valid = magnoom_validate_ovf_content(file, ctx);
			break;
		case FILE_FORMAT_VTK:
			valid = magnoom_validate_vtk_content(file, ctx);
			break;
		case FILE_FORMAT_BIN: {
			size_t cells = 0;
			long length = -1;
			if (ctx != NULL && ctx->uABC[0] > 0 && ctx->uABC[1] > 0 && ctx->uABC[2] > 0 &&
				ctx->AtomsPerBlock > 0) {
				size_t na = (size_t)ctx->uABC[0];
				size_t nb = (size_t)ctx->uABC[1];
				size_t nc = (size_t)ctx->uABC[2];
				size_t atoms = (size_t)ctx->AtomsPerBlock;
				if (na <= SIZE_MAX / nb && na*nb <= SIZE_MAX / nc &&
					na*nb*nc <= SIZE_MAX / atoms) {
					cells = na*nb*nc*atoms; /* one record per atom, not per unit cell */
					if (fseek(file, 0, SEEK_END) == 0) length = ftell(file);
				}
			}
			valid = length >= 0 && cells <= (size_t)LONG_MAX / sizeof(magnoom_bin_spin) &&
			        (size_t)length == cells*sizeof(magnoom_bin_spin);
			break;
		}
		case FILE_FORMAT_PNG: {
			unsigned char signature[sizeof(png_signature)];
			valid = fread(signature, 1, sizeof(signature), file) == sizeof(signature) &&
			        memcmp(signature, png_signature, sizeof(signature)) == 0;
			break;
		}
		default:
			break;
	}
	fclose(file);

	if (!valid) {
		if (format == FILE_FORMAT_CSV) {
			snprintf(error_message, error_capacity,
				"CSV content must contain exactly %d bounded six-value spin rows.",
				ctx != NULL ? ctx->NOS : 0);
		} else if (format == FILE_FORMAT_BIN) {
			snprintf(error_message, error_capacity,
				"BIN content size does not match the current simulation grid.");
		} else if (format == FILE_FORMAT_OVF || format == FILE_FORMAT_VTK) {
			snprintf(error_message, error_capacity,
				"The %s header, data marker, or dimensions do not match the current simulation grid.",
				fileFormatNames[format]);
		} else {
			snprintf(error_message, error_capacity,
				"The file content does not match the detected %s format.",
				format > FILE_FORMAT_UNKNOWN && format < FILE_FORMAT_COUNT ?
				fileFormatNames[format] : "unknown");
		}
		error_message[error_capacity - 1] = '\0';
	}
	return valid;
}

static bool magnoom_path_is_absolute(const char *path)
{
	if (path == NULL || path[0] == '\0') return false;
#if defined(_WIN32)
	if (magnoom_is_path_separator(path[0]) && magnoom_is_path_separator(path[1])) return true;
	return ((path[0] >= 'A' && path[0] <= 'Z') ||
	        (path[0] >= 'a' && path[0] <= 'z')) &&
	       path[1] == ':' && magnoom_is_path_separator(path[2]);
#else
	return path[0] == '/';
#endif
}

static bool magnoom_join_path(char *destination, size_t capacity,
	const char *directory, const char *filename)
{
	size_t directory_length, filename_length, required;
	bool add_separator;
	char separator;

	if (destination == NULL || capacity == 0 || directory == NULL ||
		filename == NULL || filename[0] == '\0') return false;
#if defined(_WIN32)
	if ((magnoom_is_path_separator(filename[0]) && !magnoom_is_path_separator(filename[1])) ||
		(((filename[0] >= 'A' && filename[0] <= 'Z') ||
		  (filename[0] >= 'a' && filename[0] <= 'z')) &&
		 filename[1] == ':' && !magnoom_is_path_separator(filename[2]))) return false;
#endif
	if (directory[0] == '\0' || magnoom_path_is_absolute(filename))
		return magnoom_copy_path(destination, capacity, filename);

	directory_length = strlen(directory);
	filename_length = strlen(filename);
	add_separator = !magnoom_is_path_separator(directory[directory_length - 1]);
	required = directory_length + (add_separator ? 1u : 0u) + filename_length + 1u;
	if (required > capacity) return false;

#if defined(_WIN32)
	separator = '\\';
#else
	separator = '/';
#endif
	memcpy(destination, directory, directory_length);
	if (add_separator) destination[directory_length++] = separator;
	memcpy(destination + directory_length, filename, filename_length + 1);
	return true;
}

static bool magnoom_replace_extension(char *destination, size_t capacity,
	const char *path, const char *extension)
{
	const char *basename, *dot, *cursor;
	size_t base_length, extension_length;

	if (destination == NULL || capacity == 0 || path == NULL || path[0] == '\0' ||
		extension == NULL || extension[0] == '\0') return false;
	basename = path;
	for (cursor = path; *cursor != '\0'; ++cursor) {
		if (magnoom_is_path_separator(*cursor)) basename = cursor + 1;
	}
	dot = strrchr(basename, '.');
	if (dot == basename) dot = NULL;
	base_length = dot != NULL ? (size_t)(dot - path) : strlen(path);
	extension_length = strlen(extension);
	if (base_length + extension_length + 1 > capacity) return false;
	memcpy(destination, path, base_length);
	memcpy(destination + base_length, extension, extension_length + 1);
	return true;
}

static bool magnoom_resolve_path(char *destination, size_t capacity,
	const char *directory, const char *filename)
{
	if (!magnoom_join_path(destination, capacity, directory, filename)) {
		fprintf(stderr, "Unable to resolve path from directory '%s' and file '%s': path is empty, unsupported, or too long.\n",
			directory != NULL ? directory : "(null)", filename != NULL ? filename : "(null)");
		return false;
	}
	return true;
}

static bool magnoom_resolve_path_with_extension(char *destination, size_t capacity,
	const char *directory, const char *filename, const char *extension)
{
	char path[MAGNOOM_PATH_CAPACITY];
	if (!magnoom_resolve_path(path, sizeof(path), directory, filename)) return false;
	if (!magnoom_replace_extension(destination, capacity, path, extension)) {
		fprintf(stderr, "Unable to replace the extension of '%s': path is too long.\n", path);
		return false;
	}
	return true;
}

static bool magnoom_executable_directory(char *destination, size_t capacity, const char *argv0)
{
	char executable[MAGNOOM_PATH_CAPACITY];
	bool found = false;
	char *last_separator = NULL;

#if defined(_WIN32)
	DWORD length = GetModuleFileNameA(NULL, executable, (DWORD)sizeof(executable));
	if (length > 0 && length < sizeof(executable)) {
		executable[length] = '\0';
		found = true;
	}
#elif defined(__APPLE__)
	uint32_t size = (uint32_t)sizeof(executable);
	if (_NSGetExecutablePath(executable, &size) == 0) found = true;
#elif defined(__linux__)
	ssize_t length = readlink("/proc/self/exe", executable, sizeof(executable) - 1);
	if (length > 0 && (size_t)length < sizeof(executable) - 1) {
		executable[length] = '\0';
		found = true;
	}
#endif

	if (!found && argv0 != NULL && argv0[0] != '\0')
		found = magnoom_copy_path(executable, sizeof(executable), argv0);
	if (!found) return magnoom_copy_path(destination, capacity, ".");

	for (char *cursor = executable; *cursor != '\0'; ++cursor) {
		if (magnoom_is_path_separator(*cursor)) last_separator = cursor;
	}
	if (last_separator == NULL) return magnoom_copy_path(destination, capacity, ".");
	if (last_separator == executable ||
		(last_separator == executable + 2 && executable[1] == ':')) {
		last_separator[1] = '\0';
	} else {
		*last_separator = '\0';
	}
	return magnoom_copy_path(destination, capacity, executable);
}

static bool magnoom_resolve_input_path(const magnoom_ctx *ctx,
	char *destination, size_t capacity)
{
	return magnoom_resolve_path(destination, capacity,
		ctx->input_directory, ctx->inputfilename);
}

static bool magnoom_resolve_output_path(magnoom_ctx *ctx,
	const char *filename, char *destination, size_t capacity)
{
	char directory[MAGNOOM_PATH_CAPACITY];
	pthread_mutex_lock(&ctx->record_mutex);
	magnoom_copy_path(directory, sizeof(directory), ctx->output_directory);
	pthread_mutex_unlock(&ctx->record_mutex);
	return magnoom_resolve_path(destination, capacity, directory, filename);
}

static bool magnoom_resolve_output_path_with_extension(magnoom_ctx *ctx,
	const char *extension, char *destination, size_t capacity)
{
	char directory[MAGNOOM_PATH_CAPACITY];
	pthread_mutex_lock(&ctx->record_mutex);
	magnoom_copy_path(directory, sizeof(directory), ctx->output_directory);
	pthread_mutex_unlock(&ctx->record_mutex);
	return magnoom_resolve_path_with_extension(destination, capacity,
		directory, ctx->outputfilename, extension);
}

static bool magnoom_flush_records_locked(magnoom_ctx *ctx)
{
	bool success = ctx->outFile != NULL;
	if (ctx->recordsCounter <= 0) return true;
	if (ctx->outFile != NULL) {
		for (int i = 0; i < ctx->recordsCounter; ++i) {
			snprintf(ctx->BufferString, sizeof(ctx->BufferString),
				"%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,\n",
				ctx->RecordingBuffer[0][i], ctx->RecordingBuffer[0][i]*ctx->t_step,
				ctx->RecordingBuffer[1][i], ctx->RecordingBuffer[2][i],
				ctx->RecordingBuffer[3][i], ctx->RecordingBuffer[4][i]);
			if (fputs(ctx->BufferString, ctx->outFile) == EOF) {
				success = false;
				break;
			}
		}
		if (success && fflush(ctx->outFile) == EOF) success = false;
	}
	if (success) {
		printf("Recording to file %s is done!\n", ctx->record_path);
	} else {
		fprintf(stderr, "Unable to write buffered records to '%s': %s\n",
			ctx->record_path[0] != '\0' ? ctx->record_path : "table.csv",
			ctx->outFile != NULL ? strerror(errno) : "file is not open");
	}
	ctx->recordsCounter = 0;
	return success;
}

static FILE *magnoom_open_record_stream(const char *directory,
	char *resolved_path, size_t capacity)
{
	FILE *stream;
	if (!magnoom_resolve_path(resolved_path, capacity, directory, "table.csv")) return NULL;
	stream = fopen(resolved_path, "w");
	if (stream == NULL) {
		fprintf(stderr, "Cannot open output file '%s': %s\n", resolved_path, strerror(errno));
		return NULL;
	}
	if (fputs("iter,time,Mx,My,Mz,E_tot,\n", stream) == EOF) {
		fprintf(stderr, "Cannot initialize output file '%s': %s\n", resolved_path, strerror(errno));
		fclose(stream);
		return NULL;
	}
	return stream;
}

static bool magnoom_change_output_directory(magnoom_ctx *ctx, const char *directory)
{
	char checked_directory[MAGNOOM_PATH_CAPACITY];
	char new_path[MAGNOOM_PATH_CAPACITY];
	char old_path[MAGNOOM_PATH_CAPACITY];
	FILE *new_stream;

	if (!magnoom_copy_path(checked_directory, sizeof(checked_directory), directory)) {
		fprintf(stderr, "Output directory is too long.\n");
		return false;
	}
	if (!magnoom_resolve_path(new_path, sizeof(new_path), checked_directory, "table.csv")) return false;
	if (strcmp(new_path, ctx->record_path) == 0 && ctx->outFile != NULL) {
		pthread_mutex_lock(&ctx->record_mutex);
		bool copied = magnoom_copy_path(ctx->output_directory,
			sizeof(ctx->output_directory), checked_directory);
		pthread_mutex_unlock(&ctx->record_mutex);
		return copied;
	}
	pthread_mutex_lock(&ctx->record_mutex);
	magnoom_flush_records_locked(ctx);
	magnoom_copy_path(old_path, sizeof(old_path), ctx->record_path);
	if (ctx->outFile != NULL) fclose(ctx->outFile);
	ctx->outFile = NULL;
	new_stream = magnoom_open_record_stream(checked_directory, new_path, sizeof(new_path));
	if (new_stream == NULL) {
		if (old_path[0] != '\0') {
			ctx->outFile = fopen(old_path, "a");
			if (ctx->outFile == NULL) {
				fprintf(stderr, "Cannot reopen output file '%s': %s\n", old_path, strerror(errno));
			}
		}
		pthread_mutex_unlock(&ctx->record_mutex);
		return false;
	}
	ctx->outFile = new_stream;
	magnoom_copy_path(ctx->output_directory, sizeof(ctx->output_directory), checked_directory);
	magnoom_copy_path(ctx->record_path, sizeof(ctx->record_path), new_path);
	pthread_mutex_unlock(&ctx->record_mutex);
	return true;
}

static bool magnoom_reset_record_file(magnoom_ctx *ctx)
{
	char path[MAGNOOM_PATH_CAPACITY];
	FILE *stream;

	pthread_mutex_lock(&ctx->record_mutex);
	if (ctx->outFile != NULL) fclose(ctx->outFile);
	ctx->outFile = NULL;
	ctx->record_path[0] = '\0';
	ctx->recordsCounter = 0;
	stream = magnoom_open_record_stream(ctx->output_directory, path, sizeof(path));
	if (stream != NULL) {
		ctx->outFile = stream;
		magnoom_copy_path(ctx->record_path, sizeof(ctx->record_path), path);
	}
	pthread_mutex_unlock(&ctx->record_mutex);
	return stream != NULL;
}

#ifndef MAGNOOM_NO_MAIN
static void magnoom_close_record_file(magnoom_ctx *ctx)
{
	pthread_mutex_lock(&ctx->record_mutex);
	magnoom_flush_records_locked(ctx);
	if (ctx->outFile != NULL) fclose(ctx->outFile);
	ctx->outFile = NULL;
	pthread_mutex_unlock(&ctx->record_mutex);
}
#endif

/*****************************************************************************/
/* Interleaved vector-field accessors: Sx/Sy/Sz-style per-component arrays   */
/* are being migrated to single interleaved buffers (S[3*n+0..2], one        */
/* allocation per vector field instead of three). These macros are the only */
/* sanctioned way to index such a buffer -- do not write buf[3*n+d] inline.  */
/* IMAGE_COMPONENT addresses one flat snapshot block laid out as            */
/* [imageIndex][n][d]; a specific snapshot's sub-buffer (a plain double*)   */
/* can still be read with VEC_X/VEC_Y/VEC_Z once offset to its start.       */
/*****************************************************************************/
#define VEC_COMPONENT(buf, n, d) ((buf)[3*(size_t)(n) + (size_t)(d)])
#define VEC_X(buf, n) VEC_COMPONENT(buf, n, 0)
#define VEC_Y(buf, n) VEC_COMPONENT(buf, n, 1)
#define VEC_Z(buf, n) VEC_COMPONENT(buf, n, 2)
#define IMAGE_COMPONENT(images, imageIndex, nodesNum, n, d) \
	((images)[(((size_t)(imageIndex)*(size_t)(nodesNum) + (size_t)(n))*3 + (size_t)(d))])

static void magnoom_stop_engine(magnoom_ctx *ctx)
{
	bool idle;

	ctx->Play = 0;
	pthread_mutex_lock(&ctx->culc_mutex);
	if (ctx->EngineRunState != WAIT) ctx->EngineRunState = STOP_REQUESTED;
	ctx->SleepTime = 3000;
	idle = ctx->EngineIdle;
	pthread_mutex_unlock(&ctx->culc_mutex);

	while (!idle) {
		glfwSleep(0.001);
		pthread_mutex_lock(&ctx->culc_mutex);
		idle = ctx->EngineIdle;
		pthread_mutex_unlock(&ctx->culc_mutex);
	}
}

static void magnoom_reset_solver_state(magnoom_ctx *ctx)
{
	for (int i = 0; i < ctx->NOS; ++i) {
		VEC_X(ctx->bS, i) = VEC_X(ctx->tS, i) = VEC_X(ctx->S, i);
		VEC_Y(ctx->bS, i) = VEC_Y(ctx->tS, i) = VEC_Y(ctx->S, i);
		VEC_Z(ctx->bS, i) = VEC_Z(ctx->tS, i) = VEC_Z(ctx->S, i);
	}

	if (ctx->NoiseX != NULL) memset(ctx->NoiseX, 0, (size_t)ctx->NOS*sizeof(*ctx->NoiseX));
	if (ctx->NoiseY != NULL) memset(ctx->NoiseY, 0, (size_t)ctx->NOS*sizeof(*ctx->NoiseY));
	if (ctx->NoiseZ != NULL) memset(ctx->NoiseZ, 0, (size_t)ctx->NOS*sizeof(*ctx->NoiseZ));
}

/*****************************************************************************/
/* Replace the Cartesian atom positions in the basic unit cell. This is a    */
/* startup-only operation: call it after abc/uABC are final and before       */
/* neighbor maps, spin arrays, drawing arrays, or worker threads exist.      */
/*****************************************************************************/
bool magnoom_ctx_set_block(magnoom_ctx *ctx, int atom_count, const float positions[][3])
{
	if (ctx == NULL || positions == NULL || atom_count <= 0) {
		fprintf(stderr, "magnoom_ctx_set_block: expected a context, atom positions, and a positive atom count.\n");
		return false;
	}
	if (atom_count > MAX_ATOMS_PER_BLOCK) {
		fprintf(stderr, "magnoom_ctx_set_block: at most %d atoms per block are supported.\n", MAX_ATOMS_PER_BLOCK);
		return false;
	}
	if (ctx->Neighbors.PairCount != 0 || ctx->Neighbors.AIdxBlock != NULL || ctx->S != NULL ||
		ctx->bS != NULL || ctx->PosX != NULL || ctx->vertices != NULL) {
		fprintf(stderr, "magnoom_ctx_set_block: the block can only be changed before derived geometry is allocated.\n");
		return false;
	}
	if (ctx->uABC[0] <= 0 || ctx->uABC[1] <= 0 || ctx->uABC[2] <= 0) {
		fprintf(stderr, "magnoom_ctx_set_block: lattice dimensions must be positive.\n");
		return false;
	}

	double a[3], b[3], c[3];
	for (int d = 0; d < 3; ++d) {
		a[d] = ctx->abc[0][d];
		b[d] = ctx->abc[1][d];
		c[d] = ctx->abc[2][d];
		if (!isfinite(a[d]) || !isfinite(b[d]) || !isfinite(c[d])) {
			fprintf(stderr, "magnoom_ctx_set_block: lattice vectors must be finite.\n");
			return false;
		}
	}

	double b_cross_c[3] = {
		b[1]*c[2] - b[2]*c[1],
		b[2]*c[0] - b[0]*c[2],
		b[0]*c[1] - b[1]*c[0]
	};
	double c_cross_a[3] = {
		c[1]*a[2] - c[2]*a[1],
		c[2]*a[0] - c[0]*a[2],
		c[0]*a[1] - c[1]*a[0]
	};
	double a_cross_b[3] = {
		a[1]*b[2] - a[2]*b[1],
		a[2]*b[0] - a[0]*b[2],
		a[0]*b[1] - a[1]*b[0]
	};
	double determinant = a[0]*b_cross_c[0] + a[1]*b_cross_c[1] + a[2]*b_cross_c[2];
	double a_length = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
	double b_length = sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
	double c_length = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
	double volume_scale = a_length*b_length*c_length;
	if (!isfinite(volume_scale) || volume_scale <= 0.0 ||
		fabs(determinant) <= volume_scale*(32.0*FLT_EPSILON)) {
		fprintf(stderr, "magnoom_ctx_set_block: lattice vectors do not define a valid three-dimensional unit cell.\n");
		return false;
	}

	for (int atom = 0; atom < atom_count; ++atom) {
		double p[3] = {positions[atom][0], positions[atom][1], positions[atom][2]};
		if (!isfinite(p[0]) || !isfinite(p[1]) || !isfinite(p[2])) {
			fprintf(stderr, "magnoom_ctx_set_block: atom %d has a nonfinite position.\n", atom);
			return false;
		}
		double fractional[3] = {
			(p[0]*b_cross_c[0] + p[1]*b_cross_c[1] + p[2]*b_cross_c[2])/determinant,
			(p[0]*c_cross_a[0] + p[1]*c_cross_a[1] + p[2]*c_cross_a[2])/determinant,
			(p[0]*a_cross_b[0] + p[1]*a_cross_b[1] + p[2]*a_cross_b[2])/determinant
		};
		for (int d = 0; d < 3; ++d) {
			if (!isfinite(fractional[d]) || fractional[d] < 0.0 || fractional[d] >= 1.0) {
				fprintf(stderr,
					"magnoom_ctx_set_block: atom %d at (%g, %g, %g) is outside the half-open unit cell; fractional position is (%g, %g, %g).\n",
					atom, p[0], p[1], p[2], fractional[0], fractional[1], fractional[2]);
				return false;
			}
		}
	}

	size_t na = (size_t)ctx->uABC[0];
	size_t nb = (size_t)ctx->uABC[1];
	size_t nc = (size_t)ctx->uABC[2];
	if (na > SIZE_MAX/nb || na*nb > SIZE_MAX/nc) {
		fprintf(stderr, "magnoom_ctx_set_block: lattice dimensions overflow the supported size.\n");
		return false;
	}
	size_t block_count = na*nb*nc;
	/* ARROW1 at the UI limit of 20 faces uses 540 scalar vertex components per spin. */
	if ((size_t)atom_count > SIZE_MAX/block_count ||
		(size_t)atom_count*block_count > INT_MAX/540 ||
		(size_t)atom_count > SIZE_MAX/sizeof(*ctx->Block)) {
		fprintf(stderr, "magnoom_ctx_set_block: the requested geometry exceeds the supported index range.\n");
		return false;
	}

	float (*new_block)[3] = (float (*)[3])malloc((size_t)atom_count*sizeof(*new_block));
	if (new_block == NULL) {
		fprintf(stderr, "magnoom_ctx_set_block: unable to allocate the new block.\n");
		return false;
	}
	memcpy(new_block, positions, (size_t)atom_count*sizeof(*new_block));

	float (*old_block)[3] = ctx->Block;
	ctx->Block = new_block;
	ctx->AtomsPerBlock = atom_count;
	ctx->NOB = (int)block_count;
	ctx->NOB_AL = ctx->uABC[1]*ctx->uABC[2];
	ctx->NOB_BL = ctx->uABC[0]*ctx->uABC[2];
	ctx->NOB_CL = ctx->uABC[0]*ctx->uABC[1];
	ctx->NOS = atom_count*ctx->NOB;
	ctx->NOS_AL = atom_count*ctx->NOB_AL;
	ctx->NOS_BL = atom_count*ctx->NOB_BL;
	ctx->NOS_CL = atom_count*ctx->NOB_CL;
	ctx->iNOS = 1.0/ctx->NOS;
	ctx->NOSK = 0;
	ctx->GreedFilterMaxA = ctx->uABC[0]-1;
	ctx->GreedFilterMaxB = ctx->uABC[1]-1;
	ctx->GreedFilterMaxC = ctx->uABC[2]-1;
	free(old_block);
	return true;
}

/*****************************************************************************/
/* magnoom_ctx_init: sets every field to its current compile-time default,   */
/* folding in the former InitializeGlobalState()'s derived-value computation.*/
/* Returns false (matching InitializeGlobalState's contract) if the geometry */
/* allocations fail.                                                        */
/*****************************************************************************/
bool magnoom_ctx_init(magnoom_ctx *ctx)
{
	if (ctx == NULL) return false;
	memset(ctx, 0, sizeof(*ctx));

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
	// ctx->Jij[0]=1.0f; ctx->Jij[1]=0.0f; ctx->Jij[2]=0.0f; ctx->Jij[3]=0.0f;
	// ctx->Jij[4]=0.0f; ctx->Jij[5]=0.0f;
	/*EuSi v0*/
	ctx->Jij[0]= 1.0f; ctx->Jij[1]= 0.8f; ctx->Jij[2]= 0.5f; 
	ctx->Jij[3]= 0.5f; ctx->Jij[4]=-0.3f; ctx->Jij[5]=-0.174574312f;

	ctx->Bij[0]=0.0f; ctx->Bij[1]=0.0f; ctx->Bij[2]=0.0f;
	ctx->Bij[3]=0.0f; ctx->Bij[4]=0.0f; ctx->Bij[5]=0.0f;
	ctx->Dij[0]=0.0f; ctx->Dij[1]=0.0f; ctx->Dij[2]=0.0f;
	ctx->Dij[3]=0.0f; ctx->Dij[4]=0.0f; ctx->Dij[5]=0.0f;

	/* Magnetocrystalline anisotropy */
	ctx->anisotropy_mode = ANISOTROPY_GLOBAL;
	/* Tetragonal anisotropy*/
	/*EuSi v0*/
	/* E = K1 sin^2T + K2 sin^4T + K3 sin^4T cos4F*/
	double K1 = 0.0;//-0.1;
	double K2 = 0.0;//0.1;
	double K3 = 0.0;

	/* K11=K22=-K1, K33=0; K1111=K2222=-(K2+4K3), K1122=-K2/3, K3333=-3K3,
	 * K1133=K2233=-K3 -- see the wiki's Anisotropy (F6) page. */
	if (!k2_set(ctx->anisotropy_local[0].K2, 0, 0, -K1) ||
		!k2_set(ctx->anisotropy_local[0].K2, 1, 1, -K1) ||
		!k4_set(ctx->anisotropy_local[0].K4, 0, 0, 0, 0, -(K2+4*K3)) ||
		!k4_set(ctx->anisotropy_local[0].K4, 1, 1, 1, 1, -(K2+4*K3)) ||
		!k4_set(ctx->anisotropy_local[0].K4, 0, 0, 1, 1, -K2/3) ||
		!k4_set(ctx->anisotropy_local[0].K4, 2, 2, 2, 2, -3*K3) ||
		!k4_set(ctx->anisotropy_local[0].K4, 0, 0, 2, 2, -K3) ||
		!k4_set(ctx->anisotropy_local[0].K4, 1, 1, 2, 2, -K3)) {
		fprintf(stderr, "magnoom_ctx_init: invalid compiled-in default anisotropy tensor component.\n");
		return false;
	}

	for (int atom = 0; atom < MAX_ATOMS_PER_BLOCK; ++atom) {
		ctx->anisotropy_quaternion[atom][3] = 1.0;
	}
	ctx->anisotropy_rotation_axis[0] = 0.0f;
	ctx->anisotropy_rotation_axis[1] = 0.0f;
	ctx->anisotropy_rotation_axis[2] = 1.0f;
	ctx->anisotropy_rotation_angle = 0.0f;
	ctx->anisotropy_map_theta_step = PI / 40.0;
	ctx->anisotropy_map_phi_step = PI / 80.0;
	ctx->anisotropy_map_vtk_format = LEGACY_ASCII_VTK;

	/* External magnetic field: static (DC) and time-dependent (AC) components */
	ctx->BextDCDirection[0]=0.0f; ctx->BextDCDirection[1]=0.0f; ctx->BextDCDirection[2]=1.0f;
	ctx->BextACDirection[0]=0.0f; ctx->BextACDirection[1]=0.0f; ctx->BextACDirection[2]=1.0f;
	ctx->BextACOmega = 0.005;
	ctx->BextACPulseWidth = 20.0f;
	ctx->BextACWaveform = BEXT_AC_SIN;

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
	ctx->Max_Numb_Iteration = 100000;

	/* FPS & IPS / recording */
	ctx->rec_iteration = 100;
	ctx->rec_num_mode = 10;
	ctx->num_images = 32;
	ctx->SleepTime = 1000;
	ctx->input_directory[0] = '\0';
	ctx->output_directory[0] = '\0';
	strcpy(ctx->inputfilename, "input.csv");
	strcpy(ctx->outputfilename, "output.csv");

	/* Geometry */
	ctx->uABC[0]=10; ctx->uABC[1]=10; ctx->uABC[2]=10;
	ctx->ShellNumber = 1;

	/* Keep exactly one complete crystal basis below active. */
	/* Simple cubic, one atom: */
	// const float basis[][3] = {{0.5f, 0.5f, 0.5f}};
	// ctx->abc[0][0]=1.0f; ctx->abc[0][1]=0.0f; ctx->abc[0][2]=0.0f;
	// ctx->abc[1][0]=0.0f; ctx->abc[1][1]=1.0f; ctx->abc[1][2]=0.0f;
	// ctx->abc[2][0]=0.0f; ctx->abc[2][1]=0.0f; ctx->abc[2][2]=1.0f;
	// int atom_count = (int)(sizeof(basis)/sizeof(basis[0]));

	// B20 basis (u = 0.138), cubic unit cell:
	// const float uB20 = 0.138f;
	// const float basis[][3] = {
	//     {0.0f,             0.0f,             0.0f},
	//     {0.5f,             0.5f-2.0f*uB20, 1.0f-2.0f*uB20},
	//     {1.0f-2.0f*uB20, 0.5f,             0.5f-2.0f*uB20},
	//     {0.5f-2.0f*uB20, 1.0f-2.0f*uB20, 0.5f}
	// };
	// ctx->abc[0][0]=1.0f; ctx->abc[0][1]=0.0f; ctx->abc[0][2]=0.0f;
	// ctx->abc[1][0]=0.0f; ctx->abc[1][1]=1.0f; ctx->abc[1][2]=0.0f;
	// ctx->abc[2][0]=0.0f; ctx->abc[2][1]=0.0f; ctx->abc[2][2]=1.0f;
	// int atom_count = (int)(sizeof(basis)/sizeof(basis[0]));

	/* EuSi fractional coordinates converted to normalized Cartesian positions:*/
	const float c_EuSi = 3.9845f;
	const float a_EuSi = 4.6955f/c_EuSi;
	const float b_EuSi = 11.1528f/c_EuSi;
	const float u_Eu = 0.3595f;
	const float basis[][3] = {
	    {0.25f, 0.0f,          u_Eu*b_EuSi},
	    {0.75f, 0.0f,          (1.0f-u_Eu)*b_EuSi},
	    {0.75f, 0.5f*a_EuSi, (0.5f-u_Eu)*b_EuSi},
	    {0.25f, 0.5f*a_EuSi, (0.5f+u_Eu)*b_EuSi}
	};
	ctx->abc[0][0]=1.0f; ctx->abc[0][1]=0.0f;   ctx->abc[0][2]=0.0f;
	ctx->abc[1][0]=0.0f; ctx->abc[1][1]=a_EuSi; ctx->abc[1][2]=0.0f;
	ctx->abc[2][0]=0.0f; ctx->abc[2][1]=0.0f;   ctx->abc[2][2]=b_EuSi;
	int atom_count = (int)(sizeof(basis)/sizeof(basis[0]));

	// FCC2 basis, orthogonal unit cell:
	// const float sqrt2 = sqrtf(2.0f);
	// const float basis[][3] = {
	//     {0.0f, 0.0f,             0.0f},
	//     {0.0f, sqrt2/2.0f,      0.0f},
	//     {0.5f, sqrt2/4.0f,      sqrt2/4.0f},
	//     {0.5f, 3.0f*sqrt2/4.0f, sqrt2/4.0f},
	//     {0.0f, sqrt2/2.0f,      sqrt2/2.0f}
	// };
	// ctx->abc[0][0]=1.0f; ctx->abc[0][1]=0.0f;  ctx->abc[0][2]=0.0f;
	// ctx->abc[1][0]=0.0f; ctx->abc[1][1]=sqrt2; ctx->abc[1][2]=0.0f;
	// ctx->abc[2][0]=0.0f; ctx->abc[2][1]=0.0f;  ctx->abc[2][2]=sqrt2;
	// int atom_count = (int)(sizeof(basis)/sizeof(basis[0]));

	// FCC3 basis, nonorthogonal unit cell:
	// const float sqrt2 = sqrtf(2.0f);
	// const float sqrt3 = sqrtf(3.0f);
	// const float sqrt6 = sqrtf(6.0f);
	// const float basis[][3] = {{0.0f, 0.0f, 0.0f}};
	// ctx->abc[0][0]=sqrt2;      ctx->abc[0][1]=0.0f;       ctx->abc[0][2]=0.0f;
	// ctx->abc[1][0]=sqrt2/2.0f; ctx->abc[1][1]=sqrt6/2.0f; ctx->abc[1][2]=0.0f;
	// ctx->abc[2][0]=sqrt2/2.0f; ctx->abc[2][1]=sqrt6/6.0f; ctx->abc[2][2]=sqrt3/sqrt2;
	// int atom_count = (int)(sizeof(basis)/sizeof(basis[0]));

	/* magnoom.c concurrency primitives */
	ctx->EngineRunState = WAIT;
	ctx->DataTransferState = WAIT_DATA;
	ctx->EngineIdle = true;
	ctx->EngineShutdown = false;
	ctx->EngineShutdownRequested = false;

	/* visualization.c: window / GLFW / display state */
	static char window_title[64];
	snprintf(window_title, sizeof(window_title), "%s v%s", SOFTWARE_NAME, SOFTWARE_VERSION);
	ctx->WINDOWTITLE = window_title;
	ctx->window_width = 1400;
	ctx->window_height = 800;
	ctx->ContentScaleX = 1.0;
	ctx->ContentScaleY = 1.0;

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
	ctx->Scale_BextDC = 10.0f;
	ctx->axis[0]=0.7f; ctx->axis[1]=0.7f; ctx->axis[2]=0.0f;
	ctx->angle = 0.8f;

	/* OpenGL material / lighting */
	ctx->ambient[0]=0.33f; ctx->ambient[1]=0.22f; ctx->ambient[2]=0.03f; ctx->ambient[3]=0.0f;
	ctx->diffuse[0]=0.78f; ctx->diffuse[1]=0.57f; ctx->diffuse[2]=0.11f; ctx->diffuse[3]=1.0f;
	ctx->specular[0]=0.1f; ctx->specular[1]=0.1f; ctx->specular[2]=0.08f; ctx->specular[3]=1.0f;
	ctx->shininess = 128.0f;
	ctx->LightMultiplier = 1.0f;
	ctx->LightDirection[0]=0.0f; ctx->LightDirection[1]=0.0f; ctx->LightDirection[2]=1.0f;
	ctx->WhichLightingMode = LIGHT_ADAPTIVE;

	/* initial-state generation parameters */
	ctx->WhichInitialState = RND;
	ctx->WhichGeometry = DEFAULT_G;
	ctx->chSizeG = 50.0f;
	ctx->chSize = 12.0f;
	ctx->chDir[0]=0.0f; ctx->chDir[1]=1.0f; ctx->chDir[2]=0.0f;

	/* color scheme (visualization) */
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
	ctx->arrowFaces_BextDC = 30;
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
	ctx->BextDC[0] = ctx->BextDCMagnitude*ctx->BextDCDirection[0];
	ctx->BextDC[1] = ctx->BextDCMagnitude*ctx->BextDCDirection[1];
	ctx->BextDC[2] = ctx->BextDCMagnitude*ctx->BextDCDirection[2];
	ctx->BextACAmplitude = ctx->BextDCMagnitude*(float)sin(PI/180.0);
	ctx->BextACPeriod = TPI/ctx->BextACOmega;

	if (!magnoom_ctx_set_block(ctx, atom_count, basis)) {
		return false;
	}

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




void Save_OVF_b8(magnoom_ctx *ctx, double* S, const char *ovf_filename){
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
         snprintf(ctx->ShortBuffer,80,"# xmax: %.6g\n",ctx->uABC[0]*temp1*a_lattice);
         fputs (ctx->ShortBuffer,pFile);

         temp0  = ctx->abc[1][0]*ctx->abc[1][0];
         temp0 += ctx->abc[1][1]*ctx->abc[1][1];
         temp0 += ctx->abc[1][2]*ctx->abc[1][2];
         temp2 = sqrt(temp0);
         snprintf(ctx->ShortBuffer,80,"# ymax: %.6g\n",ctx->uABC[1]*temp2*a_lattice);
         fputs (ctx->ShortBuffer,pFile);

         temp0  = ctx->abc[2][0]*ctx->abc[2][0];
         temp0 += ctx->abc[2][1]*ctx->abc[2][1];
         temp0 += ctx->abc[2][2]*ctx->abc[2][2];
         temp3 = sqrt(temp0);
         snprintf(ctx->ShortBuffer,80,"# zmax: %.6g\n",ctx->uABC[2]*temp3*a_lattice);
         fputs (ctx->ShortBuffer,pFile);
         fputs ("# valuedim: 3\n",pFile);
         fputs ("# valuelabels: m_x m_y m_z\n",pFile);
         fputs ("# valueunits: 1 1 1\n",pFile);
         fputs ("# Desc: Total simulation time:  0  s\n",pFile);

         snprintf(ctx->ShortBuffer,80,"# xbase: %.6g\n",temp1*0.5*a_lattice);
         fputs (ctx->ShortBuffer,pFile);

         snprintf(ctx->ShortBuffer,80,"# ybase: %.6g\n",temp2*0.5*a_lattice);
         fputs (ctx->ShortBuffer,pFile);

         snprintf(ctx->ShortBuffer,80,"# zbase: %.6g\n",temp3*0.5*a_lattice);
         fputs (ctx->ShortBuffer,pFile);

         snprintf(ctx->ShortBuffer,80,"# xnodes: %d\n",ctx->uABC[0]);
         fputs (ctx->ShortBuffer,pFile);
         snprintf(ctx->ShortBuffer,80,"# ynodes: %d\n",ctx->uABC[1]);
         fputs (ctx->ShortBuffer,pFile);
         snprintf(ctx->ShortBuffer,80,"# znodes: %d\n",ctx->uABC[2]);
         fputs (ctx->ShortBuffer,pFile);

         snprintf(ctx->ShortBuffer,80,"# xstepsize:  %.6g\n",temp1*a_lattice);
         fputs (ctx->ShortBuffer,pFile);

         snprintf(ctx->ShortBuffer,80,"# ystepsize: %.6g\n",temp2*a_lattice);
         fputs (ctx->ShortBuffer,pFile);

         snprintf(ctx->ShortBuffer,80,"# zstepsize: %.6g\n",temp3*a_lattice);
         fputs (ctx->ShortBuffer,pFile);

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
                         double Temp[]= {VEC_X(S,N), VEC_Y(S,N), VEC_Z(S,N)};
                         fwrite (Temp , sizeof(double), 3, pFile);
                     }
                 }
             }
         }
        fputs ("# End: Data Binary 4\n",pFile);
        fputs ("# End: Segment\n",pFile);
        fclose (pFile);
		printf("Recording to the file %s is done!\n", ovf_filename);
	} else {
		fprintf(stderr, "Cannot open output file '%s': %s\n", ovf_filename, strerror(errno));
    }
}


void SaveBin(magnoom_ctx *ctx, double* S, const char *bin_filename){
    unsigned short int num = 65535;
    FILE * pFile = fopen (bin_filename,"wb");
	if (pFile == NULL) {
		fprintf(stderr, "Cannot open output file '%s': %s\n", bin_filename, strerror(errno));
		return;
	}
    for(int k = 0; k < ctx->uABC[2]; k++){
        for(int j = 0; j < ctx->uABC[1]; j++){
            for(int i = 0; i < ctx->uABC[0]; i++){
            int n = i+j*ctx->uABC[0]+k*ctx->uABC[0]*ctx->uABC[1]; /* index of the block */
            n = n*ctx->AtomsPerBlock; /* index of the first spin in the block */

            for (int atom = 0; atom < ctx->AtomsPerBlock; atom++){
                int N = n + atom;
                double nx = VEC_X(S,N), ny = VEC_Y(S,N), nz = VEC_Z(S,N);

                double T, F;
                T = acos(nz)/PI;
                F = atan2(ny,nx)/PI;
                if(F <= 0) F += 2.0;
                F /= 2;

                unsigned short int p=0, q=0;

                q = T*num;
                p = F*num;

                magnoom_bin_spin my_par = {q, p};
                fwrite(&my_par, sizeof(my_par), 1, pFile);
            }

      }
    }
  }
  fclose (pFile);
  printf("Recording to the file %s is done!\n", bin_filename);
}


void SavePng(magnoom_ctx *ctx, double* S, const char *png_filename, enSliceMode WhichSliceMode, int x1, int y1, int z1){
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
        vec[0] = VEC_X(S,n); vec[1] = VEC_Y(S,n); vec[2] = VEC_Z(S,n);
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

void Save_VTS_b4(magnoom_ctx *ctx, double* Sx, double* Sy, double* Sz, float * Px, float * Py, float * Pz, float box[][3], const char *vts_filename){
    float a_lattice = 1.0e-9;
    FILE * pFile = fopen (vts_filename,"wb");
    if(pFile!=NULL) {
        fputs ("<?xml version=\"1.0\"?>\n",pFile);
        fputs ("<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n",pFile);
        snprintf(ctx->ShortBuffer,80,"\t<StructuredGrid WholeExtent=\"0 %d 0 %d 0 %d\">\n",ctx->uABC[0]-1, ctx->uABC[1]-1, ctx->uABC[2]-1);
        fputs (ctx->ShortBuffer,pFile);
        snprintf(ctx->ShortBuffer,80,"\t\t<Piece Extent=\"0 %d 0 %d 0 %d \">\n",ctx->uABC[0]-1, ctx->uABC[1]-1, ctx->uABC[2]-1);
        fputs (ctx->ShortBuffer,pFile);
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

void Save_VTS_ascii(magnoom_ctx *ctx, double* Sx, double* Sy, double* Sz, float * Px, float * Py, float * Pz, float box[][3], const char *vts_filename){
    float a_lattice = 1.0;//.0e-9;
    FILE * pFile = fopen (vts_filename,"wb");
    if(pFile!=NULL) {
        fputs ("<?xml version=\"1.0\"?>\n",pFile);
        fputs ("<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n",pFile);
        snprintf(ctx->ShortBuffer,80,"\t<StructuredGrid WholeExtent=\"0 %d 0 %d 0 %d\">\n",ctx->uABC[0]-1, ctx->uABC[1]-1, ctx->uABC[2]-1);
        fputs (ctx->ShortBuffer,pFile);
        snprintf(ctx->ShortBuffer,80,"\t\t<Piece Extent=\"0 %d 0 %d 0 %d \">\n",ctx->uABC[0]-1, ctx->uABC[1]-1, ctx->uABC[2]-1);
        fputs (ctx->ShortBuffer,pFile);
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
                        snprintf(ctx->ShortBuffer,80,"%.6g %.6g %.6g ",(Px[N]+Tr[0])*a_lattice, (Py[N]+Tr[1])*a_lattice, (Pz[N]+Tr[2])*a_lattice);
                        fputs (ctx->ShortBuffer,pFile);
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
                        snprintf(ctx->ShortBuffer,80,"%.6g %.6g %.6g ",Sx[N], Sy[N], Sz[N]);
                        fputs (ctx->ShortBuffer,pFile);
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

void Save_VTK(magnoom_ctx *ctx, double* S, const int mode, const char *vtk_filename)
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
        snprintf(ctx->ShortBuffer,80,"%s\n",(mode==0 ? "BINARY" : "ASCII" ));
        fputs (ctx->ShortBuffer,pFile);
        fputs ("\n",pFile);
        fputs ("DATASET STRUCTURED_POINTS\n",pFile);
        snprintf(ctx->ShortBuffer,80,"DIMENSIONS %d %d %d \n",ctx->uABC[0], ctx->uABC[1], ctx->uABC[2]);
        fputs (ctx->ShortBuffer,pFile);
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

        snprintf(ctx->ShortBuffer,80,"SPACING %.6g %.6g %.6g \n",temp1*a_lattice, temp2*a_lattice, temp3*a_lattice);
        fputs (ctx->ShortBuffer,pFile);

        snprintf(ctx->ShortBuffer,80,"POINT_DATA %d \n",ctx->uABC[0]*ctx->uABC[1]*ctx->uABC[2]);
        fputs (ctx->ShortBuffer,pFile);
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
                        temp0 += VEC_X(S,N);
                        temp1 += VEC_Y(S,N);
                        temp2 += VEC_Z(S,N);
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
                        snprintf(ctx->ShortBuffer,80,"%.6g %.6g %.6g ",temp0, temp1, temp2);
                        fputs (ctx->ShortBuffer,pFile);
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
                        temp0 += VEC_X(S,N);
                        temp1 += VEC_Y(S,N);
                        temp2 += VEC_Z(S,N);
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
                        snprintf(ctx->ShortBuffer,80,"%.6g ",temp0);
                        fputs (ctx->ShortBuffer,pFile);
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
                        temp0 += VEC_X(S,N);
                        temp1 += VEC_Y(S,N);
                        temp2 += VEC_Z(S,N);
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
                        snprintf(ctx->ShortBuffer,80,"%.6g ",temp0);
                        fputs (ctx->ShortBuffer,pFile);
                    }
                }
            }
        }
        fclose (pFile);
        printf("Recording to the file %s is done!\n", vtk_filename);
	} else {
		fprintf(stderr, "Cannot open output file '%s': %s\n", vtk_filename, strerror(errno));
    }
}

void Save_VTK_6(magnoom_ctx *ctx, double* S,
                double* dS, const int mode, const char *vtk_filename)
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
        snprintf(ctx->ShortBuffer,80,"%s\n",(mode==0 ? "BINARY" : "ASCII" ));
        fputs (ctx->ShortBuffer,pFile);
        fputs ("\n",pFile);
        fputs ("DATASET STRUCTURED_POINTS\n",pFile);
        snprintf(ctx->ShortBuffer,80,"DIMENSIONS %d %d %d \n",ctx->uABC[0], ctx->uABC[1], ctx->uABC[2]);
        fputs (ctx->ShortBuffer,pFile);
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

        snprintf(ctx->ShortBuffer,80,"SPACING %.6g %.6g %.6g \n",temp1*a_lattice, temp2*a_lattice, temp3*a_lattice);
        fputs (ctx->ShortBuffer,pFile);

        snprintf(ctx->ShortBuffer,80,"POINT_DATA %d \n",ctx->uABC[0]*ctx->uABC[1]*ctx->uABC[2]);
        fputs (ctx->ShortBuffer,pFile);
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
                        temp0 += VEC_X(S,N);
                        temp1 += VEC_Y(S,N);
                        temp2 += VEC_Z(S,N);
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
                        snprintf(ctx->ShortBuffer,80,"%.6g %.6g %.6g ",temp0, temp1, temp2);
                        fputs (ctx->ShortBuffer,pFile);
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
                        temp0 += VEC_X(dS,N);
                        temp1 += VEC_Y(dS,N);
                        temp2 += VEC_Z(dS,N);
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
                        snprintf(ctx->ShortBuffer,80,"%.6g %.6g %.6g ",temp0, temp1, temp2);
                        fputs (ctx->ShortBuffer,pFile);
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
                        temp0 += VEC_X(dS,N);
                        temp1 += VEC_Y(dS,N);
                        temp2 += VEC_Z(dS,N);
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
                        snprintf(ctx->ShortBuffer,80,"%.6g ",temp0);
                        fputs (ctx->ShortBuffer,pFile);
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
                        temp0 += VEC_X(S,N);
                        temp1 += VEC_Y(S,N);
                        temp2 += VEC_Z(S,N);
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
                        snprintf(ctx->ShortBuffer,80,"%.6g ",temp0);
                        fputs (ctx->ShortBuffer,pFile);
                    }
                }
            }
        }
        fclose (pFile);
        printf("Recording to the file %s is done!\n", vtk_filename);
	} else {
		fprintf(stderr, "Cannot open output file '%s': %s\n", vtk_filename, strerror(errno));
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

void Read_VTK(magnoom_ctx *ctx, double* S, const char *vtk_filename){
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
bool  read_ok = false;
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

                printf("%s has wrong header or wrong file format! (not vtk 2.0 file format):\n", vtk_filename);
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
                                        VEC_X(S,I)=VEC_X(ctx->bS,I)=(double)temp4_x;
                                        VEC_Y(S,I)=VEC_Y(ctx->bS,I)=(double)temp4_y;
                                        VEC_Z(S,I)=VEC_Z(ctx->bS,I)=(double)temp4_z;
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
                                        VEC_X(S,I)=VEC_X(ctx->bS,I)=temp8_x;
                                        VEC_Y(S,I)=VEC_Y(ctx->bS,I)=temp8_y;
                                        VEC_Z(S,I)=VEC_Z(ctx->bS,I)=temp8_z;
                                    }
                                }
                            }
                        }
                    }
                }
                read_ok = true;
            }
        }else{
            printf("%s cannot read data in vtk file!\n", vtk_filename);
        }
        fclose(FilePointer);
        if (read_ok) {
            printf("Reading from the file %s succeeded!\n", vtk_filename);
        } else {
            printf("Reading from the file %s failed!\n", vtk_filename);
        }
    }else{
        fprintf(stderr, "Cannot open input file '%s': %s\n", vtk_filename, strerror(errno));
        printf("Reading from the file %s failed!\n", vtk_filename);
    }
}

#include "linmath.h"		/*All global variables and constants*/

void ReallocateMemoryForImages(magnoom_ctx *ctx, int NumImages, int NOS){
    if ((size_t)NumImages > SIZE_MAX / 3 / (size_t)NOS / sizeof(double)) {
        fprintf(stderr, "ReallocateMemoryForImages: NumImages*NOS*3 would overflow, aborting allocation.\n");
        return;
    }
    size_t count = 3*(size_t)NumImages*(size_t)NOS;
    ctx->Image  = (double *) calloc(count, sizeof(double));
    ctx->dImage = (double *) calloc(count, sizeof(double));
}

#include "math_utils.c"		/*All mathematical fuctions*/
#include "lattice_geometry.c"	/*All functions salculating size and neighbors*/
#include "anisotropy.c"		/*Local and global anisotropy tensor operations*/
#include "solvers.c"		/*CALC THREAD: LLG solvers*/
#include "visualization.c"	/*VISUAL THREAD: All Visualization Functions*/
#include "initial_states.c"	/*Set of functions for initial states*/

#include <fcntl.h> 		/* For O_CREAT and other O_**** constants */

void ReallocateMemoryForSpins(magnoom_ctx *ctx, int NOS){
	if ((size_t)NOS > SIZE_MAX / 3 / sizeof(double)) {
		fprintf(stderr, "ReallocateMemoryForSpins: NOS*3 would overflow, aborting allocation.\n");
		return;
	}
	ctx->S  = (double *)calloc(3*(size_t)NOS, sizeof(double)); // <-- for 10^6 spins allocated memory for ctx->S = 24 Mega Byte
	ctx->bS = (double *)calloc(3*(size_t)NOS, sizeof(double)); // <-- + 24 Mega Byte
	ctx->Kind = (int *)calloc(NOS, sizeof(int));
}

void ReallocateMemoryForAllOther(magnoom_ctx *ctx, int NOS){
	if ((size_t)NOS > SIZE_MAX / 3 / sizeof(double)) {
		fprintf(stderr, "ReallocateMemoryForAllOther: NOS*3 would overflow, aborting allocation.\n");
		return;
	}
	ctx->tS  = (double *)calloc(3*(size_t)NOS, sizeof(double)); // <-- + 24 Mega Byte
	ctx->t2S = (double *)calloc(3*(size_t)NOS, sizeof(double)); // <-- + 24 Mega Byte
	ctx->t3S = (double *)calloc(3*(size_t)NOS, sizeof(double)); // <-- + 24 Mega Byte

	ctx->NoiseX = (float *)calloc(NOS, sizeof(float));
	ctx->NoiseY = (float *)calloc(NOS, sizeof(float));
	ctx->NoiseZ = (float *)calloc(NOS, sizeof(float));	// <-- + 12 Mega Byte

	ctx->PosX = (float *)calloc(NOS, sizeof(float));
	ctx->PosY = (float *)calloc(NOS, sizeof(float));
	ctx->PosZ = (float *)calloc(NOS, sizeof(float));	// <-- + 12 Mega Byte

	ctx->BlockPosX = (float *)calloc(ctx->NOB, sizeof(float));
	ctx->BlockPosY = (float *)calloc(ctx->NOB, sizeof(float));
	ctx->BlockPosZ = (float *)calloc(ctx->NOB, sizeof(float));	// <-- + 12 Mega Byte

	ctx->HeffX = (double *)calloc(NOS, sizeof(double));
	ctx->HeffY = (double *)calloc(NOS, sizeof(double));
	ctx->HeffZ = (double *)calloc(NOS, sizeof(double));// <-- + 12 Mega Byte

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

// Width of the table column that prints `count` integers under `header`:
// the longest printed value (including a leading '-' when needed), so every
// value lines up no matter how large or negative it gets, but never
// narrower than the header text itself.
static int IntColumnWidth(const char *header, const int *values, int count)
{
	int width = (int)strlen(header);
	for (int i = 0; i < count; i++) {
		char buffer[16];
		int length = snprintf(buffer, sizeof(buffer), "%d", values[i]);
		if (length > width) width = length;
	}
	return width;
}

// Longest printed width (including a leading '-' when needed) across
// `count` floats, each formatted with `format` exactly as it will be
// displayed - used to size a Markdown table's numeric column so every
// value in it lines up, whether or not a '-' sign appears. Never returns
// less than 1, so every column has room for at least one character.
static int MaxPrintedFloatWidth(const float *values, int count, const char *format)
{
	int width = 1;
	for (int i = 0; i < count; i++) {
		char buffer[32];
		int length = snprintf(buffer, sizeof(buffer), format, (double)values[i]);
		if (length > width) width = length;
	}
	return width;
}

// Positions of the two atoms of neighbor pair `pair` -- the central atom in
// the home block, the neighbor in its own image block -- and the distance
// between them. None of the three is stored in the map, so both the
// column-measuring pass and the printing pass below derive them here.
static float NeighborPairGeometry(magnoom_ctx *ctx, int pair, float AXYZ[3], float NXYZ[3])
{
	const NeighborMap *map = &ctx->Neighbors;
	GetPosition(ctx->abc, ctx->Block, map->AIdxBlock[pair], 0, 0, 0, AXYZ);
	GetPosition(ctx->abc, ctx->Block, map->NIdxBlock[pair],
		map->NIdxGridA[pair], map->NIdxGridB[pair], map->NIdxGridC[pair], NXYZ);
	float dx = AXYZ[0] - NXYZ[0], dy = AXYZ[1] - NXYZ[1], dz = AXYZ[2] - NXYZ[2];
	return sqrtf(dx*dx + dy*dy + dz*dz);
}

// Prints the lattice/neighbor-map diagnostics gathered during startup as
// a sequence of Markdown tables (unit cell, coordination shells,
// neighbors-per-atom, a parameter summary, and the full neighbor-pair
// map), so the output can be pasted directly into a Markdown document.
// Every numeric column is sized to the widest value actually present
// (including a '-' sign where one occurs) and to its own header text, so
// the raw text also stays aligned in a monospace terminal.
//
// Every cell is printed as "| " + a field of `width` characters + " |",
// so between two pipes a separator cell must span width + 2 characters:
// a right-aligned column prints width + 1 dashes and spends the last
// character on Markdown's ':' marker, a left-aligned one prints all
// width + 2 as dashes.
void PrintLatticeInitReport(magnoom_ctx *ctx)
{
	// Sliced with "%.*s" to emit a separator cell's dashes in one go. Far
	// longer than any column width these tables can reach, since those come
	// from "%d"/"%.6g"/"%.3f" renderings of lattice data.
	static const char dashes[] = "----------------------------------------------------------------";

	const NeighborMap *map = &ctx->Neighbors;
	int cellWidth = MaxPrintedFloatWidth((const float *)ctx->abc, 9, "%.6g");
	int blockWidth = MaxPrintedFloatWidth((const float *)ctx->Block, 3 * ctx->AtomsPerBlock, "%.6g");
	if (blockWidth > cellWidth) cellWidth = blockWidth;

	printf("\n### Translation vectors\n\n");
	const char *vectorNames[3] = { "a", "b", "c" };
	int vectorWidth = (int)strlen("Vector"); // wider than any of the "a"/"b"/"c" rows
	printf("| %-*s | %*s | %*s | %*s |\n", vectorWidth, "Vector", cellWidth, "X", cellWidth, "Y", cellWidth, "Z");
	printf("|%.*s|%.*s:|%.*s:|%.*s:|\n", vectorWidth + 2, dashes,
		cellWidth + 1, dashes, cellWidth + 1, dashes, cellWidth + 1, dashes);
	for (int v = 0; v < 3; v++) {
		printf("| %-*s | %*.6g | %*.6g | %*.6g |\n", vectorWidth, vectorNames[v],
			cellWidth, (double)ctx->abc[v][0], cellWidth, (double)ctx->abc[v][1], cellWidth, (double)ctx->abc[v][2]);
	}

	printf("\n### Atom positions in the unit cell\n\n");
	int atomIndexMax = ctx->AtomsPerBlock > 0 ? ctx->AtomsPerBlock - 1 : 0;
	int atomWidth = IntColumnWidth("Atom", &atomIndexMax, 1);
	printf("| %*s | %*s | %*s | %*s |\n", atomWidth, "Atom", cellWidth, "X", cellWidth, "Y", cellWidth, "Z");
	printf("|%.*s:|%.*s:|%.*s:|%.*s:|\n", atomWidth + 1, dashes,
		cellWidth + 1, dashes, cellWidth + 1, dashes, cellWidth + 1, dashes);
	for (int i = 0; i < ctx->AtomsPerBlock; i++) {
		printf("| %*d | %*.6g | %*.6g | %*.6g |\n", atomWidth, i,
			cellWidth, (double)ctx->Block[i][0], cellWidth, (double)ctx->Block[i][1], cellWidth, (double)ctx->Block[i][2]);
	}

	printf("\n### Coordination shells\n\n");
	int shellIndexMax = map->ShellCount > 0 ? map->ShellCount - 1 : 0;
	int shellWidth = IntColumnWidth("Shell", &shellIndexMax, 1);
	int radiusWidth = MaxPrintedFloatWidth(map->ShellRadius, map->ShellCount, "%f"); // always at least as wide as the "R" header
	printf("| %*s | %*s |\n", shellWidth, "Shell", radiusWidth, "R");
	printf("|%.*s:|%.*s:|\n", shellWidth + 1, dashes, radiusWidth + 1, dashes);
	for (int i = 0; i < map->ShellCount; i++) {
		printf("| %*d | %*f |\n", shellWidth, i, radiusWidth, map->ShellRadius[i]);
	}

	printf("\n### Neighbors per atom\n\n");
	int neighborsPerAtom[MAX_ATOMS_PER_BLOCK] = {0};
	for (int i = 0; i < map->PairCount; i++) {
		neighborsPerAtom[map->AIdxBlock[i]]++;
	}
	int neighborsWidth = IntColumnWidth("Neighbors", neighborsPerAtom, ctx->AtomsPerBlock);
	printf("| %*s | %*s |\n", atomWidth, "Atom", neighborsWidth, "Neighbors");
	printf("|%.*s:|%.*s:|\n", atomWidth + 1, dashes, neighborsWidth + 1, dashes);
	for (int i = 0; i < ctx->AtomsPerBlock; i++) {
		printf("| %*d | %*d |\n", atomWidth, i, neighborsWidth, neighborsPerAtom[i]);
	}

	printf("\n### Summary\n\n");
	const char *summaryNames[3] = { "AtomsPerBlock", "ShellNumber", "NeighborPairs" };
	int summaryValues[3] = { ctx->AtomsPerBlock, map->ShellCount, map->PairCount };
	int paramWidth = (int)strlen("Parameter");
	for (int i = 0; i < 3; i++) {
		int nameWidth = (int)strlen(summaryNames[i]);
		if (nameWidth > paramWidth) paramWidth = nameWidth;
	}
	int summaryWidth = IntColumnWidth("Value", summaryValues, 3);
	printf("| %-*s | %*s |\n", paramWidth, "Parameter", summaryWidth, "Value");
	printf("|%.*s|%.*s:|\n", paramWidth + 2, dashes, summaryWidth + 1, dashes);
	for (int i = 0; i < 3; i++) {
		printf("| %-*s | %*d |\n", paramWidth, summaryNames[i], summaryWidth, summaryValues[i]);
	}

	printf("\n### Neighbor pairs\n\n");
	// The seven integer columns, in printed order. "N" is the row number
	// rather than a stored array, so it is measured from the last index the
	// table will print; the other six are read straight from the pair map.
	const char *intHeaders[7] = { "N", "Ia", "In", "Jn", "Kn", "Ln", "Sn" };
	int lastIndex = map->PairCount > 0 ? map->PairCount - 1 : 0;
	const int *intColumns[7] = { &lastIndex, map->AIdxBlock, map->NIdxBlock,
		map->NIdxGridA, map->NIdxGridB, map->NIdxGridC, map->ShellIdx };
	int intWidths[7];
	for (int c = 0; c < 7; c++) {
		int count = (c == 0) ? 1 : map->PairCount;
		intWidths[c] = IntColumnWidth(intHeaders[c], intColumns[c], count);
	}

	// All twelve float columns share one width so the block reads as a grid.
	// Jij/Dij/D_x/D_y/D_z are stored and can be measured directly, but
	// Ax/Ay/Az/Nx/Ny/Nz/Dist are derived per row, so a first pass repeats
	// the same NeighborPairGeometry() the printing pass below performs.
	// "%.3f" is never narrower than any of these columns' headers.
	const char *floatHeaders[12] = { "Ax", "Ay", "Az", "Nx", "Ny", "Nz", "Dist", "Jij", "Dij", "D_x", "D_y", "D_z" };
	const float *floatColumns[5] = { ctx->Jexc, ctx->Dexc, ctx->VDMX, ctx->VDMY, ctx->VDMZ };
	int floatWidth = 1;
	for (int c = 0; c < 5; c++) {
		int width = MaxPrintedFloatWidth(floatColumns[c], map->PairCount, "%.3f");
		if (width > floatWidth) floatWidth = width;
	}
	for (int i = 0; i < map->PairCount; i++) {
		float AXYZ[3], NXYZ[3]; // atom / neighbor position
		float distance = NeighborPairGeometry(ctx, i, AXYZ, NXYZ);
		float rowValues[7] = { AXYZ[0], AXYZ[1], AXYZ[2], NXYZ[0], NXYZ[1], NXYZ[2], distance };
		int rowWidth = MaxPrintedFloatWidth(rowValues, 7, "%.3f");
		if (rowWidth > floatWidth) floatWidth = rowWidth;
	}

	printf("|");
	for (int c = 0; c < 7; c++) printf(" %*s |", intWidths[c], intHeaders[c]);
	for (int c = 0; c < 12; c++) printf(" %*s |", floatWidth, floatHeaders[c]);
	putchar('\n');
	printf("|");
	for (int c = 0; c < 7; c++) printf("%.*s:|", intWidths[c] + 1, dashes);
	for (int c = 0; c < 12; c++) printf("%.*s:|", floatWidth + 1, dashes);
	putchar('\n');

	for (int i = 0; i < map->PairCount; i++) {
		printf("| %*d |", intWidths[0], i); // the "N" column is the row number
		for (int c = 1; c < 7; c++) printf(" %*d |", intWidths[c], intColumns[c][i]);
		float AXYZ[3], NXYZ[3]; // atom / neighbor position
		float distance = NeighborPairGeometry(ctx, i, AXYZ, NXYZ);
		printf(" %*.3f | %*.3f | %*.3f |", floatWidth, AXYZ[0], floatWidth, AXYZ[1], floatWidth, AXYZ[2]);
		printf(" %*.3f | %*.3f | %*.3f |", floatWidth, NXYZ[0], floatWidth, NXYZ[1], floatWidth, NXYZ[2]);
		printf(" %*.3f | %*.3f | %*.3f | %*.3f | %*.3f | %*.3f |\n",
			floatWidth, distance,
			floatWidth, ctx->Jexc[i], floatWidth, ctx->Dexc[i],
			floatWidth, ctx->VDMX[i], floatWidth, ctx->VDMY[i], floatWidth, ctx->VDMZ[i]);
	}
}

/*************************************************************************/
/*                        Program Main Thread                            */
/*************************************************************************/

#ifndef MAGNOOM_NO_MAIN
int
main (int argc, char **argv){
	magnoom_ctx mag_ctx = {0};
	char executable_directory[MAGNOOM_PATH_CAPACITY];
	(void)argc;
	if (!magnoom_ctx_init(&mag_ctx)) {
		fprintf(stderr, "Unable to initialize Magnoom data.\n");
		return 1;
	}
	if (!magnoom_executable_directory(executable_directory,
		sizeof(executable_directory), argv != NULL ? argv[0] : NULL)) {
		fprintf(stderr, "Unable to determine the executable directory.\n");
		return 1;
	}
	magnoom_copy_path(mag_ctx.input_directory,
		sizeof(mag_ctx.input_directory), executable_directory);
	magnoom_copy_path(mag_ctx.output_directory,
		sizeof(mag_ctx.output_directory), executable_directory);
		
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
	if (!readConfigFile(&mag_ctx)) {
		fprintf(stderr, "Unable to apply magnoom.cfg.\n");
		free(mag_ctx.Block);
		return 1;
	}
	////////////////////////////////////////////////

	if (!anisotropy_apply_config_records(&mag_ctx)) {
		fprintf(stderr, "Unable to initialize anisotropy tensors.\n");
		return 1;
	}
	if (!magnoom_ctx_build_neighbor_map(&mag_ctx)) {
		fprintf(stderr, "Unable to build the neighbor map.\n");
		return 1;
	}
	mag_ctx.Jexc = (float *)calloc(mag_ctx.Neighbors.PairCount, sizeof(float));
	mag_ctx.Bexc = (float *)calloc(mag_ctx.Neighbors.PairCount, sizeof(float));
	mag_ctx.Dexc = (float *)calloc(mag_ctx.Neighbors.PairCount, sizeof(float));
	mag_ctx.VDMX = (float *)calloc(mag_ctx.Neighbors.PairCount, sizeof(float));
	mag_ctx.VDMY = (float *)calloc(mag_ctx.Neighbors.PairCount, sizeof(float));
	mag_ctx.VDMZ = (float *)calloc(mag_ctx.Neighbors.PairCount, sizeof(float));
	SetExch1( mag_ctx.abc, mag_ctx.Block, mag_ctx.Neighbors.PairCount, mag_ctx.Jij, mag_ctx.Bij, mag_ctx.Dij,
		mag_ctx.Neighbors.AIdxBlock, mag_ctx.Neighbors.NIdxBlock, mag_ctx.Neighbors.NIdxGridA, mag_ctx.Neighbors.NIdxGridB, mag_ctx.Neighbors.NIdxGridC, mag_ctx.Neighbors.ShellIdx,
		mag_ctx.Jexc, mag_ctx.Bexc, mag_ctx.Dexc, mag_ctx.VDMX, mag_ctx.VDMY, mag_ctx.VDMZ);
	PrintLatticeInitReport(&mag_ctx);

	// Open output file:
	pthread_mutex_init(&mag_ctx.record_mutex,0);
	magnoom_reset_record_file(&mag_ctx);

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

	GetBox(&mag_ctx);
	UpdateSpinPositions(&mag_ctx);
	UpdateKind(&mag_ctx);
	InitSpinComponents( &mag_ctx, mag_ctx.PosX, mag_ctx.PosY, mag_ctx.PosZ, mag_ctx.S, 0);
	for (int i=0;i<mag_ctx.NOS;i++) { VEC_X(mag_ctx.bS,i)=VEC_X(mag_ctx.S,i); VEC_Y(mag_ctx.bS,i)=VEC_Y(mag_ctx.S,i); VEC_Z(mag_ctx.bS,i)=VEC_Z(mag_ctx.S,i);}

    // Set OpenGL context initial state.
	setupOpenGL(&mag_ctx);
    //  Allocate memory for vetices, mag_ctx.normals, mag_ctx.colors and indicies array used in drawing subrutines
	ReallocateArrayDrawing(&mag_ctx);
	// Fill array for prototype (arrow or cane) array 
	UpdateProtoVerNorInd_Spins(&mag_ctx);
	// Fill big array for indecies for all arrows, cans, cones or boxes 
	UpdateIndices(&mag_ctx);
	UpdateVerticesNormalsColors(&mag_ctx);
	InitVBOMesh(&mag_ctx.spin_mesh, GL_DYNAMIC_DRAW);
	mag_ctx.spin_mesh.uses_normals = (mag_ctx.WhichVectorMode == BOX1 || mag_ctx.WhichVectorMode == ARROW1 || mag_ctx.WhichVectorMode == CONE1) ? 1 : 0;
	mag_ctx.spin_mesh.component_count = mag_ctx.VCNum;
	mag_ctx.spin_mesh.component_capacity = mag_ctx.VCNum;
	mag_ctx.spin_mesh.index_count = mag_ctx.IdNum;
	mag_ctx.spin_mesh.index_capacity = mag_ctx.IdNum;
	CreateVBOMesh(&mag_ctx.spin_mesh);
	UploadVBOMesh(&mag_ctx.spin_mesh, mag_ctx.vertices, mag_ctx.normals, mag_ctx.colors, mag_ctx.indices, VBO_UPLOAD_ALL);

    mag_ctx.BextDCDirection[0]=sin(PI*mag_ctx.BextDCTheta/180)*cos(PI*mag_ctx.BextDCPhi/180);
	mag_ctx.BextDCDirection[1]=sin(PI*mag_ctx.BextDCTheta/180)*sin(PI*mag_ctx.BextDCPhi/180);
	mag_ctx.BextDCDirection[2]=cos(PI*mag_ctx.BextDCTheta/180);

	ReallocateArrayDrawing_BextDC(&mag_ctx);
	UpdateProtoVerNorInd_BextDC(&mag_ctx);
    UpdateVerticesNormalsColors_BextDC(&mag_ctx);
	InitVBOMesh(&mag_ctx.BextDC_mesh, GL_STATIC_DRAW);
	mag_ctx.BextDC_mesh.uses_normals = 1;
	mag_ctx.BextDC_mesh.component_count = mag_ctx.VCNum_BextDC;
	mag_ctx.BextDC_mesh.component_capacity = mag_ctx.VCNum_BextDC;
	mag_ctx.BextDC_mesh.index_count = mag_ctx.IdNum_BextDC;
	mag_ctx.BextDC_mesh.index_capacity = mag_ctx.IdNum_BextDC;
	CreateVBOMesh(&mag_ctx.BextDC_mesh);
    UploadVBOMesh(&mag_ctx.BextDC_mesh, mag_ctx.vertices_BextDC, mag_ctx.normals_BextDC, mag_ctx.colors_BextDC, mag_ctx.indices_BextDC, VBO_UPLOAD_ALL);

	/* --- BOX mesh: static, startup-only --- */
	ReallocateArrayDrawing_BOX(&mag_ctx);
	UpdateVerticesNormalsColors_BOX(&mag_ctx);
	InitVBOMesh(&mag_ctx.box_mesh, GL_STATIC_DRAW);
	mag_ctx.box_mesh.uses_normals = 1;
	mag_ctx.box_mesh.component_count = mag_ctx.VCNum_BOX;
	mag_ctx.box_mesh.component_capacity = mag_ctx.VCNum_BOX;
	mag_ctx.box_mesh.index_count = mag_ctx.IdNum_BOX;
	mag_ctx.box_mesh.index_capacity = mag_ctx.IdNum_BOX;
	CreateVBOMesh(&mag_ctx.box_mesh);
	UploadVBOMesh(&mag_ctx.box_mesh, mag_ctx.vertices_BOX, mag_ctx.normals_BOX, mag_ctx.colors_BOX, mag_ctx.indices_BOX, VBO_UPLOAD_ALL);

	/* --- BASIS mesh: static, startup-only --- */
	ReallocateArrayDrawing_BASIS(&mag_ctx);
	UpdateVerticesNormalsColors_BASIS(&mag_ctx);
	InitVBOMesh(&mag_ctx.basis_mesh, GL_STATIC_DRAW);
	mag_ctx.basis_mesh.uses_normals = 1;
	mag_ctx.basis_mesh.component_count = mag_ctx.VCNum_BASIS;
	mag_ctx.basis_mesh.component_capacity = mag_ctx.VCNum_BASIS;
	mag_ctx.basis_mesh.index_count = mag_ctx.IdNum_BASIS;
	mag_ctx.basis_mesh.index_capacity = mag_ctx.IdNum_BASIS;
	CreateVBOMesh(&mag_ctx.basis_mesh);
	UploadVBOMesh(&mag_ctx.basis_mesh, mag_ctx.vertices_BASIS, mag_ctx.normals_BASIS, mag_ctx.colors_BASIS, mag_ctx.indices_BASIS, VBO_UPLOAD_ALL);

	/* --- PBC meshes: three static axis meshes, startup-only --- */
	ReallocateArrayDrawing_PBC(&mag_ctx);
	InitVBOMesh(&mag_ctx.pbc_mesh[0], GL_STATIC_DRAW);
	mag_ctx.pbc_mesh[0].uses_normals = 1;
	mag_ctx.pbc_mesh[0].component_count = mag_ctx.VCNum_PBC;
	mag_ctx.pbc_mesh[0].component_capacity = mag_ctx.VCNum_PBC;
	mag_ctx.pbc_mesh[0].index_count = mag_ctx.IdNum_PBC;
	mag_ctx.pbc_mesh[0].index_capacity = mag_ctx.IdNum_PBC;
	InitVBOMesh(&mag_ctx.pbc_mesh[1], GL_STATIC_DRAW);
	mag_ctx.pbc_mesh[1].uses_normals = 1;
	mag_ctx.pbc_mesh[1].component_count = mag_ctx.VCNum_PBC;
	mag_ctx.pbc_mesh[1].component_capacity = mag_ctx.VCNum_PBC;
	mag_ctx.pbc_mesh[1].index_count = mag_ctx.IdNum_PBC;
	mag_ctx.pbc_mesh[1].index_capacity = mag_ctx.IdNum_PBC;
	InitVBOMesh(&mag_ctx.pbc_mesh[2], GL_STATIC_DRAW);
	mag_ctx.pbc_mesh[2].uses_normals = 1;
	mag_ctx.pbc_mesh[2].component_count = mag_ctx.VCNum_PBC;
	mag_ctx.pbc_mesh[2].component_capacity = mag_ctx.VCNum_PBC;
	mag_ctx.pbc_mesh[2].index_count = mag_ctx.IdNum_PBC;
	mag_ctx.pbc_mesh[2].index_capacity = mag_ctx.IdNum_PBC;
	UpdateVerticesNormalsColors_PBC(&mag_ctx, 0);
	UpdateVerticesNormalsColors_PBC(&mag_ctx, 1);
	UpdateVerticesNormalsColors_PBC(&mag_ctx, 2);
	CreateVBOMesh(&mag_ctx.pbc_mesh[0]);
	UploadVBOMesh(&mag_ctx.pbc_mesh[0], mag_ctx.vertices_PBC_A, mag_ctx.normals_PBC_A, mag_ctx.colors_PBC_A, mag_ctx.indices_PBC_A, VBO_UPLOAD_ALL);
	CreateVBOMesh(&mag_ctx.pbc_mesh[1]);
	UploadVBOMesh(&mag_ctx.pbc_mesh[1], mag_ctx.vertices_PBC_B, mag_ctx.normals_PBC_B, mag_ctx.colors_PBC_B, mag_ctx.indices_PBC_B, VBO_UPLOAD_ALL);
	CreateVBOMesh(&mag_ctx.pbc_mesh[2]);
	UploadVBOMesh(&mag_ctx.pbc_mesh[2], mag_ctx.vertices_PBC_C, mag_ctx.normals_PBC_C, mag_ctx.colors_PBC_C, mag_ctx.indices_PBC_C, VBO_UPLOAD_ALL);

	// Explicit GLFW loop, matching the AntTweakBar Legacy examples.
	while (glfwGetWindowParam(GLFW_OPENED) && !mag_ctx.WindowShouldClose) {
		idle(&mag_ctx);
		Display(&mag_ctx);
		glfwSwapBuffers();
		glfwPollEvents();
		magnoom_service_modal(&mag_ctx);
	}

	mag_ctx.EngineShutdownRequested = true;
	mag_ctx.EngineRunState = DO_IT;
	for (int i = 0; i < THREADS_NUMBER; ++i) (void)pthread_join(thread_id[i], NULL);
	// Destroy GPU buffers (needs active GL context):
	DestroyVBOMesh(&mag_ctx.spin_mesh);
	DestroyVBOMesh(&mag_ctx.BextDC_mesh);
	DestroyVBOMesh(&mag_ctx.box_mesh);
	DestroyVBOMesh(&mag_ctx.basis_mesh);
	DestroyVBOMesh(&mag_ctx.pbc_mesh[0]);
	DestroyVBOMesh(&mag_ctx.pbc_mesh[1]);
	DestroyVBOMesh(&mag_ctx.pbc_mesh[2]);
	TwTerminate();
	glfwTerminate();

	// Deallocate memory before closing the program:
	free(mag_ctx.Neighbors.AIdxBlock);
	free(mag_ctx.Neighbors.NIdxBlock);
	free(mag_ctx.Neighbors.NIdxGridA);
	free(mag_ctx.Neighbors.NIdxGridB);
	free(mag_ctx.Neighbors.NIdxGridC);
	free(mag_ctx.Neighbors.ShellIdx);

	free(mag_ctx.Jexc);  			free(mag_ctx.Bexc);  			free(mag_ctx.Dexc);
	free(mag_ctx.VDMX);  			free(mag_ctx.VDMY);  			free(mag_ctx.VDMZ);

	free(mag_ctx.S);     			free(mag_ctx.bS);
	free(mag_ctx.tS);    			free(mag_ctx.t2S);   			free(mag_ctx.t3S);
	/* mag_ctx.Image/dImage intentionally not freed here -- pre-existing     */
	/* (unfixed) leak, preserved as-is.                                      */
	free(mag_ctx.HeffX); 			free(mag_ctx.HeffY); 			free(mag_ctx.HeffZ);
	free(mag_ctx.NoiseX);   			free(mag_ctx.NoiseY);   			free(mag_ctx.NoiseZ);
	free(mag_ctx.PosX);    			free(mag_ctx.PosY);    			free(mag_ctx.PosZ);
	free(mag_ctx.BlockPosX);   			free(mag_ctx.BlockPosY);   			free(mag_ctx.BlockPosZ);
	free(mag_ctx.RHue);  			free(mag_ctx.GHue);  			free(mag_ctx.BHue);
	free(mag_ctx.vertices);			free(mag_ctx.vertices_BextDC); 		free(mag_ctx.vertices_BOX); 	free(mag_ctx.vertices_PBC_A);
	free(mag_ctx.normals);			free(mag_ctx.normals_BextDC); 		free(mag_ctx.normals_BOX); 		free(mag_ctx.normals_PBC_A);
	free(mag_ctx.colors);			free(mag_ctx.colors_BextDC); 		free(mag_ctx.colors_BOX); 		free(mag_ctx.colors_PBC_A);
	free(mag_ctx.indices);			free(mag_ctx.indices_BextDC); 		free(mag_ctx.indices_BOX); 		free(mag_ctx.indices_PBC_A);
	free(mag_ctx.vertexProto);		free(mag_ctx.vertexProto_BextDC);
	free(mag_ctx.normalProto);		free(mag_ctx.normalProto_BextDC);
	free(mag_ctx.indicesProto);		free(mag_ctx.indicesProto_BextDC);

	free(mag_ctx.Block);

	free(mag_ctx.Proj);

	magnoom_close_record_file(&mag_ctx);
	pthread_mutex_destroy(&mag_ctx.record_mutex);

	return 0;
}
#endif
