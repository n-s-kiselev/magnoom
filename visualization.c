enum WindowMouseButton { WINDOW_MOUSE_LEFT, WINDOW_MOUSE_MIDDLE, WINDOW_MOUSE_RIGHT };
enum WindowButtonState { WINDOW_BUTTON_DOWN, WINDOW_BUTTON_UP };
enum WindowSpecialKey {
	WINDOW_KEY_UP = 1, WINDOW_KEY_DOWN, WINDOW_KEY_F1, WINDOW_KEY_F2,
	WINDOW_KEY_F3, WINDOW_KEY_F4, WINDOW_KEY_F5, WINDOW_KEY_F6, WINDOW_KEY_F12
};

// which button:
typedef enum {XUP, YUP, ZUP, ADD_SET, RESET, QUIT, PLAY, RECORD} enButton;

void ChangeVectorMode( magnoom_ctx *, int );
void ChangeColorMap( magnoom_ctx *, int );
void ChangeInitialState( magnoom_ctx * );
void ExecuteCommand( magnoom_ctx *, int );
void HandleKeyDown( magnoom_ctx *, unsigned char );
void HandleKeyUp( magnoom_ctx *, unsigned char );
void HandleSpecialKey( magnoom_ctx *, int );
void HandleMouseButton( magnoom_ctx *, int, int, int, int );
void HandleMouseDrag( magnoom_ctx *, int, int );
float ElapsedSeconds( );
// color control functions
void HSVtoRGB(magnoom_ctx *, float[3], float [3] , int, int);
void InitRGB(float* , float* , float* , int*);
// VBO mesh descriptor helpers
void InitVBOMesh(vbo_mesh *, GLenum);
void CreateVBOMesh(vbo_mesh *);
void UploadVBOMesh(vbo_mesh *, const float *, const float *, const float *, const GLuint *, unsigned int);
void DrawVBOMeshIndexed(const vbo_mesh *, GLenum, GLsizei, GLsizei);
void DestroyVBOMesh(vbo_mesh *);
// VBO array preparing functions
void ReallocateArrayDrawing(magnoom_ctx *);
void UpdateProtoVerNorInd_Spins(magnoom_ctx *);
void UpdateProtoVerNorInd_BextDC(magnoom_ctx *);

void UpdateSpinComponents(float * , float * , float * , int);
void UpdateSpinPositions(magnoom_ctx *);
void InitSpinComponents(magnoom_ctx *, float * , float * , float * , double * , int N);
void UpdateIndices(magnoom_ctx *);
void UpdateVerticesNormalsColors(magnoom_ctx *);
void UpdateVerticesNormalsColors_BextDC(magnoom_ctx *);
void UpdateKind(magnoom_ctx *);
// drawing functions
void GetBox(magnoom_ctx *);
void drawVBO(magnoom_ctx *);


void idle(magnoom_ctx *);
void setupTweakBar(magnoom_ctx *);
void GLFWKeyCallback(int, int);
void GLFWCharCallback(int, int);
void GLFWMouseButtonCallback(int, int);
void GLFWCursorPosCallback(int, int);
void GLFWScrollCallback(int);
void GLFWWindowSizeCallback(int, int);

// GLFW2's callbacks take no user-data parameter (there is no window handle
// to attach one to, unlike GLFW3) - this is the one deliberate global
// needed so the callbacks below can reach the program's context.
static magnoom_ctx *g_ctx;



// return the number of seconds since the start of the program:
float ElapsedSeconds( )	{
	return (float)glfwGetTime();
}

// GLFW2 has no notion of a separate framebuffer size (no HiDPI awareness):
// the window's real backing-store size is read straight from the OpenGL
// viewport and compared to the window's logical size, matching
// TwSimpleGLFW.c's ComputeHiDPIScale(). On a standard (non-HiDPI) display
// this is 1.0.
//
// Measured more than once during startup rather than a single time right
// after glfwOpenWindow() returns: Cocoa has not necessarily finished
// establishing the window's real Retina backing store at that exact
// moment, so glGetIntegerv(GL_VIEWPORT) can still report the pre-Retina
// (1x) size then - a single measurement taken too early can permanently
// latch onto that wrong value. Matching this project's former shim
// (vendor/glfw2/TwGLFW2.h's tw_glfw2_measure_scale()), measuring is
// retried a couple of times while window setup is still settling, then
// permanently locked in (see contentScaleLocked below and setupOpenGL()'s
// use of it) so that GLFWWindowSizeCallback's later, real resize calls
// reuse that fixed ratio instead of re-deriving it from glGetIntegerv
// (GL_VIEWPORT): that query only ever reflects whatever this file itself
// last passed to glViewport() (ApplyFramebufferSize()), one resize
// callback in arrears, so re-measuring on every resize is self-referential
// and freezes the real viewport at its first size forever.
static int contentScaleLocked = 0;
#ifdef _WIN32
// GLFW2's Win32 backend (vendor/glfw2/lib/win32/) has no HiDPI awareness of
// its own, and this process never declared itself DPI-aware, so Windows
// handed it a virtualized, always-96-DPI coordinate space and stretched the
// final image to fit the real screen -- blurry, but glfwGetWindowSize() and
// glGetIntegerv(GL_VIEWPORT) below always agreed in that virtualized space
// regardless, so ContentScaleX/Y still came out at the correct 1.0 (there
// was never a second, framebuffer-vs-window discrepancy for this file's
// scale-ratio comparison to detect or correct here, unlike Cocoa's real
// points-vs-pixels split below -- see MeasureContentScale()).
// EnableWindowsDpiAwareness() (called once, before any window is created --
// see setupOpenGL()) opts out of that virtualization so the window is
// created, and rendered, at its real physical pixel size instead.
//
// Resolves its Win32 API dynamically via GetProcAddress rather than linking
// directly or depending on the toolchain's headers declaring it
// (SetProcessDpiAwarenessContext only exists on Windows 10 1703+) -- this
// avoids any new compile-time dependency and degrades gracefully (falls
// back to the older, Vista+ SetProcessDPIAware, or does nothing) rather than
// failing to build or link. The DPI_AWARENESS_CONTEXT type and the
// PER_MONITOR_AWARE_V2 constant are defined locally at their documented,
// stable values instead of assuming a specific header declares them.
#ifndef DPI_AWARENESS_CONTEXT_PER_MONITOR_AWARE_V2
typedef void *MAGNOOM_DPI_AWARENESS_CONTEXT;
#define DPI_AWARENESS_CONTEXT_PER_MONITOR_AWARE_V2 ((MAGNOOM_DPI_AWARENESS_CONTEXT)(-4))
#else
typedef DPI_AWARENESS_CONTEXT MAGNOOM_DPI_AWARENESS_CONTEXT;
#endif

static void EnableWindowsDpiAwareness(void)
{
	HMODULE user32 = GetModuleHandleA("user32.dll");
	if (user32 == NULL) return;

	typedef BOOL (WINAPI *SetProcessDpiAwarenessContextProc)(MAGNOOM_DPI_AWARENESS_CONTEXT);
	SetProcessDpiAwarenessContextProc setContext =
		(SetProcessDpiAwarenessContextProc)GetProcAddress(user32, "SetProcessDpiAwarenessContext");
	if (setContext != NULL && setContext(DPI_AWARENESS_CONTEXT_PER_MONITOR_AWARE_V2)) return;

	typedef BOOL (WINAPI *SetProcessDPIAwareProc)(void);
	SetProcessDPIAwareProc setAware =
		(SetProcessDPIAwareProc)GetProcAddress(user32, "SetProcessDPIAware");
	if (setAware != NULL) setAware();
}
#endif

// GLFW2 has no notion of a separate framebuffer size on any platform (no
// HiDPI awareness of its own): the window's real backing-store size is read
// straight from the OpenGL viewport and compared to the window's logical
// size, matching TwSimpleGLFW.c's ComputeHiDPIScale(). This only ever finds
// a real discrepancy on macOS, where Cocoa keeps a genuine points-vs-pixels
// split even for a window created through this old, Retina-unaware GLFW2
// (see setupOpenGL()'s DPI comment for Windows; GLFW2's X11 backend has no
// separate framebuffer concept either, so this is 1.0 there too). On a
// standard (non-HiDPI) display, or on Windows/Linux, this is 1.0.
static void MeasureContentScale(magnoom_ctx *ctx)
{
	if (contentScaleLocked) return;
	int width = 1, height = 1;
	GLint viewport[4];
	glfwGetWindowSize(&width, &height);
	glGetIntegerv(GL_VIEWPORT, viewport);
	if (width > 0 && height > 0 && viewport[2] > 0 && viewport[3] > 0) {
		ctx->ContentScaleX = (double)viewport[2] / width;
		ctx->ContentScaleY = (double)viewport[3] / height;
	}
}

void ApplyFramebufferSize( magnoom_ctx *ctx, int window_width, int window_height)
{
	ctx->asp_rat     = (float)((   ((double)window_width)/((double)window_height)   ));
	ctx->asp_rat_inv = (float)((   ((double)window_height)/((double)window_width)   ));

	glViewport(0, 0, window_width, window_height);
	    // Send the new window size to AntTweakBar
    TwWindowSize(window_width, window_height);
}

void Zup( magnoom_ctx *ctx )
{
    ctx->q_Rotation[0]=ctx->q_Rotation[1]=ctx->q_Rotation[2]=0.;
    ctx->q_Rotation[3]=1.;
    ctx->Rot[0] = ctx->Rot[1] = ctx->Rot[2] = ctx->TransXYZ[0] = ctx->TransXYZ[1] = 0.;
}

void Xup( magnoom_ctx *ctx )
{   float quat1[4];
    Zup(ctx);
    SetQuaternionFromAxisAngle(ctx->axisZ, -D2R*90, quat1);
    MultiplyQuaternions( quat1, ctx->q_Rotation, ctx->q_Rotation);
    SetQuaternionFromAxisAngle(ctx->axisX, -D2R*90, quat1);
    MultiplyQuaternions( quat1, ctx->q_Rotation, ctx->q_Rotation);
    GetEulerFromQuaternion(ctx->q_Rotation, ctx->Rot);
    ctx->TransXYZ[0] = ctx->TransXYZ[1] = 0.;
}

void Yup( magnoom_ctx *ctx )//metka
{   float quat1[4];
    Zup(ctx);
    SetQuaternionFromAxisAngle(ctx->axisZ, D2R*180, quat1);
    MultiplyQuaternions( quat1, ctx->q_Rotation, ctx->q_Rotation);
    SetQuaternionFromAxisAngle(ctx->axisX, -D2R*90, quat1);
    MultiplyQuaternions( quat1, ctx->q_Rotation, ctx->q_Rotation);
    GetEulerFromQuaternion(ctx->q_Rotation, ctx->Rot);
	ctx->TransXYZ[0] = ctx->TransXYZ[1] = 0.;
}



void GetBox(magnoom_ctx *ctx)
{
	const float (*abc)[3] = ctx->abc;
	const int *uABC = ctx->uABC;
	float (*box)[3] = ctx->Box;
	//origine of the box is (0,0,0)
	//three vectors, b[1]+b[2]+b[2] define main diagonal of the box:
	box[0][0] = (uABC[0])*abc[0][0];
	box[0][1] = (uABC[0])*abc[0][1];
	box[0][2] = (uABC[0])*abc[0][2];

	box[1][0] = (uABC[1])*abc[1][0];
	box[1][1] = (uABC[1])*abc[1][1];
	box[1][2] = (uABC[1])*abc[1][2];

	box[2][0] = (uABC[2])*abc[2][0];
	box[2][1] = (uABC[2])*abc[2][1];
	box[2][2] = (uABC[2])*abc[2][2];
}

void ChangeBoxSize(magnoom_ctx *ctx, int Na, int Nb, int Nc){
	if (ctx->Play==1) ctx->Play=0;
}


void Display (magnoom_ctx *ctx)
{
	float mat[4*4]; // rotation matrix
    float   axisX[] = { 1, 0, 0 };
    float   axisY[] = { 0, 1, 0 };
    float   axisZ[] = { 0, 0, 1 };
    float   quat1[4];//NSK
    float   quat2[4];//NSK 
    float   quat3[4];//NSK 

	GLdouble Hight;
	float Vtemp[3];

	glClearColor( ctx->BackgroundColors[ctx->WhichBackgroundColor][0], ctx->BackgroundColors[ctx->WhichBackgroundColor][1], ctx->BackgroundColors[ctx->WhichBackgroundColor][2], 0. );// setup the clear values
	glDrawBuffer( GL_BACK );// erase the background
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_ACCUM_BUFFER_BIT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	if( ctx->WhichProjection == PERSP ) {
		mat4x4 projection;
		mat4x4_perspective(projection, D2R*ctx->PerspSet[0], ctx->asp_rat, ctx->PerspSet[2], ctx->PerspSet[3]);
		glMultMatrixf(&projection[0][0]);
	}
	else{
		Vtemp[0] = ctx->CameraEye[0] - ctx->CameraC[0];
		Vtemp[1] = ctx->CameraEye[1] - ctx->CameraC[1];
		Vtemp[2] = ctx->CameraEye[2] - ctx->CameraC[2]+ctx->TransXYZ[2];
		Hight = Unitf(Vtemp,Vtemp); //distance to the point which camera looks at
		Hight = -Hight*tan(D2R*ctx->PerspSet[0]/2.f); // hight of the view frame at the plane perp to cam view line
		glOrtho( Hight*ctx->asp_rat, -Hight*ctx->asp_rat, Hight,  -Hight,  ctx->PerspSet[2], ctx->PerspSet[3] );
	//         (     left,          right,      bottom,   top,       near,        far     )
	}

	// set the eye position, look-at position, and up-vector:
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	mat4x4 view;
	mat4x4_look_at(view, ctx->CameraEye, ctx->CameraC, ctx->CameraUp);
	glMultMatrixf(&view[0][0]);

	glEnable(GL_LIGHT0);
	glEnable(GL_NORMALIZE);
	GLfloat light_ambient[4] = {0.0f, 0.0f, 0.0f, 1.0f};
	GLfloat light_color[4] = {
		ctx->LightMultiplier,
		ctx->LightMultiplier,
		ctx->LightMultiplier,
		1.0f
	};
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_color);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_color);

	GLfloat light_position[4];
	if (ctx->WhichLightingMode == LIGHT_ADAPTIVE) {
		light_position[0] = ctx->CameraEye[0];
		light_position[1] = ctx->CameraEye[1];
		light_position[2] = ctx->CameraEye[2];
		light_position[3] = 1.0f;
	} else {
		light_position[0] = ctx->LightDirection[0];
		light_position[1] = ctx->LightDirection[1];
		light_position[2] = ctx->LightDirection[2];
		light_position[3] = 0.0f;
	}
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, 1.0f);
	glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.0f);
	glLightf(GL_LIGHT0, GL_QUADRATIC_ATTENUATION, 0.0f);

	if (ctx->WhichLightingMode == LIGHT_OFF) glDisable(GL_LIGHTING);
	else                                     glEnable(GL_LIGHTING);

	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, ctx->specular);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, ctx->shininess);
	glPushMatrix();

	// translate the scene:
	ctx->TransXYZ[0]+=ctx->dTransXYZ[0];
	ctx->TransXYZ[1]+=ctx->dTransXYZ[1];
	ctx->TransXYZ[2]+=ctx->dTransXYZ[2];

	glTranslatef( (GLfloat)ctx->TransXYZ[0], (GLfloat)ctx->TransXYZ[1], -(GLfloat)ctx->TransXYZ[2] );
    //new value for the Euler angles due to mose rotation:
    GetEulerFromQuaternion(ctx->q_Rotation, ctx->Rot); 
    //new directions of Cartesian axes:
    RotateVectorByQuaternion(axisX, ctx->q_Rotation);
    RotateVectorByQuaternion(axisY, ctx->q_Rotation);
    RotateVectorByQuaternion(axisZ, ctx->q_Rotation);
    //adding the rotation about each axes (from keyboard), to the rotation by mouse:
    SetQuaternionFromAxisAngle(axisX, D2R*ctx->dRot[0], quat1);
    SetQuaternionFromAxisAngle(axisY, D2R*ctx->dRot[1], quat2);
    SetQuaternionFromAxisAngle(axisZ, D2R*ctx->dRot[2], quat3);
    //combining all rotations by means of quaternion multiplication:
    MultiplyQuaternions( quat2, ctx->q_Rotation, ctx->q_Rotation);
    MultiplyQuaternions( quat3, ctx->q_Rotation, ctx->q_Rotation);
    MultiplyQuaternions( quat1, ctx->q_Rotation, ctx->q_Rotation);
    //
	ConvertQuaternionToMatrix(ctx->q_Rotation, mat);
    //
    glMultMatrixf(mat);


	drawVBO(ctx); // Draw VBO for spins
	DrawVBOMeshIndexed(&ctx->BextDC_mesh, GL_TRIANGLES, ctx->BextDC_mesh.index_count, 0);
	
	// possibly draw the box and periodic boundary condition :
	if( ctx->BoxOn != 0 ) {	
		DrawVBOMeshIndexed(&ctx->box_mesh, GL_TRIANGLES, ctx->box_mesh.index_count, 0);
		if(ctx->Boundary[0]!=0) DrawVBOMeshIndexed(&ctx->pbc_mesh[0], GL_TRIANGLES, ctx->pbc_mesh[0].index_count, 0);
		if(ctx->Boundary[1]!=0) DrawVBOMeshIndexed(&ctx->pbc_mesh[1], GL_TRIANGLES, ctx->pbc_mesh[1].index_count, 0);
		if(ctx->Boundary[2]!=0) DrawVBOMeshIndexed(&ctx->pbc_mesh[2], GL_TRIANGLES, ctx->pbc_mesh[2].index_count, 0);
	}
	// possibly draw the axes:
	if( ctx->AxesOn != 0 ) DrawVBOMeshIndexed(&ctx->basis_mesh, GL_TRIANGLES, ctx->basis_mesh.index_count, 0);

	glPopMatrix();

	glDisable(GL_LIGHTING);

	// Draw tweak bars
	TwDraw();
}

void setupOpenGL (magnoom_ctx *ctx)
{
#ifdef _WIN32
	// Must run before any top-level window is created (glfwOpenWindow()
	// below), or Windows has already committed to treating this process as
	// DPI-unaware for the window's lifetime.
	EnableWindowsDpiAwareness();
#endif
	if (!glfwInit()) {
		fprintf(stderr, "GLFW initialization failed\n");
		exit(1);
	}
	glfwOpenWindowHint(GLFW_FSAA_SAMPLES, ctx->N_Multisample);
	if (!glfwOpenWindow(ctx->window_width, ctx->window_height, 8, 8, 8, 8, 24, 0, GLFW_WINDOW)) {
		fprintf(stderr, "Cannot open GLFW window\n");
		glfwTerminate();
		exit(1);
	}
	glfwSetWindowTitle(ctx->WINDOWTITLE);
	glfwEnable(GLFW_KEY_REPEAT);
	glfwDisable(GLFW_AUTO_POLL_EVENTS);
	g_ctx = ctx;
	// Not gladLoadGLLoader(glfwGetProcAddress): GLFW2's own glfwGetProcAddress
	// (unlike GLFW3's) is a thin wglGetProcAddress() wrapper with no fallback
	// to GetProcAddress() on opengl32.dll, so it fails to resolve OpenGL 1.1
	// core functions on Windows (wglGetProcAddress is documented by Microsoft
	// as unreliable for those). gladLoadGL() is self-contained and already
	// has that fallback built in (see get_proc() in vendor/glad/src/glad.c) -
	// same fix already applied and confirmed working in the standalone
	// AntTweakBar-Legacy project's own GLFW examples.
	if (!gladLoadGL()) {
		fprintf(stderr, "Failed to initialize GLAD\n");
		glfwTerminate();
		exit(1);
	}

	for (int i=0;i<CAMERA_POSITION_SLOTS;i++){
		ctx->CameraPosition[i][0]=0;
		ctx->CameraPosition[i][1]=0;
		ctx->CameraPosition[i][2]=0;
		ctx->CameraPosition[i][3]=0;
		ctx->CameraPosition[i][4]=0;
		ctx->CameraPosition[i][5]=0;
		ctx->CameraPosition[i][6]=ctx->PerspSet[0];
	}

	glfwSetKeyCallback(GLFWKeyCallback);
	glfwSetCharCallback(GLFWCharCallback);
	glfwSetMouseButtonCallback(GLFWMouseButtonCallback);
	glfwSetMousePosCallback(GLFWCursorPosCallback);
	glfwSetMouseWheelCallback(GLFWScrollCallback);
	glfwSetWindowSizeCallback(GLFWWindowSizeCallback);

	// Best-effort early measurement, needed because AntTweakBar's
	// "fontscaling" below must be set before TwInit() and so cannot wait
	// for GLFWWindowSizeCallback's own re-measurement below.
	MeasureContentScale(ctx);

	// Initialize AntTweakBar
	char fontScalingDef[64];
	snprintf(fontScalingDef, sizeof(fontScalingDef), "GLOBAL fontscaling=%g", ctx->ContentScaleX);
	TwDefine(fontScalingDef);
	if (!TwInit(TW_OPENGL, NULL)) {
		fprintf(stderr, "TwInit failed: %s\n", TwGetLastError());
		glfwTerminate();
		exit(1);
	}
	{
		int windowWidth = 0, windowHeight = 0;
		glfwGetWindowSize(&windowWidth, &windowHeight);
		GLFWWindowSizeCallback(windowWidth, windowHeight);
	}
	// Startup settling is over: lock ContentScaleX/Y so every later, real
	// GLFWWindowSizeCallback (the user resizing the window) reuses this
	// fixed ratio instead of re-measuring it against glGetIntegerv
	// (GL_VIEWPORT) -- see MeasureContentScale()'s comment for why
	// re-measuring on a real resize is self-defeating.
	contentScaleLocked = 1;
	setupTweakBar(ctx);

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	GLfloat global_ambient[4] = {0.15f, 0.15f, 0.15f, 1.0f};
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);
	glShadeModel(GL_SMOOTH); //Enable smooth shading

	glCullFace(GL_FRONT);//GL_FRONT//GL_FRONT_AND_BACK
	glEnable(GL_CULL_FACE);

	glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);

	glEnable(GL_LINE_SMOOTH);               
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

	glEnable(GL_POLYGON_SMOOTH);
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_MULTISAMPLE);
}

void idle (magnoom_ctx *ctx)
{  
	ctx->currentTime = (int)(glfwGetTime() * 1000.0);
	ctx->timeInterval = ctx->currentTime - ctx->previousTime;

	if((ctx->timeInterval > 40 && ctx->Play!=0)||ctx->SpecialEvent!=0)//40ms gives approximately 25 FPS +/-1 if the engine works faster then 25 IPS
	{
		if( ctx->DataTransferState==TAKE_DATA || ctx->SpecialEvent!=0)
		{
			switch (ctx->SpecialEvent)
			{
				case 0:
				ChangeVectorMode(ctx, 1);
				case 1:
				ctx->totalEnergy = GetTotalEnergy(ctx);
				// totalEnergy = 0;
				ctx->mtot[0] = ctx->Mtot[0]/ctx->NOS;
				ctx->mtot[1] = ctx->Mtot[1]/ctx->NOS;
				ctx->mtot[2] = ctx->Mtot[2]/ctx->NOS;
                // int  k=123;
                // mtot[2] = sqrt(bSx[k]*bSx[k]+bSy[k]*bSy[k]+bSz[k]*bSz[k]);//metka
				ctx->perSpEnergy = ctx->totalEnergy/ctx->NOS;
				ctx->totalEnergyFerro = GetTotalEnergyFerro(ctx);
				// totalEnergyFerro = 0;
				ctx->totalEnergyFerro = ctx->totalEnergyFerro/ctx->NOS;
				ctx->perSpEnergyMinusFerro = ctx->perSpEnergy - ctx->totalEnergyFerro;
				ctx->SpecialEvent=0;
			}
			if (ctx->timeInterval > 500){//~0.5 second
				ctx->FPS = ctx->frameCount / (ctx->timeInterval * 0.002f);
				ctx->previousTime = ctx->currentTime;
				ctx->frameCount = 0;	
				if (ctx->Play!=0){
					ctx->IPS = (ctx->currentIteration - ctx->previousIteration)/ (ctx->timeInterval * 0.002f);;
					ctx->previousIteration = ctx->currentIteration;
				}
			}
			pthread_mutex_lock( &ctx->show_mutex);
				ctx->DataTransferState = WAIT_DATA; // meaning that OpenGL is waiting for new data from engine
			pthread_mutex_unlock( &ctx->show_mutex);
		} 

		ctx->frameCount++;
	}
}


void TW_CALL CB_SetN_Multisample(const void *value, void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	(void)ctx; // unused: this callback is disabled (see the commented-out registration below)
    // N_Multisample = 16;//*( int *)value;
    // N_Multisample = *( int *)value;
	// ChangeVectorMode (&mag_ctx, 0);
}


// void TW_CALL CB_GetN_Multisample(void *value, void *clientData)
// {
//     magnoom_ctx *ctx = (magnoom_ctx *)clientData;
//     *(int *)value = ctx->N_Multisample;
// }


void TW_CALL CB_Set_Run( const void *value, void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	ctx->Play = *( int *)value;
    if (ctx->Play!=0){
        pthread_mutex_lock(&ctx->culc_mutex);
			ctx->EngineIdle=false;
            ctx->EngineRunState=DO_IT;
            ctx->SleepTime=100;
        pthread_mutex_unlock(&ctx->culc_mutex);
    }else{
		magnoom_stop_engine(ctx);
    }
    // printf("IchBinHier!\n");

}

void TW_CALL CB_Get_Run(void *value, void *clientData)
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	*(float *)value = ctx->Play;
}


void TW_CALL CB_SetRotX(const void *value, void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    GetEulerFromQuaternion(ctx->q_Rotation, ctx->Rot); 
    ctx->Rot[0] = *( float *)value;
    GetQuaternionFromEuler(ctx->q_Rotation, ctx->Rot);
}

void TW_CALL CB_GetRotX(void *value, void *clientData)
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(float *)value = ctx->Rot[0];
}


void TW_CALL CB_SetRotY(const void *value, void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    GetEulerFromQuaternion(ctx->q_Rotation, ctx->Rot); 
    ctx->Rot[1] = *( float *)value;
    GetQuaternionFromEuler(ctx->q_Rotation, ctx->Rot);
}

void TW_CALL CB_GetRotY(void *value, void *clientData)
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(float *)value = ctx->Rot[1];
}

void TW_CALL CB_SetRotZ(const void *value, void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    GetEulerFromQuaternion(ctx->q_Rotation, ctx->Rot); 
    ctx->Rot[2] = *( float *)value;
    GetQuaternionFromEuler(ctx->q_Rotation, ctx->Rot);
}

void TW_CALL CB_GetRotZ(void *value, void *clientData)
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(float *)value = ctx->Rot[2];
}

void TW_CALL CB_SetScale(const void *value, void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    ctx->Scale = *( float *)value; // copy value to Scale
	ChangeVectorMode (ctx, 0);
}

//  Callback function called by the tweak bar to get the 'AutoRotate' value
void TW_CALL CB_GetScale(void *value, void *clientData)
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(float *)value = ctx->Scale; // just copy Scale to value
}


void TW_CALL CB_SetScaleBextDC(const void *value, void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    ctx->Scale_BextDC = *( float *)value;
	UpdateVerticesNormalsColors_BextDC(ctx);
    UploadVBOMesh(&ctx->BextDC_mesh, ctx->vertices_BextDC, ctx->normals_BextDC, ctx->colors_BextDC,
                  ctx->indices_BextDC, VBO_UPLOAD_ALL);
}

//  Callback function called by the tweak bar to get the 'AutoRotate' value
void TW_CALL CB_GetScaleBextDC(void *value, void *clientData)
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(float *)value = ctx->Scale_BextDC; // just copy Scale to value
}



void TW_CALL CB_SetVectorMode(const void *value, void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    ctx->WhichVectorMode = *( enVectorMode *)value; // copy value to WhichVectorMode
    ChangeVectorMode (ctx, 0);
}


void TW_CALL CB_GetVectorMode(void *value, void *clientData)
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(int *)value = ctx->WhichVectorMode; // just copy WhichVectorMode to value
}


void TW_CALL CB_SetFaces(const void *value, void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    ctx->arrowFaces = *( int *)value; // copy value to arrowFaces
	ChangeVectorMode (ctx, 0);
}


void TW_CALL CB_GetFaces(void *value, void *clientData)
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(int *)value = ctx->arrowFaces; // just copy arrowFaces to value
}

void TW_CALL CB_SetPivot(const void *value, void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    ctx->Pivot = *( float *)value; // copy value to Pivot
	ChangeVectorMode (ctx, 0);
}

void TW_CALL CB_GetPivot(void *value, void *clientData)
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(float *)value = ctx->Pivot; // just copy Pivot to value
}

void TW_CALL CB_SaveCameraPosition ( void *clientData )//metka has to be fixed beacuse of new quaternion rotation
{
	magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	ctx->CameraPosition[ctx->CurrentCameraPositionBank][0]=ctx->Rot[0];
	ctx->CameraPosition[ctx->CurrentCameraPositionBank][1]=ctx->Rot[1];
	ctx->CameraPosition[ctx->CurrentCameraPositionBank][2]=ctx->Rot[2];
	ctx->CameraPosition[ctx->CurrentCameraPositionBank][3]=ctx->TransXYZ[0];
	ctx->CameraPosition[ctx->CurrentCameraPositionBank][4]=ctx->TransXYZ[1];
	ctx->CameraPosition[ctx->CurrentCameraPositionBank][5]=ctx->TransXYZ[2];
	ctx->CameraPosition[ctx->CurrentCameraPositionBank][6]=ctx->PerspSet[0];
}

void TW_CALL CB_GetCameraPosition ( void *clientData )//metka has to be fixed beacuse of new quaternion rotation
{
	magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	ctx->Rot[0]=ctx->CameraPosition[ctx->CurrentCameraPositionBank][0];
	ctx->Rot[1]=ctx->CameraPosition[ctx->CurrentCameraPositionBank][1];
	ctx->Rot[2]=ctx->CameraPosition[ctx->CurrentCameraPositionBank][2];
	ctx->TransXYZ[0]=ctx->CameraPosition[ctx->CurrentCameraPositionBank][3];
	ctx->TransXYZ[1]=ctx->CameraPosition[ctx->CurrentCameraPositionBank][4];
	ctx->TransXYZ[2]=ctx->CameraPosition[ctx->CurrentCameraPositionBank][5];
	ctx->PerspSet[0]=ctx->CameraPosition[ctx->CurrentCameraPositionBank][6];
}

void TW_CALL CB_SetColorShift(const void *value, void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    ctx->ColorShift = *( int *)value; // copy value to ColorShift
		// if(WhichColorScheme==RGB)
		// {
		// HueMap[0]=HueMapRGB[0]+ColorShift;;
		// HueMap[1]=HueMapRGB[1]+ColorShift;;
		// HueMap[2]=HueMapRGB[2]+ColorShift;;
		// HueMap[3]=HueMapRGB[3]+ColorShift;;
		// HueMap[4]=HueMapRGB[4]+ColorShift;;
		// HueMap[5]=HueMapRGB[5]+ColorShift;;
		// } 	
		// else
		// {
		// HueMap[0]=HueMapRYGB[0]+ColorShift;;
		// HueMap[1]=HueMapRYGB[1]+ColorShift;;
		// HueMap[2]=HueMapRYGB[2]+ColorShift;;
		// HueMap[3]=HueMapRYGB[3]+ColorShift;;
		// HueMap[4]=HueMapRYGB[4]+ColorShift;;
		// HueMap[5]=HueMapRYGB[5]+ColorShift;;
		// }
		// InitRGB(RHue, GHue, BHue, HueMap);
		ChangeVectorMode (ctx, 1);
}

void TW_CALL CB_GetColorShift(void *value, void *clientData)
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(int *)value = ctx->ColorShift; // just copy ColorShift to value
}

void TW_CALL CB_SetInvHue(const void *value, void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    ctx->InvertHue = *( int *)value; // copy value to InvertHue
	ChangeVectorMode (ctx, 1);
}

void TW_CALL CB_GetInvHue(void *value, void *clientData)
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(int *)value = ctx->InvertHue; // just copy InvertHue to value
}

void TW_CALL CB_SetInvVal(const void *value, void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    ctx->InvertValue = *( int *)value; // copy value to InvertValue
	ChangeVectorMode (ctx, 1);
}

void TW_CALL CB_GetInvVal(void *value, void *clientData)
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(int *)value = ctx->InvertValue; // just copy InvertValue to value
}


void TW_CALL CB_SetBextDCMagnitude(const void *value, void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	ctx->BextDCMagnitude = *( float *)value;
    // metka
    // if (BextDCMagnitude>1.1*fabs(ctx->Jij[0])) BextDCMagnitude=1.1*fabs(ctx->Jij[0]);
    ctx->BextDC[0] = ctx->BextDCMagnitude * ctx->BextDCDirection[0];
    ctx->BextDC[1] = ctx->BextDCMagnitude * ctx->BextDCDirection[1];
    ctx->BextDC[2] = ctx->BextDCMagnitude * ctx->BextDCDirection[2];
	UpdateVerticesNormalsColors_BextDC(ctx);
	UploadVBOMesh(&ctx->BextDC_mesh, ctx->vertices_BextDC, ctx->normals_BextDC, ctx->colors_BextDC, ctx->indices_BextDC, VBO_UPLOAD_ALL);
}

void TW_CALL CB_GetBextDCMagnitude(void *value, void *clientData)
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	*(float *)value = ctx->BextDCMagnitude;
}


void TW_CALL CB_SetBextDCTheta(const void *value, void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    ctx->BextDCTheta = *(float*)value;
    ctx->BextDCDirection[0]=sin(PI*ctx->BextDCTheta/180)*cos(PI*ctx->BextDCPhi/180);
	ctx->BextDCDirection[1]=sin(PI*ctx->BextDCTheta/180)*sin(PI*ctx->BextDCPhi/180);
	ctx->BextDCDirection[2]=cos(PI*ctx->BextDCTheta/180);
	ctx->BextDC[0]=ctx->BextDCMagnitude*ctx->BextDCDirection[0];
	ctx->BextDC[1]=ctx->BextDCMagnitude*ctx->BextDCDirection[1];
	ctx->BextDC[2]=ctx->BextDCMagnitude*ctx->BextDCDirection[2];
	UpdateVerticesNormalsColors_BextDC(ctx);
	UploadVBOMesh(&ctx->BextDC_mesh, ctx->vertices_BextDC, ctx->normals_BextDC, ctx->colors_BextDC, ctx->indices_BextDC, VBO_UPLOAD_ALL);
}

void TW_CALL CB_GetBextDCTheta(void *value, void *clientData)
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(float*)value = ctx->BextDCTheta;
}

void TW_CALL CB_SetBextDCPhi(const void *value, void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    ctx->BextDCPhi = *(float*)value;
	ctx->BextDCDirection[0]=sin(PI*ctx->BextDCTheta/180)*cos(PI*ctx->BextDCPhi/180);
	ctx->BextDCDirection[1]=sin(PI*ctx->BextDCTheta/180)*sin(PI*ctx->BextDCPhi/180);
	ctx->BextDCDirection[2]=cos(PI*ctx->BextDCTheta/180);
	ctx->BextDC[0]=ctx->BextDCMagnitude*ctx->BextDCDirection[0];
	ctx->BextDC[1]=ctx->BextDCMagnitude*ctx->BextDCDirection[1];
	ctx->BextDC[2]=ctx->BextDCMagnitude*ctx->BextDCDirection[2];
	UpdateVerticesNormalsColors_BextDC(ctx);
	UploadVBOMesh(&ctx->BextDC_mesh, ctx->vertices_BextDC, ctx->normals_BextDC, ctx->colors_BextDC, ctx->indices_BextDC, VBO_UPLOAD_ALL);
}

void TW_CALL CB_GetBextDCPhi(void *value, void *clientData)
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(float*)value = ctx->BextDCPhi;
}

static const char *anisotropy_k2_control_names[ANISOTROPY_K2_COMPONENT_COUNT] = {
	"K11", "K12", "K13", "K22", "K23", "K33"
};

static const char *anisotropy_k4_control_names[ANISOTROPY_K4_COMPONENT_COUNT] = {
	"K1111", "K1112", "K1113", "K1122", "K1123", "K1133", "K1222", "K1223",
	"K1233", "K1333", "K2222", "K2223", "K2233", "K2333", "K3333"
};

static void anisotropy_refresh_controls(magnoom_ctx *ctx)
{
	if (ctx == NULL || ctx->anisotropy_bar == NULL) return;
	int visible = ctx->anisotropy_mode == ANISOTROPY_INDIVIDUAL;
	TwSetParam(ctx->anisotropy_bar, "Atom", "visible", TW_PARAM_INT32, 1, &visible);
	TwRefreshBar(ctx->anisotropy_bar);
}

void TW_CALL CB_SetAnisotropyMode(const void *value, void *clientData)
{
	magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	int mode = *(const int *)value;
	if (mode != ANISOTROPY_GLOBAL && mode != ANISOTROPY_INDIVIDUAL) return;
	ctx->anisotropy_mode = (AnisotropyMode)mode;
	anisotropy_refresh_controls(ctx);
}

void TW_CALL CB_GetAnisotropyMode(void *value, void *clientData)
{
	*(int *)value = ((magnoom_ctx *)clientData)->anisotropy_mode;
}

void TW_CALL CB_SetAnisotropyAtom(const void *value, void *clientData)
{
	magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	int atom = *(const int *)value;
	if (atom < 0 || atom >= ctx->AtomsPerBlock) return;
	ctx->anisotropy_selected_atom = atom;
	anisotropy_refresh_controls(ctx);
}

void TW_CALL CB_GetAnisotropyAtom(void *value, void *clientData)
{
	*(int *)value = ((magnoom_ctx *)clientData)->anisotropy_selected_atom;
}

void TW_CALL CB_SetAnisotropyQuaternion(const void *value, void *clientData)
{
	magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	int atom = anisotropy_site_index(ctx, ctx->anisotropy_selected_atom);
	(void)anisotropy_set_quaternion(ctx, atom, (const double *)value);
	anisotropy_refresh_controls(ctx);
}

void TW_CALL CB_GetAnisotropyQuaternion(void *value, void *clientData)
{
	magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	int atom = anisotropy_site_index(ctx, ctx->anisotropy_selected_atom);
	memcpy(value, ctx->anisotropy_quaternion[atom], sizeof(ctx->anisotropy_quaternion[atom]));
}

void TW_CALL CB_AnisotropyResetQuaternion(void *clientData)
{
	magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	int atom = anisotropy_site_index(ctx, ctx->anisotropy_selected_atom);
	(void)anisotropy_reset_quaternion(ctx, atom);
	anisotropy_refresh_controls(ctx);
}

void TW_CALL CB_AnisotropyApplyAxisAngle(void *clientData)
{
	magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	int atom = anisotropy_site_index(ctx, ctx->anisotropy_selected_atom);
	double axis[3] = {
		(double)ctx->anisotropy_rotation_axis[0],
		(double)ctx->anisotropy_rotation_axis[1],
		(double)ctx->anisotropy_rotation_axis[2]
	};
	(void)anisotropy_compose_axis_angle(ctx, atom, axis, (double)ctx->anisotropy_rotation_angle);
	anisotropy_refresh_controls(ctx);
}

void TW_CALL CB_SetAnisotropyComponent(const void *value, void *clientData)
{
	AnisotropyComponentControl *control = (AnisotropyComponentControl *)clientData;
	magnoom_ctx *ctx = control->ctx;
	double component_value = *(const float *)value;
	int atom = anisotropy_site_index(ctx, ctx->anisotropy_selected_atom);
	bool changed;
	if (control->kind == ANISOTROPY_COMPONENT_K2) {
		const int *index = anisotropy_k2_components[control->component];
		changed = anisotropy_set_k2(ctx, atom, index[0], index[1], component_value);
	} else {
		const int *index = anisotropy_k4_components[control->component];
		changed = anisotropy_set_k4(ctx, atom, index[0], index[1], index[2], index[3],
			component_value);
	}
	if (changed) anisotropy_rotate_site(ctx, atom);
}

void TW_CALL CB_GetAnisotropyComponent(void *value, void *clientData)
{
	AnisotropyComponentControl *control = (AnisotropyComponentControl *)clientData;
	magnoom_ctx *ctx = control->ctx;
	int atom = anisotropy_site_index(ctx, ctx->anisotropy_selected_atom);
	const AnisotropyTensor *tensor = &ctx->anisotropy_local[atom];
	if (control->kind == ANISOTROPY_COMPONENT_K2) {
		const int *index = anisotropy_k2_components[control->component];
		*(float *)value = (float)tensor->K2[index[0]][index[1]];
	} else {
		const int *index = anisotropy_k4_components[control->component];
		*(float *)value = (float)tensor->K4[index[0]][index[1]][index[2]][index[3]];
	}
}

void TW_CALL CB_CopyAnisotropyAtom0(void *clientData)
{
	magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	if (!anisotropy_copy_atom0_tensors(ctx)) return;
	anisotropy_refresh_controls(ctx);
}

void TW_CALL CB_ExportAnisotropyMap(void *clientData)
{
	magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	if (!anisotropy_export_energy_maps(ctx)) {
		fprintf(stderr, "Anisotropy map export failed.\n");
	}
}



// void TW_CALL CB_SetBextDCDirection(const void *value, void *clientData )
// {
// 	(void)clientData; // unused
//  //    BextDCTheta = *(float*)value;
//  //    BextDCDirection[0]=sin(PI*BextDCTheta/180)*cos(PI*BextDCPhi/180);
// 	// BextDCDirection[1]=sin(PI*BextDCTheta/180)*sin(PI*BextDCPhi/180);
// 	// BextDCDirection[2]=cos(PI*BextDCTheta/180);
// }

// void TW_CALL CB_GetBextDCDirection(void *value, void *clientData)
// {
//     (void)clientData; 
//     // *(float*)value = BextDCTheta; //metka check this comment!
// }


void TW_CALL CB_SetNumImages(const void *value, void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	ctx->num_images = *( int *)value;
    ReallocateMemoryForImages(ctx, ctx->num_images, ctx->NOS);
}

void TW_CALL CB_GetNumImages(void *value, void *clientData)
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	*(int *)value = ctx->num_images;
}

void TW_CALL CB_SetBextACPeriod(const void *value, void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    ctx->BextACPeriod = *( double *)value; // copy value to BextACPeriod
    ctx->BextACOmega = TPI/ctx->BextACPeriod;
}

void TW_CALL CB_GetBextACPeriod(void *value, void *clientData)
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(double *)value = ctx->BextACPeriod; // just copy BextACPeriod to value
}

void TW_CALL CB_SetBextACOmega(const void *value, void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	ctx->BextACOmega = *( double *)value;
    ctx->BextACPeriod = TPI/ctx->BextACOmega;
}

void TW_CALL CB_GetBextACOmega(void *value, void *clientData)
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(double *)value = ctx->BextACOmega; // just copy BextACOmega to value
}

void TW_CALL CB_SetInitial( void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	ChangeInitialState(ctx);
}

void TW_CALL CB_SetShape( void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	UpdateKind(ctx);
	ChangeVectorMode(ctx, 1);
}

void TW_CALL CB_RotateAllSpins( void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	if(fabs(ctx->chDir[0])+fabs(ctx->chDir[1])+fabs(ctx->chDir[2])!=0)
	{
        double tmp[3];
		for (int i=0; i<ctx->NOS; i++) {
			RotateVector(VEC_X(ctx->S,i), VEC_Y(ctx->S,i), VEC_Z(ctx->S,i), ctx->chDir[0], ctx->chDir[1], ctx->chDir[2], ctx->RotateAllSpins, tmp);
			VEC_X(ctx->bS,i) = VEC_X(ctx->S,i) = tmp[0];
			VEC_Y(ctx->bS,i) = VEC_Y(ctx->S,i) = tmp[1];
			VEC_Z(ctx->bS,i) = VEC_Z(ctx->S,i) = tmp[2];
		}
		ChangeVectorMode(ctx, 1);		
	}
}

void TW_CALL CB_InvertX( void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	for (int i=0; i<ctx->NOS; i++) {VEC_X(ctx->S,i) = -VEC_X(ctx->S,i); VEC_X(ctx->bS,i) = -VEC_X(ctx->bS,i);}
	ChangeVectorMode(ctx, 1);
}

void TW_CALL CB_InvertY( void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	for (int i=0; i<ctx->NOS; i++) {VEC_Y(ctx->S,i) = -VEC_Y(ctx->S,i); VEC_Y(ctx->bS,i) = -VEC_Y(ctx->bS,i);}
	ChangeVectorMode(ctx, 1);
}

void TW_CALL CB_InvertZ( void *clientData ){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	for (int i=0; i<ctx->NOS; i++) {VEC_Z(ctx->S,i) = -VEC_Z(ctx->S,i); VEC_Z(ctx->bS,i) = -VEC_Z(ctx->bS,i);}
	ChangeVectorMode(ctx, 1);
}

void TW_CALL CB_CleanSxSySzFile( void *clientData ){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	magnoom_reset_record_file(ctx);
}

void TW_CALL CB_SetInputDirectory(const void *value, void *clientData)
{
	magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	if (!magnoom_copy_path(ctx->input_directory,
		sizeof(ctx->input_directory), (const char *)value)) {
		fprintf(stderr, "Input directory is too long.\n");
	}
}

void TW_CALL CB_GetInputDirectory(void *value, void *clientData)
{
	magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	magnoom_copy_path((char *)value, sizeof(ctx->input_directory), ctx->input_directory);
}

void TW_CALL CB_SetOutputDirectory(const void *value, void *clientData)
{
	magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	magnoom_change_output_directory(ctx, (const char *)value);
}

void TW_CALL CB_GetOutputDirectory(void *value, void *clientData)
{
	magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	magnoom_copy_path((char *)value, sizeof(ctx->output_directory), ctx->output_directory);
}

void TW_CALL CB_SetSliceMode(const void *value, void *clientData ){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    ctx->WhichSliceMode = *( enSliceMode *)value; // copy value to ctx->WhichSliceMode
    ChangeVectorMode(ctx, 0);
}


void TW_CALL CB_GetSliceMode(void *value, void *clientData){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(int *)value = ctx->WhichSliceMode; // just copy ctx->WhichSliceMode to value
}

void TW_CALL CB_GetThetaMax1(void *value, void *clientData){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(int *)value = ctx->theta_max1; 
}

void TW_CALL CB_SetThetaMax1(const void *value, void *clientData ){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	int test= *( int *)value; 
	if (test>=ctx->theta_min1 ){
        ctx->theta_max1 = test;
        ctx->Sz_min1=cos(ctx->theta_max1*PI/180.0);
        ChangeVectorMode(ctx, 0);
	}
}

void TW_CALL CB_GetThetaMax2(void *value, void *clientData){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(int *)value = ctx->theta_max2; 
}

void TW_CALL CB_SetThetaMax2(const void *value, void *clientData ){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	int test= *( int *)value; 
	if (test>=ctx->theta_min2 ){
        ctx->theta_max2 = test;
        ctx->Sz_min2=cos(ctx->theta_max2*PI/180.0);
        ChangeVectorMode(ctx, 0);
	}
}

void TW_CALL CB_GetThetaMax3(void *value, void *clientData){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(int *)value = ctx->theta_max3; 
}

void TW_CALL CB_SetThetaMax3(const void *value, void *clientData ){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	int test= *( int *)value; 
	if (test>=ctx->theta_min3 ){
        ctx->theta_max3 = test;
        ctx->Sz_min3=cos(ctx->theta_max3*PI/180.0);
        ChangeVectorMode(ctx, 0);
	}
}

void TW_CALL CB_GetThetaMin1(void *value, void *clientData){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(int *)value = ctx->theta_min1; 
}

void TW_CALL CB_SetThetaMin1(const void *value, void *clientData ){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	int test= *( int *)value; 
	if (test<=ctx->theta_max1 ){
        ctx->theta_min1 = test; 
        ctx->Sz_max1=cos(ctx->theta_min1*PI/180.0);
        ChangeVectorMode(ctx, 0);
	}
}

void TW_CALL CB_GetThetaMin2(void *value, void *clientData){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(int *)value = ctx->theta_min2; 
}

void TW_CALL CB_SetThetaMin2(const void *value, void *clientData ){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	int test= *( int *)value; 
	if (test<=ctx->theta_max2 ){
        ctx->theta_min2 = test; 
        ctx->Sz_max2=cos(ctx->theta_min2*PI/180.0);
        ChangeVectorMode(ctx, 0);
	}
}

void TW_CALL CB_GetThetaMin3(void *value, void *clientData){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(int *)value = ctx->theta_min3; 
}

void TW_CALL CB_SetThetaMin3(const void *value, void *clientData ){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	int test= *( int *)value; 
	if (test<=ctx->theta_max3 ){
        ctx->theta_min3 = test; 
        ctx->Sz_max3=cos(ctx->theta_min3*PI/180.0);
        ChangeVectorMode(ctx, 0);
	}
}

void TW_CALL CB_GetPhiMax1(void *value, void *clientData){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(int *)value = ctx->phi_max1; 
}


void TW_CALL CB_SetPhiMax1(const void *value, void *clientData ){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	int test= *( int *)value; 
	if (test>=ctx->phi_min1 ){
        ctx->phi_max1 = test;
        ChangeVectorMode(ctx, 0);
	}
}

void TW_CALL CB_GetPhiMin1(void *value, void *clientData){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(int *)value = ctx->phi_min1; 
}


void TW_CALL CB_SetPhiMin1(const void *value, void *clientData ){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	int test= *( int *)value; 
	if (test<=ctx->phi_max1 ){
        ctx->phi_min1 = test;
        ChangeVectorMode(ctx, 0);
	}
}

void TW_CALL CB_GetPhiMax2(void *value, void *clientData){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(int *)value = ctx->phi_max2; 
}


void TW_CALL CB_SetPhiMax2(const void *value, void *clientData ){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	int test= *( int *)value; 
	if (test>=ctx->phi_min2 ){
        ctx->phi_max2 = test;
        ChangeVectorMode(ctx, 0);
	}
}

void TW_CALL CB_GetPhiMin2(void *value, void *clientData){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(int *)value = ctx->phi_min2; 
}


void TW_CALL CB_SetPhiMin2(const void *value, void *clientData ){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	int test= *( int *)value; 
	if (test<=ctx->phi_max2 ){
        ctx->phi_min2 = test;
        ChangeVectorMode(ctx, 0);
	}
}

void TW_CALL CB_GetPhiMax3(void *value, void *clientData){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(int *)value = ctx->phi_max3; 
}


void TW_CALL CB_SetPhiMax3(const void *value, void *clientData ){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	int test= *( int *)value; 
	if (test>=ctx->phi_min3 ){
        ctx->phi_max3 = test;
        ChangeVectorMode(ctx, 0);
	}
}

void TW_CALL CB_GetPhiMin3(void *value, void *clientData){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(int *)value = ctx->phi_min3; 
}

void TW_CALL CB_SetPhiMin3(const void *value, void *clientData ){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	int test= *( int *)value; 
	if (test<=ctx->phi_max3 ){
        ctx->phi_min3 = test;
        ChangeVectorMode(ctx, 0);
	}
}
/**********************************************************************/
void TW_CALL CB_GetGreedMaxA(void *value, void *clientData){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(int *)value = ctx->GreedFilterMaxA; 
}

void TW_CALL CB_SetGreedMaxA(const void *value, void *clientData ){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	int test= *( int *)value; 
	if (test<ctx->uABC[0] && test>ctx->GreedFilterMinA){
        ctx->GreedFilterMaxA = test;
        ChangeVectorMode(ctx, 0);
	}
}

void TW_CALL CB_GetGreedMinA(void *value, void *clientData){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(int *)value = ctx->GreedFilterMinA; 
}

void TW_CALL CB_SetGreedMinA(const void *value, void *clientData ){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	int test= *( int *)value; 
	if (test<ctx->GreedFilterMaxA && test>=0){
        ctx->GreedFilterMinA = test;
        ChangeVectorMode(ctx, 0);
	}
}
/*************/
void TW_CALL CB_GetGreedMaxB(void *value, void *clientData){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(int *)value = ctx->GreedFilterMaxB; 
}

void TW_CALL CB_SetGreedMaxB(const void *value, void *clientData ){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	int test= *( int *)value; 
	if (test<ctx->uABC[1] && test>ctx->GreedFilterMinB){
        ctx->GreedFilterMaxB = test;
        ChangeVectorMode(ctx, 0);
	}
}

void TW_CALL CB_GetGreedMinB(void *value, void *clientData){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(int *)value = ctx->GreedFilterMinB; 
}

void TW_CALL CB_SetGreedMinB(const void *value, void *clientData ){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	int test= *( int *)value; 
	if (test<ctx->GreedFilterMaxB && test>=0){
        ctx->GreedFilterMinB = test;
        ChangeVectorMode(ctx, 0);
	}
}
/*************/
void TW_CALL CB_GetGreedMaxC(void *value, void *clientData){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(int *)value = ctx->GreedFilterMaxC; 
}

void TW_CALL CB_SetGreedMaxC(const void *value, void *clientData ){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	int test= *( int *)value; 
	if (test<ctx->uABC[2] && test>ctx->GreedFilterMinC){
        ctx->GreedFilterMaxC = test;
        ChangeVectorMode(ctx, 0);
	}
}

void TW_CALL CB_GetGreedMinC(void *value, void *clientData){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(int *)value = ctx->GreedFilterMinC; 
}

void TW_CALL CB_SetGreedMinC(const void *value, void *clientData ){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	int test= *( int *)value; 
	if (test<ctx->GreedFilterMaxC && test>=0){
        ctx->GreedFilterMinC = test;
        ChangeVectorMode(ctx, 0);
	}
}

void TW_CALL CB_GetGreedFilterInvert(void *value, void *clientData){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(bool *)value = ctx->GreedFilterInvert; 
}

void TW_CALL CB_SetGreedFilterInvert(const void *value, void *clientData ){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	ctx->GreedFilterInvert= *( bool *)value; 
    ChangeVectorMode(ctx, 0);
}

void TW_CALL CB_GetSpinFilter1(void *value, void *clientData){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(bool *)value = ctx->SpinFilter1; 
}

void TW_CALL CB_SetSpinFilter1(const void *value, void *clientData ){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	ctx->SpinFilter1= *( bool *)value; 
    ChangeVectorMode(ctx, 0);
}

void TW_CALL CB_GetSpinFilter2(void *value, void *clientData){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(bool *)value = ctx->SpinFilter2; 
}

void TW_CALL CB_SetSpinFilter2(const void *value, void *clientData ){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	ctx->SpinFilter2= *( bool *)value; 
    ChangeVectorMode(ctx, 0);
}

void TW_CALL CB_GetSpinFilter3(void *value, void *clientData){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    *(bool *)value = ctx->SpinFilter3; 
}

void TW_CALL CB_SetSpinFilter3(const void *value, void *clientData ){
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	ctx->SpinFilter3= *( bool *)value; 
    ChangeVectorMode(ctx, 0);
}

void TW_CALL CB_GetGreedFilter(void *value, void *clientData)
{
	magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	*(bool *)value = ctx->GreedFilter; 
}

void TW_CALL CB_SetGreedFilter(const void *value, void *clientData )
{
	magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	ctx->GreedFilter= *( bool *)value; 
	ChangeVectorMode(ctx, 0);
}


void TW_CALL CB_ResetIterations( void *clientData )
{
	magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	ctx->ITERATION=0;
	ctx->currentIteration=0;
	ctx->SpecialEvent=1;
}



void TW_CALL CB_SaveCSV( void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    char output_path[MAGNOOM_PATH_CAPACITY];
    FILE * pFile;
    if (!magnoom_resolve_output_path(ctx, ctx->outputfilename,
        output_path, sizeof(output_path))) return;
    pFile = fopen(output_path,"w");
   
    if (pFile!=NULL)
    { 	
        // fputs ("px,py,pz,nx,ny,nz,\n",pFile);
        // for (int i=0;i<ctx->NOS;i++)
        // {
        //     snprintf(shortBufer,200,"%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,\n",ctx->PosX[i],ctx->PosY[i],ctx->PosZ[i],bSx[i],bSy[i],bSz[i]);
        //     fputs (shortBufer,pFile);
        // }
        // fclose (pFile);

    	int an,bn,cn,atom,n,N;
    	int anini=0;
    	int anfin=ctx->uABC[0];
    	int bnini=0;
    	int bnfin=ctx->uABC[1];
    	int cnini=0;
    	int cnfin=ctx->uABC[2];
    	if (ctx->save_slice==1){
    		switch( ctx->WhichSliceMode){
    			case A_AXIS:
    				anini=ctx->A_layer_min-1;
    		        anfin=ctx->A_layer_max;
    			break;
    			case B_AXIS:
    				bnini=ctx->B_layer_min-1;
    		        bnfin=ctx->B_layer_max;
    			break;
    			case C_AXIS:
    				cnini=ctx->C_layer_min-1;
    		        cnfin=ctx->C_layer_max;
    			break;
    			default:
    			break;
    		}
    	}

    	if (ctx->WhichAverageMode==ALONG_0){
    		printf("Write to csv-file!\n");

            // for export from x,y,z--> x,z,Nc-y
            // from Filipp's code
            // for (bn = bnini; bn<bnfin; bn++) {
            //     for (cn = cnfin-1; cn>=cnini; cn--) {
            //         for (an = anini; an<anfin; an++) {
            //             n = an+bn*ctx->uABC[0]+cn*ctx->uABC[0]*ctx->uABC[1];// index of the block
            //             n = n*AtomsPerBlock;//index of the first spin in the block
            //             for (atom=0; atom<AtomsPerBlock; atom++){
            //                 N = n + atom;
            //                 snprintf(shortBufer,200,"%.14g,%.14g,%.14g,\n", ctx->Sx[N],-ctx->Sz[N],ctx->Sy[N]);
            //                 fputs (shortBufer,pFile); 
            //             }
            //         }
            //     }
            // }
            // fputs ("px,py,pz,nx,ny,nz,\n",pFile);//metka
            // fputs ("px,py,pz,nx,ny,nz,T,F\n",pFile);
            for (cn = cnini; cn<cnfin; cn++) {
    			for (bn = bnini; bn<bnfin; bn++) {
    				// for (an = anini; an<anfin; an++) {
    				for (an = anfin-1; an>=anini; an--) {
    					n = an+bn*ctx->uABC[0]+cn*ctx->uABC[0]*ctx->uABC[1];// index of the block
    					n = n*ctx->AtomsPerBlock;//index of the first spin in the block
    					for (atom=0; atom<ctx->AtomsPerBlock; atom++){
    					    N = n + atom;
                            //metka
    						snprintf(ctx->ShortBuffer,200,"%2.5f,%2.5f,%2.5f,%0.15f,%0.15f,%0.15f,\n",ctx->PosX[N],ctx->PosY[N],ctx->PosZ[N],VEC_X(ctx->bS,N),VEC_Y(ctx->bS,N),VEC_Z(ctx->bS,N));
                            // snprintf(shortBufer,200,"%.14g,%.14g,%.14g,%.14g,%.14g,%.14g,%.14g,%.14g,\n",
                                // ctx->PosX[N],ctx->PosY[N],ctx->PosZ[N],bSx[N],bSy[N],bSz[N],acos(bSz[N])*R2D,atan2 (bSy[N],bSx[N])*R2D);
    						fputs (ctx->ShortBuffer,pFile); 
    					}
    				}
    			}
    		}
		}else{// WhichAverageMode==ALONG_A or ALONG_B or ALONG_C
    		double tSpin[3];
    		float tPositin[3]; 
    		fputs ("px,py,pz,<nx>,<ny>,<nz>,|<n>|\n",pFile);
    		if (ctx->WhichAverageMode==ALONG_C){
    			for (bn = bnini; bn<bnfin; bn++){
    				for (an = anini; an<anfin; an++){
    					tSpin[0]=0;
    					tSpin[1]=0;
    					tSpin[2]=0;
    					tPositin[0]=ctx->abc[0][0]*0.5+ctx->abc[1][0]*0.5+ctx->abc[2][0]*0.5;
    					tPositin[1]=ctx->abc[0][1]*0.5+ctx->abc[1][1]*0.5+ctx->abc[2][1]*0.5;
    					tPositin[2]=ctx->abc[0][2]*0.5+ctx->abc[1][2]*0.5+ctx->abc[2][2]*0.5;
    					for (cn = cnini; cn<cnfin; cn++){
    						n = an+bn*ctx->uABC[0]+cn*ctx->uABC[0]*ctx->uABC[1];// index of the block
    						n = n*ctx->AtomsPerBlock;//index of the first spin in the block
    						for (atom=0; atom<ctx->AtomsPerBlock; atom++){
    						    N = n + atom;
    						    tSpin[0]+=VEC_X(ctx->bS,N);
    						    tSpin[1]+=VEC_Y(ctx->bS,N);
    						    tSpin[2]+=VEC_Z(ctx->bS,N);
    						}
    					}
    					tPositin[0]+=ctx->abc[0][0]*an+ctx->abc[1][0]*bn;
    					tPositin[1]+=ctx->abc[0][1]*an+ctx->abc[1][1]*bn;
    					tPositin[2]+=ctx->abc[0][2]*an+ctx->abc[1][2]*bn;
    					double modulus=sqrt(tSpin[0]*tSpin[0] + tSpin[1]*tSpin[1] + tSpin[2]*tSpin[2]);
    					int cN = cnfin - cnini;
    					snprintf(ctx->ShortBuffer,200,"%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,\n",tPositin[0],tPositin[1],tPositin[2],tSpin[0]/cN,tSpin[1]/cN,tSpin[2]/cN,modulus/cN);
    					fputs (ctx->ShortBuffer,pFile);
    				}
    			}
    		}else if(ctx->WhichAverageMode==ALONG_B){
    			for (an = anini; an<anfin; an++){
    				for (cn = cnini; cn<cnfin; cn++){
    					tSpin[0]=0;
    					tSpin[1]=0;
    					tSpin[2]=0;
    					tPositin[0]=ctx->abc[0][0]*0.5+ctx->abc[1][0]*0.5+ctx->abc[2][0]*0.5;
    					tPositin[1]=ctx->abc[0][1]*0.5+ctx->abc[1][1]*0.5+ctx->abc[2][1]*0.5;
    					tPositin[2]=ctx->abc[0][2]*0.5+ctx->abc[1][2]*0.5+ctx->abc[2][2]*0.5;
    					for (bn = bnini; bn<bnfin; bn++){
    						n = an+bn*ctx->uABC[0]+cn*ctx->uABC[0]*ctx->uABC[1];// index of the block
    						n = n*ctx->AtomsPerBlock;//index of the first spin in the block
    						for (atom=0; atom<ctx->AtomsPerBlock; atom++){
    						    N = n + atom;
    						    tSpin[0]+=VEC_X(ctx->bS,N);
    						    tSpin[1]+=VEC_Y(ctx->bS,N);
    						    tSpin[2]+=VEC_Z(ctx->bS,N);
    						}
    					}
    					tPositin[0]+=ctx->abc[0][0]*an+ctx->abc[2][0]*cn;
    					tPositin[1]+=ctx->abc[0][1]*an+ctx->abc[2][1]*cn;
    					tPositin[2]+=ctx->abc[0][2]*an+ctx->abc[2][2]*cn;
    					double modulus=sqrt(tSpin[0]*tSpin[0] + tSpin[1]*tSpin[1] + tSpin[2]*tSpin[2]);
    					int bN = bnfin - bnini;
    					snprintf(ctx->ShortBuffer,200,"%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,\n",tPositin[0],tPositin[1],tPositin[2],tSpin[0]/bN,tSpin[1]/bN,tSpin[2]/bN,modulus/bN);
    					fputs (ctx->ShortBuffer,pFile);
    				}
    			}
    		}else if(ctx->WhichAverageMode==ALONG_A){
    			for (cn = cnini; cn<cnfin; cn++){
    				for (bn = bnini; bn<bnfin; bn++){
    					tSpin[0]=0;
    					tSpin[1]=0;
    					tSpin[2]=0;
    					tPositin[0]=ctx->abc[0][0]*0.5+ctx->abc[1][0]*0.5+ctx->abc[2][0]*0.5;
    					tPositin[1]=ctx->abc[0][1]*0.5+ctx->abc[1][1]*0.5+ctx->abc[2][1]*0.5;
    					tPositin[2]=ctx->abc[0][2]*0.5+ctx->abc[1][2]*0.5+ctx->abc[2][2]*0.5;
    					for (an = anini; an<anfin; an++){
    						n = an+bn*ctx->uABC[0]+cn*ctx->uABC[0]*ctx->uABC[1];// index of the block
    						n = n*ctx->AtomsPerBlock;//index of the first spin in the block
    						for (atom=0; atom<ctx->AtomsPerBlock; atom++){
    						    N = n + atom;
    						    tSpin[0]+=VEC_X(ctx->bS,N);
    						    tSpin[1]+=VEC_Y(ctx->bS,N);
    						    tSpin[2]+=VEC_Z(ctx->bS,N);
    						}
    					}
    					tPositin[0]+=ctx->abc[2][0]*cn+ctx->abc[1][0]*bn;
    					tPositin[1]+=ctx->abc[2][1]*cn+ctx->abc[1][1]*bn;
    					tPositin[2]+=ctx->abc[2][2]*cn+ctx->abc[1][2]*bn;
    					double modulus=sqrt(tSpin[0]*tSpin[0] + tSpin[1]*tSpin[1] + tSpin[2]*tSpin[2]);
    					int aN = anfin - anini;
    					snprintf(ctx->ShortBuffer,200,"%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,\n",tPositin[0],tPositin[1],tPositin[2],tSpin[0]/aN,tSpin[1]/aN,tSpin[2]/aN,modulus/aN);
    					fputs (ctx->ShortBuffer,pFile);
    				}
			}
		printf("averaged output --> %s done!\n", output_path);
		}
		}
        fclose (pFile);
	} else {
		fprintf(stderr, "Cannot open output file '%s': %s\n", output_path, strerror(errno));
    }
}

void TW_CALL CB_ReadCSV( void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
    char input_path[MAGNOOM_PATH_CAPACITY];
    if (!magnoom_resolve_input_path(ctx, input_path, sizeof(input_path))) return;
    FILE * pFile = fopen(input_path, "r");
    if(pFile) {
		char c;
		char line[120];
		float px,py,pz,sx,sy,sz;
		int pos=0;
        //unkoment if the csv file has headers of the columns 
        // do{ // read titles of the columns 
        //  	c = (char)fgetc(pFile);//get char and move pointer to the next position
        //  	if (c != EOF) {
        //  		 	line[pos++] = c;	
        //  	}
        //  }while(c != EOF && c != '\n');
		//printf("%s\n", line);
    	int i=0;//line (spin) index
		do { // read all lines in file
    		pos = 0;//initial position in line
		    do{ // read one line
		    	c = (char)fgetc(pFile);//get char and move pointer to the next position
		    	if (c != EOF) {
		    			if (c==',') c = ' ';
		    		 	line[pos++] = c;	
		    	}
		    }while(c != EOF && c != '\n');

		    line[pos] = 0;//add a line end.
			sscanf(line, "%f %f %f %f %f %f ", &px,&py,&pz,&sx,&sy,&sz);
            // sscanf(line, "%f %f %f ", &sx,&sy,&sz);
			//printf("%f,%f,%f,%f,%f,%f\n", px,py,pz,sx,sy,sz);
			if (i<ctx->NOS) 
			{
				// ctx->PosX[i]=px;
				// ctx->PosY[i]=py;
				// ctx->PosZ[i]=pz;
				VEC_X(ctx->bS,i)=VEC_X(ctx->S,i)=sx;
				VEC_Y(ctx->bS,i)=VEC_Y(ctx->S,i)=sy;
				VEC_Z(ctx->bS,i)=VEC_Z(ctx->S,i)=sz;
		 	}
		 	i++;
		}while(c != EOF); 
		fclose(pFile);
		magnoom_reset_solver_state(ctx);
	} else {
		fprintf(stderr, "Cannot open input file '%s': %s\n", input_path, strerror(errno));
	}
	ChangeVectorMode(ctx, 1);
}



void TW_CALL CB_ReadOVF( void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	char input_path[MAGNOOM_PATH_CAPACITY];
	char  line[256];//whole line of header should be not longer then 256 characters
	int   lineLength=0;
	int   valuedim=3;
	int   xnodes;
	int   ynodes;
	int   znodes;
    char  keyW1 [256];//key word 1
    char  keyW2 [256];//key word 2
    char  keyW3 [256];//key word 3
    int   binType = 4;
	float temp4_x, temp4_y, temp4_z;
	double temp8_x, temp8_y, temp8_z;
	if (!magnoom_resolve_input_path(ctx, input_path, sizeof(input_path))) return;
    FILE * FilePointer = fopen(input_path, "rb");
	if(FilePointer!=NULL) {	
		lineLength=ReadHeaderLine(FilePointer, line);//read and check the first nonempty line which starts with '#'
		if (lineLength==-1) {// if there are no one line which starts with '#'
			printf("%s has a wrong file format! \n", input_path);
		}else{
		    sscanf(line, "# %s %s %s", keyW1, keyW2, keyW3 );
		    if(strncmp(keyW1, "OOMMF",5)!=0 || strncmp(keyW2, "OVF",  3)!=0 || strncmp(keyW3, "2.0",  3)!=0){
		        //if the first line isn't "OOMMF OFV 2.0"
			printf("%s has wrong header of wrong file format! \n", input_path);
		    	lineLength=-1;
		    }
		}
		//READING HEADER
		if (lineLength!=-1){
			do{
				lineLength = ReadHeaderLine(FilePointer, line);
				sscanf(line, "# %s %s %s", keyW1, keyW2, keyW3 );
				//printf("%s %s %s\n", keyW1, keyW2, keyW3);
				if (strncmp(keyW1, "valuedim:",9)==0) {
					sscanf(keyW2, "%d", &valuedim );
					printf("valuedim=%d\n", valuedim);					
				}else if (strncmp(keyW1, "xnodes:",7)==0) {
					sscanf(keyW2, "%d", &xnodes );
					printf("xnodes=%d\n", xnodes);
				}else if (strncmp(keyW1, "ynodes:",7)==0) {
					sscanf(keyW2, "%d", &ynodes );
					printf("ynodes=%d\n", ynodes);					
				}else if (strncmp(keyW1, "znodes:",7)==0) {
					sscanf(keyW2, "%d", &znodes );
					printf("znodes=%d\n", znodes);					
				} 
			}while(!(strncmp(keyW1, "Begin:",6)==0 && strncmp(keyW2, "Data",4)==0) && lineLength != -1 );
		}
        //READING DATA
		if (valuedim!=0 && xnodes!=0 && ynodes!=0 && znodes!=0){
			sscanf(line, "#%*s %s %s %s", keyW1, keyW2, keyW3 );
			//int imax,jmax,kmax;
			//if (xnodes>ctx->uABC[0]) {imax = ctx->uABC[0];}else{imax = xnodes;}
			//if (ynodes>ctx->uABC[1]) {jmax = ctx->uABC[1];}else{jmax = ynodes;}
			//if (znodes>ctx->uABC[2]) {kmax = ctx->uABC[2];}else{kmax = znodes;}
			int n;
			if (strncmp(keyW2, "Text",4)==0){
				//Text data format
				printf("...reading data in text format: %s \n", input_path);
				for (int k=0; k<znodes; k++){
					for (int j=0; j<ynodes; j++){
						for (int i=0; i<xnodes; i++){
							ReadDataLine(FilePointer, line);
							if (k<ctx->uABC[2] && j<ctx->uABC[1] && i<ctx->uABC[0]){
								n = i + j*ctx->uABC[0] + k*ctx->uABC[0]*ctx->uABC[1];
								sscanf(line, "%lf %lf %lf", &VEC_X(ctx->bS,n),&VEC_Y(ctx->bS,n),&VEC_Z(ctx->bS,n));
								VEC_X(ctx->S,n)=VEC_X(ctx->bS,n); VEC_Y(ctx->S,n)=VEC_Y(ctx->bS,n);VEC_Z(ctx->S,n)=VEC_Z(ctx->bS,n);
							}
						}
					}
				}
			}else if (strncmp(keyW2, "Binary",6)==0){
				if(strncmp(keyW3, "4",1)==0){
					binType = 4;
				}else if (strncmp(keyW3, "8",1)==0){
					binType = 8;
				}
				//Binary data format
				printf("...reading data of binary (%d) format: %s \n", binType, input_path);
				// fread (&bSx[0],binType,1,FilePointer);
				// //printf("%f\n",nx[0]);
				// for (int k=0; k<znodes; k++){
				// 	for (int j=0; j<ynodes; j++){
				// 		for (int i=0; i<xnodes; i++){
				// 			int n = i + j*xnodes + k*xnodes*ynodes;
				// 			fread (&bSx[n],binType,1,FilePointer);
				// 			fread (&bSy[n],binType,1,FilePointer);
				// 			fread (&bSz[n],binType,1,FilePointer);
				// 			ctx->Sx[n]=bSx[n]; ctx->Sy[n]=bSy[n];ctx->Sz[n]=bSz[n];
				// 		}
				// 	}
				// }
				if(fread(&VEC_X(ctx->bS,0),binType,1,FilePointer)) {
				//printf("%f \n",bSx[0]);	
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
											VEC_X(ctx->S,I)=VEC_X(ctx->bS,I)=(double)temp4_x;
											VEC_Y(ctx->S,I)=VEC_Y(ctx->bS,I)=(double)temp4_y;
											VEC_Z(ctx->S,I)=VEC_Z(ctx->bS,I)=(double)temp4_z;
										}
									}else{
										if(!fread(&temp8_x,binType,1,FilePointer)) break;
										if(!fread(&temp8_y,binType,1,FilePointer)) break;
										if(!fread(&temp8_z,binType,1,FilePointer)) break;
										for (int t=0; t<ctx->AtomsPerBlock; t++){
											int I=n*ctx->AtomsPerBlock+t;
											VEC_X(ctx->S,I)=VEC_X(ctx->bS,I)=temp8_x;
											VEC_Y(ctx->S,I)=VEC_Y(ctx->bS,I)=temp8_y;
											VEC_Z(ctx->S,I)=VEC_Z(ctx->bS,I)=temp8_z;
										}
									}
								}	
							}
						}
					}
				}else{printf("problem\n");}
			}else{
				printf("Do not know what to do with \"%s\" data format in %s\n", keyW2, input_path);
			}
		}else{
			printf("%s has wrong data format or dimentionality!\n", input_path);
		}       
		// when everything is done
		printf("Done!\n");
		fclose(FilePointer);
		magnoom_reset_solver_state(ctx);
	}else{fprintf(stderr, "Cannot open input file '%s': %s\n", input_path, strerror(errno));}
    //metka dlya schiutyvaniya equilibrium state for dm
    for (int i=0; i<ctx->NOS; i++){
        VEC_X(ctx->t3S,i)=VEC_X(ctx->S,i);
        VEC_Y(ctx->t3S,i)=VEC_Y(ctx->S,i);
        VEC_Z(ctx->t3S,i)=VEC_Z(ctx->S,i);                       
    }
	ChangeVectorMode(ctx, 1);
}

void TW_CALL CB_ReadBIN( void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	char input_path[MAGNOOM_PATH_CAPACITY];
	

	unsigned short int num = 65535;



	int Nx = ctx->uABC[0], Ny = ctx->uABC[1], Nz = ctx->uABC[2];
	if (!magnoom_resolve_input_path(ctx, input_path, sizeof(input_path))) return;
	FILE * FilePointer = fopen(input_path, "rb");
	if (FilePointer == NULL) {
		fprintf(stderr, "Cannot open input file '%s': %s\n", input_path, strerror(errno));
		return;
	}
  	for(int k = 0; k<Nz; k++){
	    for(int j = 0; j<Ny; j++){
			for(int i = 0; i <Nx;i++){
				magnoom_bin_spin my_par_red;
				if (fread(&my_par_red, sizeof(my_par_red), 1, FilePointer)){
					double nx,ny,nz;
					unsigned short int p=my_par_red.t, q=my_par_red.f;

					double T = (double)(p+0.5)*PI/num;
					double F = (double)2*(q+0.5)*PI/num;

					nx = sin(T)*cos(F);
					ny = sin(T)*sin(F);
					nz = cos(T);
					int n = (i)+(j)*Nx+k*Nx*Ny;
					VEC_X(ctx->S,n) = nx; VEC_Y(ctx->S,n) = ny; VEC_Z(ctx->S,n) = nz;
					VEC_X(ctx->bS,n) = nx; VEC_Y(ctx->bS,n) = ny; VEC_Z(ctx->bS,n) = nz;
				}
			}
	    }
	}
	fclose (FilePointer);
	magnoom_reset_solver_state(ctx);
	ChangeVectorMode(ctx, 1);
}


void TW_CALL CB_ReadVTK( void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	char input_path[MAGNOOM_PATH_CAPACITY];
	if (!magnoom_resolve_input_path(ctx, input_path, sizeof(input_path))) return;
    Read_VTK(ctx, ctx->S, input_path);
	magnoom_reset_solver_state(ctx);
    ChangeVectorMode(ctx, 1);
}


void TW_CALL CB_Save_OVF_b8( void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	char output_path[MAGNOOM_PATH_CAPACITY];
	if (!magnoom_resolve_output_path_with_extension(ctx, ".ovf",
		output_path, sizeof(output_path))) return;
	Save_OVF_b8(ctx, ctx->bS, output_path);
}

void TW_CALL CB_Save_VTK_b4( void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	char output_path[MAGNOOM_PATH_CAPACITY];
	if (!magnoom_resolve_output_path_with_extension(ctx, ".vtk",
		output_path, sizeof(output_path))) return;
	Save_VTK(ctx, ctx->bS, 0, output_path);//metka 0->1

        // Save_VTS_b4(bSx, bSy, bSz, ctx->PosX, ctx->PosY, ctx->PosZ, Box, vts_filename);
        // Save_VTS_ascii(bSx, bSy, bSz, ctx->PosX, ctx->PosY, ctx->PosZ, Box, vts_filename);

        // float * Spins_xyz;
        // Spins_xyz = (float *)calloc(ctx->NOS*3, sizeof(float));   
        // for (int n=0; n<ctx->NOS; n++){
        //         Spins_xyz[3*n+0]=ctx->Sx[n];
        //         Spins_xyz[3*n+1]=ctx->Sy[n];
        //         Spins_xyz[3*n+2]=ctx->Sz[n];         
        // }
        //     save_vtk(vts_filename,"name",3,Spins_xyz,"special_flag",1,ctx->Kind,ctx->uABC[0],ctx->uABC[1],ctx->uABC[2],ctx->uABC[0],ctx->uABC[1],ctx->uABC[2],1);
        // 
}

void TW_CALL CB_Save_BIN( void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	char output_path[MAGNOOM_PATH_CAPACITY];
	if (!magnoom_resolve_output_path_with_extension(ctx, ".bin",
		output_path, sizeof(output_path))) return;
	SaveBin(ctx, ctx->bS, output_path);//metka 0->1
}

void TW_CALL CB_Save_PNG( void *clientData )
{
    magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	char output_path[MAGNOOM_PATH_CAPACITY];
	if (!magnoom_resolve_output_path_with_extension(ctx, ".png",
		output_path, sizeof(output_path))) return;
	SavePng(ctx, ctx->bS, output_path, ctx->WhichSliceMode, ctx->A_layer_min-1, ctx->B_layer_min-1, ctx->C_layer_min-1);//metka 0->1
}

static void magnoom_request_error_dialog(magnoom_ctx *ctx, const char *message)
{
	if (ctx->modal_bar != NULL || ctx->modal_open_requested) return;
	if (!magnoom_copy_path(ctx->modal_message, sizeof(ctx->modal_message), message))
		magnoom_copy_path(ctx->modal_message, sizeof(ctx->modal_message), "The file operation failed.");
	ctx->modal_close_requested = false;
	ctx->modal_open_requested = true;
}

#ifndef MAGNOOM_NO_MAIN
static void TW_CALL CB_CloseErrorDialog(void *clientData)
{
	magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	ctx->modal_close_requested = true;
}

static void magnoom_set_bar_visibility(TwBar *bar, int visible)
{
	TwSetParam(bar, NULL, "visible", TW_PARAM_INT32, 1, &visible);
}

static void magnoom_restore_modal_bars(magnoom_ctx *ctx)
{
	for (int i = 0; i < ctx->modal_saved_count; ++i)
		magnoom_set_bar_visibility(ctx->modal_saved_bars[i], ctx->modal_saved_visibility[i]);
	ctx->modal_saved_count = 0;
}

static void magnoom_open_error_dialog(magnoom_ctx *ctx)
{
	char definition[256];
	int framebuffer_width, framebuffer_height;
	int bar_width, bar_height, position_x, position_y;
	TwBar *bar;

	if (ctx->modal_bar != NULL) return;
	{
		int window_width = 0, window_height = 0;
		glfwGetWindowSize(&window_width, &window_height);
		framebuffer_width = (int)(window_width * ctx->ContentScaleX + 0.5);
		framebuffer_height = (int)(window_height * ctx->ContentScaleY + 0.5);
	}
	bar_width = framebuffer_width < 560 ? framebuffer_width : 560;
	bar_height = 200;
	if (bar_width < 220) bar_width = 220;
	position_x = (framebuffer_width - bar_width)/2;
	position_y = (framebuffer_height - bar_height)/2;
	if (position_x < 0) position_x = 0;
	if (position_y < 0) position_y = 0;

	bar = TwNewBar("FileError");
	if (bar == NULL) {
		fprintf(stderr, "%s\n", ctx->modal_message);
		return;
	}
	TwDefine(" FileError color='180 35 45' alpha=200");
	TwDefine(" FileError resizable=false movable=false iconifiable=false");
	snprintf(definition, sizeof(definition),
		" FileError label='File error' size='%d %d' position='%d %d' ",
		bar_width, bar_height, position_x, position_y);
	if (!TwDefine(definition) ||
		!TwAddButton(bar, "Message", NULL, NULL, "label='File operation failed.'") ||
		!TwSetParam(bar, "Message", "label", TW_PARAM_CSTRING, 1, ctx->modal_message) ||
		!TwAddButton(bar, "OK", CB_CloseErrorDialog, ctx, "label='OK'")) {
		TwDeleteBar(bar);
		fprintf(stderr, "%s\n", ctx->modal_message);
		return;
	}

	ctx->modal_bar = bar;
	ctx->modal_saved_count = 0;
	int bar_count = TwGetBarCount();
	for (int i = 0; i < bar_count && ctx->modal_saved_count < MAGNOOM_MODAL_BAR_CAPACITY; ++i) {
		TwBar *current = TwGetBarByIndex(i);
		int visible = 1;
		if (current == bar) continue;
		TwGetParam(current, NULL, "visible", TW_PARAM_INT32, 1, &visible);
		ctx->modal_saved_bars[ctx->modal_saved_count] = current;
		ctx->modal_saved_visibility[ctx->modal_saved_count] = visible;
		ctx->modal_saved_count++;
		magnoom_set_bar_visibility(current, 0);
	}
	TwSetTopBar(bar);
}

static void magnoom_service_modal(magnoom_ctx *ctx)
{
	if (ctx->modal_close_requested && ctx->modal_bar != NULL) {
		TwDeleteBar(ctx->modal_bar);
		ctx->modal_bar = NULL;
		ctx->modal_close_requested = false;
		magnoom_restore_modal_bars(ctx);
	}
	if (ctx->modal_open_requested && ctx->modal_bar == NULL) {
		ctx->modal_open_requested = false;
		magnoom_open_error_dialog(ctx);
	}
}
#endif

static bool magnoom_get_requested_file_format(magnoom_ctx *ctx, const char *filename,
	bool importing, FileFormatEnum *format)
{
	const char *extension = magnoom_get_file_extension(filename);
	char message[MAGNOOM_MODAL_MESSAGE_CAPACITY];

	if (extension == NULL) {
		snprintf(message, sizeof(message), "The file name has no extension.");
		magnoom_request_error_dialog(ctx, message);
		return false;
	}
	*format = GetFileFormatFromExtension(filename);
	if (*format == FILE_FORMAT_UNKNOWN) {
		snprintf(message, sizeof(message), "The file extension is not supported.");
		magnoom_request_error_dialog(ctx, message);
		return false;
	}
	if (importing && !magnoom_file_format_can_import(*format)) {
		snprintf(message, sizeof(message), "%s files are export-only and cannot be imported.",
			fileFormatNames[*format]);
		magnoom_request_error_dialog(ctx, message);
		return false;
	}
	if (!importing && !magnoom_file_format_can_export(*format)) {
		snprintf(message, sizeof(message), "%s files cannot be exported.", fileFormatNames[*format]);
		magnoom_request_error_dialog(ctx, message);
		return false;
	}
	return true;
}

void TW_CALL CB_Import(void *clientData)
{
	magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	char input_path[MAGNOOM_PATH_CAPACITY];
	char error_message[MAGNOOM_MODAL_MESSAGE_CAPACITY];
	FileFormatEnum format;

	magnoom_stop_engine(ctx);
	if (!magnoom_get_requested_file_format(ctx, ctx->inputfilename, true, &format)) return;
	if (!magnoom_resolve_input_path(ctx, input_path, sizeof(input_path))) {
		magnoom_request_error_dialog(ctx, "The input path is empty, unsupported, or too long.");
		return;
	}
	if (!ValidateFileFormat(ctx, input_path, format, error_message, sizeof(error_message))) {
		magnoom_request_error_dialog(ctx, error_message);
		return;
	}

	switch (format) {
		case FILE_FORMAT_CSV: CB_ReadCSV(ctx); break;
		case FILE_FORMAT_OVF: CB_ReadOVF(ctx); break;
		case FILE_FORMAT_VTK: CB_ReadVTK(ctx); break;
		case FILE_FORMAT_BIN: CB_ReadBIN(ctx); break;
		default: break;
	}
}

void TW_CALL CB_Export(void *clientData)
{
	magnoom_ctx *ctx = (magnoom_ctx *)clientData;
	FileFormatEnum format;

	if (!magnoom_get_requested_file_format(ctx, ctx->outputfilename, false, &format)) return;
	switch (format) {
		case FILE_FORMAT_CSV: CB_SaveCSV(ctx); break;
		case FILE_FORMAT_OVF: CB_Save_OVF_b8(ctx); break;
		case FILE_FORMAT_VTK: CB_Save_VTK_b4(ctx); break;
		case FILE_FORMAT_BIN: CB_Save_BIN(ctx); break;
		case FILE_FORMAT_PNG: CB_Save_PNG(ctx); break;
		default: break;
	}
}

void UpdateKind(magnoom_ctx *ctx)
{
	int *Kind = ctx->Kind;
	const float *Px = ctx->PosX;
	const float *Py = ctx->PosY;
	const float *Pz = ctx->PosZ;
	const int NOS = ctx->NOS;
	float dist, dist_max = ctx->chSizeG * ctx->chSizeG;
	ctx->NOSK = 0;
	switch(ctx->WhichGeometry){
		case CILINDER_G:
			for (int i=0; i<NOS; i++){
				dist = Px[i]*Px[i]+Py[i]*Py[i];
				if (dist>dist_max){
					Kind[i] = 0;
				}else{
					Kind[i] = 1;
					ctx->NOSK++;
				}
			}
		break;

		case SPHERE_G:
			for (int i=0; i<NOS; i++){
				dist = Px[i]*Px[i]+Py[i]*Py[i]+Pz[i]*Pz[i];
				if (dist>dist_max){
					Kind[i] = 0;
				}else{
					Kind[i] = 1;
					ctx->NOSK++;
				}
			}
		break;

		default:
			for (int i=0; i<NOS; i++){
					Kind[i] = 1;
				}
			ctx->NOSK = NOS;
		break;
	}

    for (int n=0; n<NOS; n++)
    {   
        VEC_X(ctx->S,n)*= Kind[n];
        VEC_Y(ctx->S,n)*= Kind[n];
        VEC_Z(ctx->S,n)*= Kind[n];

        VEC_X(ctx->bS,n)*= Kind[n];
        VEC_Y(ctx->bS,n)*= Kind[n];
        VEC_Z(ctx->bS,n)*= Kind[n];
        
        VEC_X(ctx->tS,n)*= Kind[n];
        VEC_Y(ctx->tS,n)*= Kind[n];
        VEC_Z(ctx->tS,n)*= Kind[n]; 
        
        VEC_X(ctx->t2S,n)*= Kind[n];
        VEC_Y(ctx->t2S,n)*= Kind[n];
        VEC_Z(ctx->t2S,n)*= Kind[n];  

        VEC_X(ctx->t3S,n)*= Kind[n];
        VEC_Y(ctx->t3S,n)*= Kind[n];
        VEC_Z(ctx->t3S,n)*= Kind[n]; 
    }
}

static bool ParseConfigFloat(const char *text, float *result)
{
	char *end = NULL;
	errno = 0;
	float value = strtof(text, &end);
	if (end == text || *end != '\0' || errno == ERANGE || !isfinite(value)) return false;
	*result = value;
	return true;
}

static bool ParseConfigInt(const char *text, int *result)
{
	char *end = NULL;
	errno = 0;
	long value = strtol(text, &end, 10);
	if (end == text || *end != '\0' || errno == ERANGE || value < INT_MIN || value > INT_MAX) return false;
	*result = (int)value;
	return true;
}

static bool ParseConfigDouble(const char *text, double *result)
{
	char *end = NULL;
	errno = 0;
	double value = strtod(text, &end);
	if (end == text || *end != '\0' || errno == ERANGE || !isfinite(value)) return false;
	*result = value;
	return true;
}

/* Trim leading/trailing whitespace from a NUL-terminated string in place. */
static char *ConfigTrim(char *text)
{
	while (isspace((unsigned char)*text)) ++text;
	if (*text == '\0') return text;
	char *end = text + strlen(text) - 1;
	while (end > text && isspace((unsigned char)*end)) --end;
	end[1] = '\0';
	return text;
}

/* Split a trimmed, non-empty, non-comment line into a trimmed key/value pair
 * around the first '='. Returns false if there is no '=' or either side is
 * empty after trimming. */
static bool SplitConfigKeyValue(char *line, char **key, char **value)
{
	char *equals = strchr(line, '=');
	if (equals == NULL) return false;
	*equals = '\0';
	*key = ConfigTrim(line);
	*value = ConfigTrim(equals + 1);
	return (*key)[0] != '\0' && (*value)[0] != '\0';
}

static bool ConfigError(const char *filename, int line_number, const char *format, ...)
{
	fprintf(stderr, "%s:%d: ", filename, line_number);
	va_list args;
	va_start(args, format);
	vfprintf(stderr, format, args);
	va_end(args);
	fprintf(stderr, "\n");
	return false;
}

typedef enum { CONFIG_FIELD_FLOAT, CONFIG_FIELD_INT } ConfigFieldType;

typedef struct ConfigScalarField {
	const char *key;
	ConfigFieldType type;
	size_t offset;
} ConfigScalarField;

static const ConfigScalarField CONFIG_SCALAR_FIELDS[] = {
	{"ax", CONFIG_FIELD_FLOAT, offsetof(magnoom_ctx, abc[0][0])},
	{"ay", CONFIG_FIELD_FLOAT, offsetof(magnoom_ctx, abc[0][1])},
	{"az", CONFIG_FIELD_FLOAT, offsetof(magnoom_ctx, abc[0][2])},
	{"bx", CONFIG_FIELD_FLOAT, offsetof(magnoom_ctx, abc[1][0])},
	{"by", CONFIG_FIELD_FLOAT, offsetof(magnoom_ctx, abc[1][1])},
	{"bz", CONFIG_FIELD_FLOAT, offsetof(magnoom_ctx, abc[1][2])},
	{"cx", CONFIG_FIELD_FLOAT, offsetof(magnoom_ctx, abc[2][0])},
	{"cy", CONFIG_FIELD_FLOAT, offsetof(magnoom_ctx, abc[2][1])},
	{"cz", CONFIG_FIELD_FLOAT, offsetof(magnoom_ctx, abc[2][2])},
	{"Na", CONFIG_FIELD_INT, offsetof(magnoom_ctx, uABC[0])},
	{"Nb", CONFIG_FIELD_INT, offsetof(magnoom_ctx, uABC[1])},
	{"Nc", CONFIG_FIELD_INT, offsetof(magnoom_ctx, uABC[2])},
	{"Shells", CONFIG_FIELD_INT, offsetof(magnoom_ctx, ShellNumber)},
	{"BCa", CONFIG_FIELD_INT, offsetof(magnoom_ctx, Boundary[0])},
	{"BCb", CONFIG_FIELD_INT, offsetof(magnoom_ctx, Boundary[1])},
	{"BCc", CONFIG_FIELD_INT, offsetof(magnoom_ctx, Boundary[2])},
};

static bool ApplyConfigScalarField(magnoom_ctx *ctx, const ConfigScalarField *field,
	const char *value_text)
{
	void *target = (char *)ctx + field->offset;
	if (field->type == CONFIG_FIELD_FLOAT) return ParseConfigFloat(value_text, (float *)target);
	return ParseConfigInt(value_text, (int *)target);
}

/* Match a leading "Atom<N>" or "AtomAll" prefix. On success, *atom is the
 * 0-based atom index (-1 for "All") and *rest points just past the prefix
 * (e.g. "x" for "Atom0x", "_K11" for "Atom0_K11"). */
static bool ParseAtomKeyPrefix(const char *key, int *atom, const char **rest)
{
	if (strncmp(key, "Atom", 4) != 0) return false;
	const char *after = key + 4;
	if (strncmp(after, "All", 3) == 0) {
		*atom = -1;
		*rest = after + 3;
		return true;
	}
	if (!isdigit((unsigned char)*after)) return false;
	char *end = NULL;
	long value = strtol(after, &end, 10);
	if (end == after || value < 0 || value >= MAX_ATOMS_PER_BLOCK) return false;
	*atom = (int)value;
	*rest = end;
	return true;
}

/* Handles one already-split "key = value" pair: scalar fields, Atom<N>x/y/z
 * positions (staged into `positions`, applied to ctx after the whole file is
 * read), and Atom<N|All>_K.. / Atom<N|All>_Q. anisotropy records (staged
 * directly into ctx->anisotropy_config_records[], same as before). */
static bool HandleConfigKeyValue(magnoom_ctx *ctx, const char *filename, int line_number,
	const char *key, const char *value,
	float positions[MAX_ATOMS_PER_BLOCK][3], bool position_seen[MAX_ATOMS_PER_BLOCK][3],
	int *max_atom_seen, bool *any_position_seen,
	bool quaternion_seen[MAX_ATOMS_PER_BLOCK + 1][4])
{
	for (size_t i = 0; i < sizeof(CONFIG_SCALAR_FIELDS)/sizeof(CONFIG_SCALAR_FIELDS[0]); ++i) {
		if (strcmp(key, CONFIG_SCALAR_FIELDS[i].key) == 0) {
			if (!ApplyConfigScalarField(ctx, &CONFIG_SCALAR_FIELDS[i], value))
				return ConfigError(filename, line_number, "invalid value for %s", key);
			return true;
		}
	}

	int atom;
	const char *rest;
	if (!ParseAtomKeyPrefix(key, &atom, &rest))
		return ConfigError(filename, line_number, "unknown configuration key '%s'", key);

	if (strcmp(rest, "x") == 0 || strcmp(rest, "y") == 0 || strcmp(rest, "z") == 0) {
		if (atom < 0)
			return ConfigError(filename, line_number,
				"AtomAll cannot be used for atom positions; use Atom<N>%s", rest);
		float parsed;
		if (!ParseConfigFloat(value, &parsed))
			return ConfigError(filename, line_number, "invalid value for %s", key);
		int component = rest[0] - 'x';
		positions[atom][component] = parsed;
		position_seen[atom][component] = true;
		*any_position_seen = true;
		if (atom > *max_atom_seen) *max_atom_seen = atom;
		return true;
	}

	if (rest[0] != '_')
		return ConfigError(filename, line_number, "unknown configuration key '%s'", key);
	const char *component = rest + 1;
	size_t component_len = strlen(component);

	AnisotropyConfigRecord record;
	memset(&record, 0, sizeof(record));
	record.line = line_number;
	record.atom = atom;

	if (component[0] == 'K' && (component_len == 3 || component_len == 5)) {
		int digit_count = (int)component_len - 1;
		for (int i = 0; i < digit_count; ++i) {
			char digit = component[1 + i];
			if (digit < '1' || digit > '3')
				return ConfigError(filename, line_number, "invalid tensor component '%s'", key);
			record.index[i] = digit - '1';
		}
		record.kind = digit_count == 2 ? ANISOTROPY_RECORD_K2 : ANISOTROPY_RECORD_K4;
		if (!ParseConfigDouble(value, &record.value))
			return ConfigError(filename, line_number, "invalid value for %s", key);
	} else if (strcmp(component, "Qx") == 0 || strcmp(component, "Qy") == 0 ||
		strcmp(component, "Qz") == 0 || strcmp(component, "Qs") == 0) {
		record.kind = ANISOTROPY_RECORD_QUATERNION;
		record.index[0] = component[1] == 'x' ? 0 : component[1] == 'y' ? 1 :
			component[1] == 'z' ? 2 : 3;
		if (!ParseConfigDouble(value, &record.value))
			return ConfigError(filename, line_number, "invalid value for %s", key);
		size_t slot = atom < 0 ? (size_t)MAX_ATOMS_PER_BLOCK : (size_t)atom;
		quaternion_seen[slot][record.index[0]] = true;
	} else {
		return ConfigError(filename, line_number, "unknown configuration key '%s'", key);
	}

	if (ctx->anisotropy_config_record_count >= MAX_ANISOTROPY_CONFIG_RECORDS)
		return ConfigError(filename, line_number, "too many anisotropy records (max %d)",
			MAX_ANISOTROPY_CONFIG_RECORDS);
	ctx->anisotropy_config_records[ctx->anisotropy_config_record_count++] = record;
	return true;
}

/* Write a fresh magnoom.cfg reflecting ctx's current (default) values, in the
 * key = value syntax readConfigFile() understands. */
bool writeDefaultConfigFile(const magnoom_ctx *ctx, const char *filename)
{
	FILE *FilePointer = fopen(filename, "wb");
	if (FilePointer == NULL) return false;

	fprintf(FilePointer, "# Lattice vectors\n");
	fprintf(FilePointer, "ax = %g\nay = %g\naz = %g\n",
		(double)ctx->abc[0][0], (double)ctx->abc[0][1], (double)ctx->abc[0][2]);
	fprintf(FilePointer, "bx = %g\nby = %g\nbz = %g\n",
		(double)ctx->abc[1][0], (double)ctx->abc[1][1], (double)ctx->abc[1][2]);
	fprintf(FilePointer, "cx = %g\ncy = %g\ncz = %g\n",
		(double)ctx->abc[2][0], (double)ctx->abc[2][1], (double)ctx->abc[2][2]);

	fprintf(FilePointer, "\n# Number of cells\n");
	fprintf(FilePointer, "Na = %d\nNb = %d\nNc = %d\n",
		ctx->uABC[0], ctx->uABC[1], ctx->uABC[2]);
	fprintf(FilePointer, "Shells = %d\n", ctx->ShellNumber);

	fprintf(FilePointer, "\n# Boundary conditions\n");
	fprintf(FilePointer, "BCa = %d\nBCb = %d\nBCc = %d\n",
		ctx->Boundary[0], ctx->Boundary[1], ctx->Boundary[2]);

	fprintf(FilePointer, "\n# Atom positions (fractional coordinates in the unit cell)\n");
	for (int atom = 0; atom < ctx->AtomsPerBlock; ++atom) {
		fprintf(FilePointer, "Atom%dx = %g\nAtom%dy = %g\nAtom%dz = %g\n",
			atom, (double)ctx->Block[atom][0],
			atom, (double)ctx->Block[atom][1],
			atom, (double)ctx->Block[atom][2]);
	}

	bool ok = ferror(FilePointer) == 0;
	if (fclose(FilePointer) != 0) ok = false;
	return ok;
}

bool readConfigFile(magnoom_ctx *ctx)
{
	const char configfilename[] = "magnoom.cfg";
	char raw_line[256];
	int line_number = 0;
	ctx->anisotropy_config_record_count = 0;

	FILE *FilePointer = fopen(configfilename, "rb");
	if (FilePointer == NULL) {
		if (writeDefaultConfigFile(ctx, configfilename)) {
			printf("%s not found; a default configuration file has been created.\n",
				configfilename);
		} else {
			printf("%s not found and could not be created; proceeding with built-in defaults.\n",
				configfilename);
		}
		return true;
	}

	float positions[MAX_ATOMS_PER_BLOCK][3];
	bool position_seen[MAX_ATOMS_PER_BLOCK][3];
	memset(position_seen, 0, sizeof(position_seen));
	bool any_position_seen = false;
	int max_atom_seen = -1;
	bool quaternion_seen[MAX_ATOMS_PER_BLOCK + 1][4];
	memset(quaternion_seen, 0, sizeof(quaternion_seen));

	bool first_line = true;
	while (fgets(raw_line, sizeof(raw_line), FilePointer) != NULL) {
		++line_number;
		if (strchr(raw_line, '\n') == NULL && !feof(FilePointer)) {
			ConfigError(configfilename, line_number, "line longer than %zu characters",
				sizeof(raw_line) - 1);
			fclose(FilePointer);
			return false;
		}

		char *line = ConfigTrim(raw_line);
		if (first_line) {
			first_line = false;
			if (strcmp(line, "# begin magnoom config") == 0) {
				ConfigError(configfilename, line_number,
					"uses the old comment-based syntax, which is no longer supported; "
					"convert it to 'key = value' form (see README.md)");
				fclose(FilePointer);
				return false;
			}
		}
		if (line[0] == '\0' || line[0] == '#') continue;

		char *key, *value;
		if (!SplitConfigKeyValue(line, &key, &value)) {
			ConfigError(configfilename, line_number, "malformed line (expected 'key = value')");
			fclose(FilePointer);
			return false;
		}

		if (!HandleConfigKeyValue(ctx, configfilename, line_number, key, value,
			positions, position_seen, &max_atom_seen, &any_position_seen, quaternion_seen)) {
			fclose(FilePointer);
			return false;
		}
	}

	if (ctx->uABC[0] <= 0 || ctx->uABC[1] <= 0 || ctx->uABC[2] <= 0 ||
		ctx->ShellNumber <= 0 || ctx->ShellNumber > MAX_SHELL_NUM) {
		fprintf(stderr, "%s must define positive dimensions and between 1 and %d shells.\n",
			configfilename, MAX_SHELL_NUM);
		fclose(FilePointer);
		return false;
	}

	/* magnoom_ctx_set_block() is also what recomputes NOB/NOS and the other
	 * uABC/abc-derived geometry fields, so it must run again whenever this
	 * file could have changed Na/Nb/Nc or the lattice vectors -- even if it
	 * didn't specify any Atom<N>x/y/z keys, in which case the atom positions
	 * already in ctx (magnoom_ctx_init()'s default block) are reapplied
	 * against the now-current uABC/abc. Skipping this call after a uABC
	 * change would leave NOB/NOS stale relative to the actual grid size. */
	int atom_count = ctx->AtomsPerBlock;
	const float (*basis_positions)[3] = ctx->Block;
	if (any_position_seen) {
		atom_count = max_atom_seen + 1;
		basis_positions = positions;
		for (int atom = 0; atom < atom_count; ++atom) {
			if (!position_seen[atom][0] || !position_seen[atom][1] || !position_seen[atom][2]) {
				fprintf(stderr,
					"%s: Atom%d is missing one of its x/y/z components; atom indices must be "
					"contiguous starting at Atom0.\n", configfilename, atom);
				fclose(FilePointer);
				return false;
			}
		}
	}
	if (!magnoom_ctx_set_block(ctx, atom_count, basis_positions)) {
		fprintf(stderr, "%s defines an invalid set of atom positions.\n", configfilename);
		fclose(FilePointer);
		return false;
	}

	for (int slot = 0; slot <= MAX_ATOMS_PER_BLOCK; ++slot) {
		bool any = quaternion_seen[slot][0] || quaternion_seen[slot][1] ||
			quaternion_seen[slot][2] || quaternion_seen[slot][3];
		bool all = quaternion_seen[slot][0] && quaternion_seen[slot][1] &&
			quaternion_seen[slot][2] && quaternion_seen[slot][3];
		if (any && !all) {
			if (slot == MAX_ATOMS_PER_BLOCK)
				fprintf(stderr, "%s: AtomAll is missing one or more of its Qx/Qy/Qz/Qs "
					"components; all four must be given together.\n", configfilename);
			else
				fprintf(stderr, "%s: Atom%d is missing one or more of its Qx/Qy/Qz/Qs "
					"components; all four must be given together.\n", configfilename, slot);
			fclose(FilePointer);
			return false;
		}
	}

	printf("Done!\n");
	fclose(FilePointer);
	return true;
}


// AntTweakBar lays out bars in raw pixel units with no DPI awareness;
// these convert logical sizes/positions to pixel units via the
// once-computed ContentScaleX/Y, replacing this project's former shim's
// tw_glfw2_set_bar_size/position/values_width helpers.
static void SetBarSize(magnoom_ctx *ctx, TwBar *bar, int width, int height)
{
	int size[2] = {
		(int)(width * ctx->ContentScaleX + 0.5),
		(int)(height * ctx->ContentScaleY + 0.5)
	};
	TwSetParam(bar, NULL, "size", TW_PARAM_INT32, 2, size);
}

static void SetBarPosition(magnoom_ctx *ctx, TwBar *bar, int x, int y)
{
	int position[2] = {
		(int)(x * ctx->ContentScaleX + 0.5),
		(int)(y * ctx->ContentScaleY + 0.5)
	};
	TwSetParam(bar, NULL, "position", TW_PARAM_INT32, 2, position);
}

static void SetBarValuesWidth(magnoom_ctx *ctx, TwBar *bar, int width)
{
	int scaled_width = (int)(width * ctx->ContentScaleX + 0.5);
	TwSetParam(bar, NULL, "valueswidth", TW_PARAM_INT32, 1, &scaled_width);
}

void setupTweakBar(magnoom_ctx *ctx)
{
/*  Global settings for the bar-menu  */
	TwDefine(" GLOBAL iconpos=topleft "); // icons go to top-left corner of the window
	TwDefine(" GLOBAL iconalign=horizontal "); // icons will be aligned horizontally
	TwDefine(" GLOBAL contained=true "); // bars cannot move outside of the window
	char iconMarginDef[64];
	snprintf(iconMarginDef, sizeof(iconMarginDef), "GLOBAL iconmargin='%d %d'",
	         (int)(ctx->ContentScaleX + 0.5), (int)(8.0 * ctx->ContentScaleY + 0.5));
	TwDefine(iconMarginDef);
/*  Help Bar F1 */
	ctx->help_bar = TwGetBarByIndex(0);
	TwDefine(" TW_HELP color='70 100 100'");
	SetBarSize(ctx, ctx->help_bar, 440, 530);
	TwDefine(" TW_HELP help='F1: show/hide (this) Help bar' "); // change default tweak bar size and color

/*  View Bar F2 */
    ctx->view_bar = TwNewBar("View");
    TwDefine(" View iconified=true "); 
    TwDefine(" View color='100 100 70' alpha=200 "); // change default tweak bar color
    SetBarSize(ctx, ctx->view_bar, 220, 530);
    TwDefine(" View help='F2: show/hide View bar' "); // change default tweak bar size and color

	// TwAddVarCB(ctx->view_bar, "Multisampling", TW_TYPE_INT32, CB_SetN_Multisample, CB_GetN_Multisample, ctx, " label='Multisamples' min=1 max=32 step=1 help='Multisampling' group='Camera'");
	{
	TwEnumVal		enProjectionsTw[] = { {ORTHO, "Orthogonal"}, {PERSP, "Perspective"} };
	TwType			TW_TYPE_PROJ = TwDefineEnum("ProjectionType", enProjectionsTw, 2);
	TwAddVarRW(ctx->view_bar, "Projection", TW_TYPE_PROJ, &ctx->WhichProjection, "keyIncr='p' help='Type of 3D projection' group='Camera'");
	}

    TwAddVarRW(ctx->view_bar, "ObjRotation", TW_TYPE_QUAT4F, &ctx->q_Rotation, " label='Scene rotation' opened=true help='Change the 3D scene orientation.' ");

	TwAddVarRW(ctx->view_bar, "CamAng", TW_TYPE_FLOAT, &ctx->PerspSet[0], " label='camera angle' min=1 max=120 help='camera angle' group='Camera'");
	TwAddVarRW(ctx->view_bar, "PosX", TW_TYPE_FLOAT, &ctx->TransXYZ[0], " label='position in x' min=-1000 max=1000 help='camera position along X-axis' group='Camera'");
	TwAddVarRW(ctx->view_bar, "PosY", TW_TYPE_FLOAT, &ctx->TransXYZ[1], " label='position in y' min=-1000 max=1000 help='camera position along Y-axis' group='Camera'");
	TwAddVarRW(ctx->view_bar, "PosZ", TW_TYPE_FLOAT, &ctx->TransXYZ[2], " label='position in z' min=-1000 max=1000 help='camera position along Z-axis' group='Camera'");

    TwAddVarCB(ctx->view_bar, "RotX", TW_TYPE_FLOAT, CB_SetRotX, CB_GetRotX, ctx, " label='turn around X' help='rotate camera around X-axis' group='Camera'");

    TwAddVarCB(ctx->view_bar, "RotY", TW_TYPE_FLOAT, CB_SetRotY, CB_GetRotY, ctx, " label='turn around Y' help='rotate camera around Y-axis' group='Camera'");

    TwAddVarCB(ctx->view_bar, "RotZ", TW_TYPE_FLOAT, CB_SetRotZ, CB_GetRotZ, ctx, " label='turn around Z' help='rotate camera around Z-axis' group='Camera'");

	TwAddVarRW(ctx->view_bar, "RotSpeed", TW_TYPE_FLOAT, &ctx->RotSpeed, " label='rotation speed' min=0 max=10 step=1 help='speed of rotation around any axis' group='Camera'");
	TwAddVarRW(ctx->view_bar, "TransSpeed", TW_TYPE_FLOAT, &ctx->TransSpeed, " label='translation speed' min=0 max=10 step=1 help='speed of translation along any axis' group='Camera'");

	TwDefine(" View/Camera opened=false ");

	// TwAddVarRW(view_bar, "CamBank", TW_TYPE_INT32, &CurrentCameraPositionBank, " label='Current camera' min=0 max=4 group='CameraRW'");
	// TwAddButton(view_bar, "Read Camera", CB_GetCameraPosition, ctx, "label='read camera pos.' group='CameraRW'");
	// TwAddButton(view_bar, "Write Camera", CB_SaveCameraPosition, ctx, "label='save camera pos.' group='CameraRW'");
	// TwDefine(" View/CameraRW opened=false ");

	/****** Light ******/
	{
	TwEnumVal enLightingModeTw[] = {
		{LIGHT_OFF, "Off"}, {LIGHT_FIXED, "Fixed"}, {LIGHT_ADAPTIVE, "Adaptive"}
	};
	TwType TW_TYPE_LIGHTING_MODE = TwDefineEnum("LightingMode", enLightingModeTw, 3);
	TwAddVarRW(ctx->view_bar, "Lighting", TW_TYPE_LIGHTING_MODE, &ctx->WhichLightingMode,
		" label='Lighting' keyIncr='l' help='Cycle lighting: Off, Fixed, Adaptive.' group='Light'");
	}
	TwAddVarRW(ctx->view_bar, "Intensity", TW_TYPE_FLOAT, &ctx->LightMultiplier, " label='Light intensity' min=0.1 max=2 step=0.02 help='Increase/decrease the light power.' group='Light' ");
	TwAddVarRW(ctx->view_bar, "LightDir", TW_TYPE_DIR3F, &ctx->LightDirection, " label='Fixed light direction' opened=false help='Change the direction used by Fixed lighting.' group='Light'");
	ctx->temp_color[0] = 230;
	ctx->temp_color[1] = 230;
	ctx->temp_color[2] = 255;
	TwSetParam(ctx->view_bar, "LightDir", "arrowcolor", TW_PARAM_INT32, 3, ctx->temp_color);


    TwDefine(" View/Light opened=false ");

	{
	TwEnumVal		enColorsTw[] = { {WHITE,"White"}, {BLACK, "Black"}, {RED, "Red"}, {GREEN, "Green"}, {BLUE, "Blue"}, {MANUAL, "Manual"} };
	TwType			TW_TYPE_COLOR = TwDefineEnum("BG_Color", enColorsTw, 6);
	TwAddVarRW(ctx->view_bar, "Choose_background", TW_TYPE_COLOR, &ctx->WhichBackgroundColor, "help='Background color for 3D scene' group='Background'");
	}
	TwAddVarRW(ctx->view_bar, "red", TW_TYPE_FLOAT, &ctx->BackgroundColors[5][0], 
	" min=0 max=1 step=0.01 group='Background'");
	TwAddVarRW(ctx->view_bar, "green", TW_TYPE_FLOAT, &ctx->BackgroundColors[5][1], 
	" min=0 max=1 step=0.01 group='Background'");
	TwAddVarRW(ctx->view_bar, "blue", TW_TYPE_FLOAT, &ctx->BackgroundColors[5][2], 
	" min=0 max=1 step=0.01 group='Background'");
	TwDefine(" View/Background opened=false ");


	{
	TwEnumVal		enVectorModeTw[] = {{ARROW1, "Arrows"} 
										,{CONE1,   "Cones"} 
										,{CANE,    "Canes"} 
										,{uPOINT,  "Points"}
										,{BOX1,    "Boxes"}
									};
	TwType			TW_TYPE_VEC_MOD = TwDefineEnum("Type_of_vectors", enVectorModeTw, 5);
	TwAddVarCB(ctx->view_bar, "Type of vectors", TW_TYPE_VEC_MOD, CB_SetVectorMode, CB_GetVectorMode, ctx, "keyIncr='v' keyDecr='V' help='Type of 3D vectors' group='Appearance' ");
	}

	TwAddVarCB(ctx->view_bar, "Pivot", TW_TYPE_FLOAT, CB_SetPivot, CB_GetPivot,  ctx, " min=0 max=1 step=0.01 help='Pivot of 3D arrow.' group='Appearance' ");
	TwAddVarCB(ctx->view_bar, "Faces", TW_TYPE_INT32, CB_SetFaces, CB_GetFaces,  ctx, " min=3 max=20 step=1 help='Number of faces for 3D arrow.' group='Appearance' ");
	TwAddVarCB(ctx->view_bar, "Scale", TW_TYPE_FLOAT, CB_SetScale, CB_GetScale,  ctx, " min=0.1 max=10 step=0.01 keyIncr='+' keyDecr='-' help='Scale the vectors.' group='Appearance' ");
	TwAddVarCB(ctx->view_bar, "BextDCArrowScale", TW_TYPE_FLOAT, CB_SetScaleBextDC, CB_GetScaleBextDC, ctx,
		"label='Bext DC arrow scale' min=0.1 max=10 step=0.01 keyIncr='+' keyDecr='-' help='Scale the rendered static external-field arrow.' group='Appearance'");


	TwAddVarRW(ctx->view_bar, "Show basis", TW_TYPE_BOOL32, &ctx->AxesOn, " key=CTRL+o  group='Appearance'");
	TwAddVarRW(ctx->view_bar, "Show box", TW_TYPE_BOOL32, &ctx->BoxOn, " key=CTRL+b group='Appearance'");
	TwDefine(" View/Appearance opened=false ");

	{
	TwEnumVal		enSliceModeTw[] = { 
										{A_AXIS, "a-axis"}, 
										{B_AXIS, "b-axis"}, 
										{C_AXIS, "c-axis"},
										{FILTER, "filter"}
									  };
	TwType			TV_TYPE_VEC_MOD = TwDefineEnum("Slicing", enSliceModeTw, 4);
	TwAddVarCB(ctx->view_bar, "Slicing mode", TV_TYPE_VEC_MOD, CB_SetSliceMode, CB_GetSliceMode, ctx, "keyIncr='/' keyDecr='?' help='Slising plane perpenticulat to the choosen axis' group='Filters&Slices' ");
	}
	TwAddVarCB(ctx->view_bar, "GreedFilter", TW_TYPE_BOOL32, CB_SetGreedFilter, CB_GetGreedFilter, ctx, 
		" label='Greed filter' true='On' false='Off' group='Filters&Slices' ");
	TwAddVarCB(ctx->view_bar, "Spin_filter1", TW_TYPE_BOOL32, CB_SetSpinFilter1, CB_GetSpinFilter1, ctx, 
		" label='Spin filter N1' true='On' false='Off' group='Filters&Slices' ");
	TwAddVarCB(ctx->view_bar, "Spin_filter2", TW_TYPE_BOOL32, CB_SetSpinFilter2, CB_GetSpinFilter2, ctx, 
		" label='Spin filter N2' true='On' false='Off' group='Filters&Slices' ");
	TwAddVarCB(ctx->view_bar, "Spin_filter3", TW_TYPE_BOOL32, CB_SetSpinFilter3, CB_GetSpinFilter3, ctx, 
		" label='Spin filter N3' true='On' false='Off' group='Filters&Slices' ");

	TwAddVarCB(ctx->view_bar, "T_max1", TW_TYPE_INT32, CB_SetThetaMax1, CB_GetThetaMax1, ctx, 
		" group='Spin_filter_1' label='Theta max 1' min=0 max=180  help='max value for polar angle theta'");
	TwAddVarCB(ctx->view_bar, "T_min1", TW_TYPE_INT32, CB_SetThetaMin1, CB_GetThetaMin1, ctx, 
		" group='Spin_filter_1' label='Theta min 1' min=0 max=180  help='min value for polar angle theta'");
	TwAddVarCB(ctx->view_bar, "F_max1", TW_TYPE_INT32, CB_SetPhiMax1, CB_GetPhiMax1, ctx, 
		" group='Spin_filter_1' label='Phi max 1' min=0 max=360 step=1 help='max value for azimuthal angle phi'");
	TwAddVarCB(ctx->view_bar, "F_min1", TW_TYPE_INT32, CB_SetPhiMin1, CB_GetPhiMin1, ctx, 
		" group='Spin_filter_1' label='Phi min 1' min=0 max=360 step=1 help='min value for azimuthal angle phi'");

	TwAddVarCB(ctx->view_bar, "T_max2", TW_TYPE_INT32, CB_SetThetaMax2, CB_GetThetaMax2, ctx, 
		" group='Spin_filter_2' label='Theta max 2' min=0 max=180  help='max value for polar angle theta'");
	TwAddVarCB(ctx->view_bar, "T_min2", TW_TYPE_INT32, CB_SetThetaMin2, CB_GetThetaMin2, ctx, 
		" group='Spin_filter_2' label='Theta min 2' min=0 max=180  help='min value for polar angle theta'");	
	TwAddVarCB(ctx->view_bar, "F_max2", TW_TYPE_INT32, CB_SetPhiMax2, CB_GetPhiMax2, ctx, 
		" group='Spin_filter_2' label='Phi max 2' min=0 max=360 step=1 help='max value for azimuthal angle phi'");
	TwAddVarCB(ctx->view_bar, "F_min2", TW_TYPE_INT32, CB_SetPhiMin2, CB_GetPhiMin2, ctx, 
		" group='Spin_filter_2' label='Phi min 2' min=0 max=360 step=1 help='min value for azimuthal angle phi'");

	TwAddVarCB(ctx->view_bar, "T_max3", TW_TYPE_INT32, CB_SetThetaMax3, CB_GetThetaMax3, ctx, 
		" group='Spin_filter_3' label='Theta max 3' min=0 max=180  help='max value for polar angle theta'");
	TwAddVarCB(ctx->view_bar, "T_min3", TW_TYPE_INT32, CB_SetThetaMin3, CB_GetThetaMin3, ctx, 
		" group='Spin_filter_3' label='Theta min 3' min=0 max=180  help='min value for polar angle theta'");	
	TwAddVarCB(ctx->view_bar, "F_max3", TW_TYPE_INT32, CB_SetPhiMax3, CB_GetPhiMax3, ctx, 
		" group='Spin_filter_3' label='Phi max 3' min=0 max=360 step=1 help='max value for azimuthal angle phi'");
	TwAddVarCB(ctx->view_bar, "F_min3", TW_TYPE_INT32, CB_SetPhiMin3, CB_GetPhiMin3, ctx, 
		" group='Spin_filter_3' label='Phi min 3' min=0 max=360 step=1 help='min value for azimuthal angle phi'");

	TwDefine(" View/Spin_filter_1 opened=false group='Filters&Slices'");
	TwDefine(" View/Spin_filter_2 opened=false group='Filters&Slices'");
	TwDefine(" View/Spin_filter_3 opened=false group='Filters&Slices'");

	TwAddVarCB(ctx->view_bar, "GFmaxA", TW_TYPE_INT32, CB_SetGreedMaxA, CB_GetGreedMaxA, ctx, " group='Greed_filter' label='max na' ");
	TwAddVarCB(ctx->view_bar, "GFminA", TW_TYPE_INT32, CB_SetGreedMinA, CB_GetGreedMinA, ctx, " group='Greed_filter' label='min na' ");
	TwAddVarCB(ctx->view_bar, "GFmaxB", TW_TYPE_INT32, CB_SetGreedMaxB, CB_GetGreedMaxB, ctx, " group='Greed_filter' label='max nb' ");
	TwAddVarCB(ctx->view_bar, "GFminB", TW_TYPE_INT32, CB_SetGreedMinB, CB_GetGreedMinB, ctx, " group='Greed_filter' label='min nb' ");
	TwAddVarCB(ctx->view_bar, "GFmaxC", TW_TYPE_INT32, CB_SetGreedMaxC, CB_GetGreedMaxC, ctx, " group='Greed_filter' label='max nc' ");
	TwAddVarCB(ctx->view_bar, "GFminC", TW_TYPE_INT32, CB_SetGreedMinC, CB_GetGreedMinC, ctx, " group='Greed_filter' label='min nc' ");

    TwAddVarCB(ctx->view_bar, "GreedFilterInvert", TW_TYPE_BOOL32, CB_SetGreedFilterInvert, CB_GetGreedFilterInvert, ctx, " label='Invert G filter' true='On' false='Off' group='Greed_filter' ");


	TwDefine(" View/Greed_filter opened=false group='Filters&Slices'");
	TwDefine(" View/Filters&Slices opened=false ");


	TwAddVarCB(ctx->view_bar, "ColorShift", TW_TYPE_INT32, CB_SetColorShift, CB_GetColorShift, ctx, " label='Rotate hue' min=0 max=360 help='rotate color hue in xy-plane' group='HSV to RGB'");
	TwAddVarCB(ctx->view_bar, "InvHue", TW_TYPE_BOOL32, CB_SetInvHue, CB_GetInvHue, ctx, " label='Invert hue' help='invert RGB to RBG color hue' group='HSV to RGB'");
	TwAddVarCB(ctx->view_bar, "InvVal", TW_TYPE_BOOL32, CB_SetInvVal, CB_GetInvVal, ctx, " label='Invert value' help='invert black to white' group='HSV to RGB'");
	TwDefine(" View/'HSV to RGB' opened=false ");

/*  Hamiltonian parameters&controls F3 */
	ctx->control_bar = TwNewBar("Parameters&Controls");
    TwDefine(" Parameters&Controls iconified=true "); 
    TwDefine(" Parameters&Controls color='100 70 100' alpha=200 "); // change default tweak bar color
    SetBarSize(ctx, ctx->control_bar, 220, 530);
    TwDefine(" Parameters&Controls help='F3: show/hide Control bar' "); // change default tweak bar size and color

	//TwAddButton(control_bar, "Run", CB_Run, ctx, "key='SPACE' label='RUN simulation' ");
	TwAddVarCB(ctx->control_bar, "Run", TW_TYPE_BOOL32, CB_Set_Run, CB_Get_Run, ctx, " keyIncr='SPACE' true='On' false='Off' label='RUN simulation' ");
	TwAddVarRW(ctx->control_bar, "Record", TW_TYPE_BOOL32, &ctx->Record, "label='Recording' true='On' false='Off' help='Recording <sx>, <sy>, <sz> in each iteration'");
	TwAddVarRW(ctx->control_bar, "Rec_Iteration", TW_TYPE_INT32, &ctx->rec_iteration, "label='Every iteration' min=1 max=1000 step=1 ");
    TwAddVarRW(ctx->control_bar, "Max_Iteration", TW_TYPE_INT32, &ctx->Max_Numb_Iteration, "label='Max. Iteration' min=1 max=100000000 step=100 ");
	TwAddButton(ctx->control_bar, "Clean the record", CB_CleanSxSySzFile, ctx, "label= 'Clean the record' help='clean the output file with <sx>, <sy>, <sz>' ");
	TwAddButton(ctx->control_bar, "Reset iterations", CB_ResetIterations, ctx, "label='Reset iterations' ");

	TwAddSeparator(ctx->control_bar, "sep-3", NULL);

	TwAddVarRW(ctx->control_bar, "BCinA", TW_TYPE_BOOL32, &ctx->Boundary[0], "label='along a' group='Boundary conditions' true='periodic' false='open' help='set boundary conditions along translation vector 'a' '");
	TwAddVarRW(ctx->control_bar, "BCinB", TW_TYPE_BOOL32, &ctx->Boundary[1], "label='along b' group='Boundary conditions' true='periodic' false='open' help='set boundary conditions along translation vector 'b' '");
	TwAddVarRW(ctx->control_bar, "BCinC", TW_TYPE_BOOL32, &ctx->Boundary[2], "label='along c' group='Boundary conditions' true='periodic' false='open' help='set boundary conditions along translation vector 'c' '");
    
	{
    TwEnumVal       enIntegrationScheme[] = {{HEUN, "Heun(1) "},
                                             {SIB,  " SIB(2) "},
                                             {RK2, " RK2 (midpoint, fixed step) "},
                                             {RK4, " RK4 (classical, fixed step) "},
                                             {RELAX, " RELAX "}};
    TwType          TV_TYPE_INTEGRATION_SCHEME = TwDefineEnum("Solver", enIntegrationScheme, 5);
    TwAddVarRW(ctx->control_bar, "Integration scheme", TV_TYPE_INTEGRATION_SCHEME, &ctx->WhichIntegrationScheme, "group='LLG' help='Choose the integration scheme'");
    }
    //TwAddVarRW(control_bar, "Preces", TW_TYPE_BOOL32, &Precession, "label='precession' group='LLG' true='On' false='Off' help='On/Off precession'");
    TwAddVarRW(ctx->control_bar, "Damping", TW_TYPE_FLOAT, &ctx->damping, "label='Damping' min=0 max=100 step=0.000001 group='LLG' ");
	TwAddVarRW(ctx->control_bar, "Time_step", TW_TYPE_FLOAT, &ctx->t_step, "label='Time step' min=0 max=10.0 step=0.000001   group='LLG' ");
	TwAddVarRW(ctx->control_bar, "temperature", TW_TYPE_FLOAT, &ctx->Temperature, "label='k_b*T' min=0 max=100 step=0.000001 group='LLG' ");

	TwAddSeparator(ctx->control_bar, "sep-4", NULL);
	TwAddVarRW(ctx->control_bar, "Xi", TW_TYPE_FLOAT, &ctx->Xi, "label='Xi'  step=0.001 help='Non-adiabaticity parameter' ");
	TwAddVarRW(ctx->control_bar, "Curr_u", TW_TYPE_FLOAT, &ctx->Curr_u, "label='Curr_u'  step=0.001 help='Current in Zhang-Li torque' ");

    TwAddSeparator(ctx->control_bar, "sep-2", NULL);
	TwAddVarCB(ctx->control_bar, "BextDCTheta", TW_TYPE_FLOAT, CB_SetBextDCTheta, CB_GetBextDCTheta, ctx, "label='Bext DC theta' step=0.0001 help='Change the static external-field direction' ");
	TwAddVarCB(ctx->control_bar, "BextDCPhi", TW_TYPE_FLOAT, CB_SetBextDCPhi, CB_GetBextDCPhi, ctx, "label='Bext DC phi' step=0.0001 help='Change the static external-field direction' ");
	TwAddVarCB(ctx->control_bar, "BextDCMagnitude", TW_TYPE_FLOAT, CB_SetBextDCMagnitude, CB_GetBextDCMagnitude, ctx, "label='Bext DC magnitude' min=0 step=0.0000001 help='Change the static external-field magnitude' ");
	// TwAddVarCB(control_bar, "BextDCDirection", TW_TYPE_DIR3F, CB_SetBextDCDirection, CB_GetBextDCDirection, BextDCDirection,
	// "label='Bext DC direction' opened=true help='Change the static external-field direction' ");
	// temp_color[0] = 55;
	// temp_color[1] = 55;
	// temp_color[2] = 155;
	// TwSetParam(control_bar, "BextDCDirection", "arrowcolor", TW_PARAM_INT32, 3, temp_color);
	// // TwAddVarRW(control_bar, "BextDCMagnitude", TW_TYPE_FLOAT, &BextDCMagnitude,
	// // "label='Bext DC magnitude' help='Static external-field magnitude' ");
	// TwAddSeparator(control_bar, "control_sep1", NULL);
	// TwAddVarCB(control_bar, "BextDCDirection", TW_TYPE_DIR3F, CB_SetBextDCDirection, CB_GetBextDCDirection, ctx,
	// "label='Bext DC direction' opened=false help='Change the static external-field direction.' ");

	TwAddSeparator(ctx->control_bar, "control_sep2", NULL);
	for(int s=0; s<ctx->ShellNumber; s++)
	{
	snprintf(ctx->ShortBuffer,200,"J%1i",s);
	TwAddVarRW(ctx->control_bar, ctx->ShortBuffer, TW_TYPE_FLOAT, &ctx->Jij[s], "help='Heisenberg exchange' ");
	TwSetParam(ctx->control_bar, ctx->ShortBuffer, "label",  TW_PARAM_CSTRING , 1, ctx->ShortBuffer);
	}
	//////////////////////////////////////////
	TwAddSeparator(ctx->control_bar, "sep2", NULL);
	//////////////////////////////////////////	
	for(int s=0; s<ctx->ShellNumber; s++)
	{
	snprintf(ctx->ShortBuffer,200,"B%1i",s);
	TwAddVarRW(ctx->control_bar, ctx->ShortBuffer, TW_TYPE_FLOAT, &ctx->Bij[s], "help='Bi-quadratic exchange' ");
	TwSetParam(ctx->control_bar, ctx->ShortBuffer, "label",  TW_PARAM_CSTRING , 1, ctx->ShortBuffer);
	}
	//////////////////////////////////////////
	TwAddSeparator(ctx->control_bar, "sep3", NULL);
	//////////////////////////////////////////
	for(int s=0; s<ctx->ShellNumber; s++)
	{
	snprintf(ctx->ShortBuffer,200,"D%1i",s);
	TwAddVarRW(ctx->control_bar, ctx->ShortBuffer, TW_TYPE_FLOAT, &ctx->Dij[s], "help='Dzyaloshinskii-Moriya' ");
	TwSetParam(ctx->control_bar, ctx->ShortBuffer, "label",  TW_PARAM_CSTRING , 1, ctx->ShortBuffer);
	}
	//////////////////////////////////////////
	TwAddSeparator(ctx->control_bar, "sep4", NULL);
	//////////////////////////////////////////
	TwAddVarRW(ctx->control_bar, "ST param.", TW_TYPE_FLOAT, &ctx->Cu, "help='Dzyaloshinskii-Moriya' ");
	TwAddVarRW(ctx->control_bar, "CurrentDir", TW_TYPE_DIR3F, &ctx->VCu, "label='cur. dir.' opened=true help='The polarization direction of electric current' ");

/*  Initial state F4 */
	ctx->initial_bar = TwNewBar("Initial_State");
	TwDefine(" Initial_State iconified=true "); 
	TwDefine(" Initial_State color='70 70 100' alpha=200"); // change default tweak bar color
	SetBarSize(ctx, ctx->initial_bar, 220, 530);
	TwDefine(" Initial_State help='F4: show/hide Initial state bar' "); // change default tweak bar size and color
/*	{
	TwEnumVal		enGeomTw[] = { 	{DEFAULT_G, 	"Default"		    }, 
									{CILINDER_G, 	"Cilinder"	        }, 
									{SPHERE_G, 		"Sphere"	        }
									};
	TwType			TV_TYPE_GEOMETRY = TwDefineEnum("DomainShape", enGeomTw, 3);
	TwAddVarRW(initial_bar, "Choose shape", TV_TYPE_GEOMETRY, &WhichGeometry, "help='Choose shape of the simulated domain'");
	}
	TwAddVarRW(initial_bar, "Size", TW_TYPE_FLOAT,  &chSizeG, 
	" min=0 max=100000 step=0.5 help='characteristic size of the shape (radius)' ");

	TwAddButton(initial_bar, "Set shape", CB_SetShape, ctx, " label='Set shape' ");

	TwAddSeparator(initial_bar, "sep00", NULL);
*/
	{
	TwEnumVal		enIniStateTw[] = { 	{RND, 		"Random"		        }, 
										{HOMO, 		"Homogeneous"	        }, 
										{SKYRM1, 	"Skyrmion Q=1"	        }, 
										{SKYRM2, 	"Skyrmion Q=2"	        }, 
										{SKYRM3, 	"Skyrmion Q=3"	        }, 
										{BOBBER_T, 	"Bobber top"	        },  
										{BOBBER_B, 	"Bobber bottom"	        },
										{BOBBER_L,	"Bobber lattice"	    }, 
										{BOBBER_L_T,"Bobber latt. top"	    },  
										{BOBBER_L_B,"Bobber latt. bottom"	}, 
										{HOPFION1, 	"Hopfion"		        }, 
										{SPIRAL, 	"Spiral"		        }, 
										{SKYRMION_L,"Sk. lattice"	        },
										{GLOBULA,   "Globula"	            },
										{MultyQ,    "Multy-Q"				},
										{NORM,      "Normalize all"			},
									};
	TwType			TV_TYPE_INI_STATE = TwDefineEnum("IniState", enIniStateTw, 16);
	TwAddVarRW(ctx->initial_bar, "Choose ini. state", TV_TYPE_INI_STATE, &ctx->WhichInitialState, "help='Choose initial spin configuration'");
	}

	TwAddVarRW(ctx->initial_bar, "chSize", TW_TYPE_FLOAT,  &ctx->chSize, " min=-100000 max=100000 step=0.5 help='characteristic size of modulated state: Skyrmion diameter or spiral period' ");

	TwAddVarRW(ctx->initial_bar, "chDir", TW_TYPE_DIR3F,  &ctx->chDir, " help='direction of modulations e.g. k-vector of the spin spiral.' ");

	TwAddButton(ctx->initial_bar, "Set initial", CB_SetInitial, ctx, "key=I label='insert initial state' ");

	TwAddSeparator(ctx->initial_bar, "sep01", NULL);
	TwAddVarRW(ctx->initial_bar, "Degrees", TW_TYPE_FLOAT,  &ctx->RotateAllSpins, " min=-360 max=360 step=1 help='Rotate all spins about characteristic direction' ");
	TwAddButton(ctx->initial_bar, "Rotate spins", CB_RotateAllSpins, ctx, "label='rotate all spins' ");

	TwAddSeparator(ctx->initial_bar, "sep1", NULL);
	TwAddButton(ctx->initial_bar, "Invert X", CB_InvertX, ctx, "label='invert n_x component' ");
	TwAddButton(ctx->initial_bar, "Invert Y", CB_InvertY, ctx, "label='invert n_y component' ");
	TwAddButton(ctx->initial_bar, "Invert Z", CB_InvertZ, ctx, "label='invert n_z component' ");

	TwAddSeparator(ctx->initial_bar, "sep02", NULL);
	TwAddVarCB(ctx->initial_bar, "Input directory:", TW_TYPE_CSSTRING(sizeof(ctx->input_directory)),
		CB_SetInputDirectory, CB_GetInputDirectory, ctx,
		"help='Directory used for relative input file names; absolute file names override it'");
	TwAddVarRW(ctx->initial_bar, "Input file name:", TW_TYPE_CSSTRING(sizeof(ctx->inputfilename)), ctx->inputfilename, "");
	TwAddButton(ctx->initial_bar, "Import", CB_Import, ctx,
		"label='Import' help='Import CSV, OVF, VTK, or BIN according to the input file extension'");

	TwAddSeparator(ctx->initial_bar, "sep2", NULL);
	TwAddVarRW(ctx->initial_bar, "Save slice", TW_TYPE_BOOL32, &ctx->save_slice, " label='Save slice' help='Save current slice only' ");
	{
	TwEnumVal		enAverage_mode[] = {{ALONG_A, 	"Along a-axis"	}, 
										{ALONG_B, 	"Along b-axis"	}, 
										{ALONG_C, 	"Along c-axis"	}, 
										{ALONG_0, 	"No averaging"	}
									};

	TwType			TV_TYPE_AVERAGE_MODE = TwDefineEnum("Average mode", enAverage_mode, 4);
	TwAddVarRW(ctx->initial_bar, "Averaging mode", TV_TYPE_AVERAGE_MODE, &ctx->WhichAverageMode, "help='Choose type of average mode'");}

	TwAddVarCB(ctx->initial_bar, "Output directory:", TW_TYPE_CSSTRING(sizeof(ctx->output_directory)),
		CB_SetOutputDirectory, CB_GetOutputDirectory, ctx,
		"help='Directory for all simulation outputs; changing it starts a new table.csv'");
	TwAddVarRW(ctx->initial_bar, "Output file name:", TW_TYPE_CSSTRING(sizeof(ctx->outputfilename)), ctx->outputfilename, ""); 
	TwAddButton(ctx->initial_bar, "Export", CB_Export, ctx,
		"label='Export' help='Export CSV, OVF, VTK, BIN, or PNG according to the output file extension'");

    TwAddSeparator(ctx->initial_bar, "sep_isoline", NULL);

/*  Time-dependent external field F5 */
	ctx->BextAC_bar = TwNewBar("BextAC");
	TwDefine(" BextAC iconified=true ");
	TwDefine(" BextAC color='100 70 70' alpha=200"); // change default tweak bar color
	SetBarSize(ctx, ctx->BextAC_bar, 220, 530);
	TwDefine(" BextAC help='F5: show/hide Bext AC bar' "); // change default tweak bar size and color
	TwAddVarRW(ctx->BextAC_bar, "BextACEnabled", TW_TYPE_BOOL32, &ctx->BextACEnabled, "keyIncr='f' label='Bext AC enabled' true='on' false='off' help='Enable the time-dependent external-field component'");
    TwAddVarRW(ctx->BextAC_bar, "BextACModeRecording", TW_TYPE_BOOL32, &ctx->BextACModeRecording, "label='Mode recording' true='on' false='off' help='Record the phase-resolved dynamical mode'");
    TwAddVarRW(ctx->BextAC_bar, "Average Dynamical Mode", TW_TYPE_INT32, &ctx->rec_num_mode, "label='Num Mode Average' help='Number over which dynamical mode is averaged'");
    TwAddVarCB(ctx->BextAC_bar, "ModeRecordingPhaseCount", TW_TYPE_INT32, CB_SetNumImages, CB_GetNumImages, ctx, "label='Phase count' min=1 help='Number of phase bins used for dynamical-mode recording' ");

	TwAddVarRW(ctx->BextAC_bar, "BextACAmplitude", TW_TYPE_FLOAT, &ctx->BextACAmplitude, "label='Bext AC amplitude' help='Set the time-dependent external-field amplitude' ");

	TwAddVarRW(ctx->BextAC_bar, "BextACDirection", TW_TYPE_DIR3F, &ctx->BextACDirection, "label='Bext AC direction' opened=true help='Change the time-dependent external-field direction' ");

	TwAddSeparator(ctx->BextAC_bar, "sep0", NULL);

	ctx->temp_color[0] = 55;
	ctx->temp_color[1] = 55;
	ctx->temp_color[2] = 155;

	TwSetParam(ctx->BextAC_bar, "BextACDirection", "arrowcolor", TW_PARAM_INT32, 3, ctx->temp_color);
	{
	TwEnumVal       enBextACWaveformTw[] = { {BEXT_AC_SIN, "Sin(omega*t)"},
	                                    {BEXT_AC_GAUSSIAN, "Gaussian pulse"},
	                                    {BEXT_AC_SINC, "Sinc(omega*[t-offset])"},
	                                    {BEXT_AC_CIRCULAR, "Circular (omega*t)"}};
	TwType          TV_TYPE_INI_STATE = TwDefineEnum("BextACWaveform", enBextACWaveformTw, 4);
	TwAddVarRW(ctx->BextAC_bar, "BextACWaveform", TV_TYPE_INI_STATE, &ctx->BextACWaveform, "label='Bext AC waveform' help='Choose the time-dependent external-field waveform'");
	}

	TwAddVarRW(ctx->BextAC_bar, "BextACTimeOffset", TW_TYPE_FLOAT, &ctx->BextACTimeOffset, "label='Time offset' min=0 help='Set the pulse time offset' ");
	TwAddVarRW(ctx->BextAC_bar, "BextACPulseWidth", TW_TYPE_FLOAT, &ctx->BextACPulseWidth, "label='Pulse width/cutoff' min=0 help='Set the Gaussian width or sinc cutoff time' ");
	TwAddVarCB(ctx->BextAC_bar, "BextACPeriod", TW_TYPE_DOUBLE, CB_SetBextACPeriod, CB_GetBextACPeriod, ctx, "label='Bext AC period' min=0 help='Set the time-dependent field period' ");
	TwAddVarCB(ctx->BextAC_bar, "BextACOmega", TW_TYPE_DOUBLE, CB_SetBextACOmega, CB_GetBextACOmega, ctx, "label='Bext AC omega' min=0 help='Set omega; omega=2*pi/period' ");


/*  General tensor anisotropy F6 */
	ctx->anisotropy_bar = TwNewBar("Anisotropy");
	TwDefine(" Anisotropy iconified=true ");
	TwDefine(" Anisotropy color='70 100 100' alpha=200 ");
	SetBarSize(ctx, ctx->anisotropy_bar, 260, 530);
	TwDefine(" Anisotropy help='F6: show/hide tensor anisotropy bar' ");
	{
		TwEnumVal modes[] = {
			{ANISOTROPY_GLOBAL, "Global"},
			{ANISOTROPY_INDIVIDUAL, "Individual"}
		};
		TwType mode_type = TwDefineEnum("AnisotropyMode", modes, 2);
		TwAddVarCB(ctx->anisotropy_bar, "Mode", mode_type,
			CB_SetAnisotropyMode, CB_GetAnisotropyMode, ctx,
			"label='Mode' help='Global uses atom 0 for every spin; Individual uses per-atom tensors'");
	}
	{
		TwEnumVal atoms[MAX_ATOMS_PER_BLOCK];
		char labels[MAX_ATOMS_PER_BLOCK][16];
		for (int atom = 0; atom < ctx->AtomsPerBlock; ++atom) {
			snprintf(labels[atom], sizeof(labels[atom]), "Atom %d", atom);
			atoms[atom].Value = atom;
			atoms[atom].Label = labels[atom];
		}
		TwType atom_type = TwDefineEnum("AnisotropyAtom", atoms, (unsigned int)ctx->AtomsPerBlock);
		TwAddVarCB(ctx->anisotropy_bar, "Atom", atom_type,
			CB_SetAnisotropyAtom, CB_GetAnisotropyAtom, ctx,
			"label='Atom' help='Unit-cell atom edited in Individual mode'");
	}
	TwAddVarCB(ctx->anisotropy_bar, "Rotation", TW_TYPE_QUAT4D,
		CB_SetAnisotropyQuaternion, CB_GetAnisotropyQuaternion, ctx,
		"label='Rotation' opened=true showval=true help='Rotate the selected atom local anisotropy tensor (all four components are editable)'");
	TwAddVarRW(ctx->anisotropy_bar, "Axis", TW_TYPE_DIR3F, &ctx->anisotropy_rotation_axis,
		" help='Unit axis of the additional rotation composed with the current quaternion' ");
	TwAddVarRW(ctx->anisotropy_bar, "Angle", TW_TYPE_FLOAT, &ctx->anisotropy_rotation_angle,
		" min=-360 max=360 step=1 help='Angle (degrees) of the additional rotation' ");
	TwAddButton(ctx->anisotropy_bar, "Compose rotation", CB_AnisotropyApplyAxisAngle, ctx,
		" label='compose axis-angle' help='Compose the current quaternion with the axis-angle rotation' ");
	TwAddButton(ctx->anisotropy_bar, "Reset rotation", CB_AnisotropyResetQuaternion, ctx,
		" label='reset to identity' help='Reset the rotation quaternion to identity (0,0,0,1)' ");
	TwAddButton(ctx->anisotropy_bar, "CopyAtom0", CB_CopyAnisotropyAtom0, ctx,
		"label='Copy atom 0 to all' help='Copy atom 0 local K2/K4 tensors to all atoms while preserving quaternions'");

	for (int component = 0; component < ANISOTROPY_K2_COMPONENT_COUNT; ++component) {
		AnisotropyComponentControl *control = &ctx->anisotropy_component_controls[component];
		const int *index = anisotropy_k2_components[component];
		char definition[160];
		control->ctx = ctx;
		control->kind = ANISOTROPY_COMPONENT_K2;
		control->component = component;
		snprintf(definition, sizeof(definition),
			"label='K_%d%d' step=0.000001 precision=9 group='Rank 2'",
			index[0] + 1, index[1] + 1);
		TwAddVarCB(ctx->anisotropy_bar, anisotropy_k2_control_names[component], TW_TYPE_FLOAT,
			CB_SetAnisotropyComponent, CB_GetAnisotropyComponent, control, definition);
	}
	for (int component = 0; component < ANISOTROPY_K4_COMPONENT_COUNT; ++component) {
		AnisotropyComponentControl *control =
			&ctx->anisotropy_component_controls[ANISOTROPY_K2_COMPONENT_COUNT + component];
		const int *index = anisotropy_k4_components[component];
		char definition[160];
		control->ctx = ctx;
		control->kind = ANISOTROPY_COMPONENT_K4;
		control->component = component;
		snprintf(definition, sizeof(definition),
			"label='K_%d%d%d%d' step=0.000001 precision=9 group='Rank 4'",
			index[0] + 1, index[1] + 1, index[2] + 1, index[3] + 1);
		TwAddVarCB(ctx->anisotropy_bar, anisotropy_k4_control_names[component], TW_TYPE_FLOAT,
			CB_SetAnisotropyComponent, CB_GetAnisotropyComponent, control, definition);
	}
	anisotropy_refresh_controls(ctx);
	TwAddButton(ctx->anisotropy_bar, "ExportAnisotropyMap", CB_ExportAnisotropyMap, ctx,
		"label='Anisotropy map to file' help='Sample the anisotropy energy over (theta, phi) and write CSV + legacy VTK surfaces for every atom and the total'");


/*  Info bar F12 */
	ctx->info_bar = TwNewBar("Info");
	TwDefine(" Info refresh=0.5 ");
	TwDefine(" Info iconified = false movable = false alwaysbottom=true resizable=false fontstyle=default fontsize=2"); 
	TwDefine(" Info help='F12: show/hide info-bar' "); // change default tweak bar size and color
	TwDefine(" Info color='10 10 10' alpha=0 "); // change default tweak bar size and color
	SetBarSize(ctx, ctx->info_bar, 170, 620);
	SetBarPosition(ctx, ctx->info_bar, 1, 30);
	SetBarValuesWidth(ctx, ctx->info_bar, 130);
	TwAddVarRO(ctx->info_bar, "Run/Stop", TW_TYPE_BOOL32,  &ctx->Play, "true='RUNING' false='STOPED' ");
	TwAddVarRO(ctx->info_bar, "RECORD", TW_TYPE_BOOL32,  &ctx->Record, "true='On' false='Off' ");
	TwAddSeparator(ctx->info_bar, "sep+21", NULL);
	//TwAddVarRO(info_bar, "BextACEnabled", TW_TYPE_BOOL32, &BextACEnabled, "true='On' false='Off'");
	TwAddButton(ctx->info_bar, "Bext AC", NULL, NULL, " ");
	TwAddVarRO(ctx->info_bar, " BextAC x ", TW_TYPE_DOUBLE, &ctx->BextAC[0], " ");
	TwAddVarRO(ctx->info_bar, " BextAC y ", TW_TYPE_DOUBLE, &ctx->BextAC[1], " ");
	TwAddVarRO(ctx->info_bar, " BextAC z ", TW_TYPE_DOUBLE, &ctx->BextAC[2], " ");

	TwAddSeparator(ctx->info_bar, "sep+11", NULL);

	//TwAddVarRO(info_bar, "Bext DC", TW_TYPE_FLOAT, NULL, "label='Bext DC'");
	TwAddButton(ctx->info_bar, "Bext DC", NULL, NULL, " ");
	TwAddVarRO(ctx->info_bar, " BextDC x", TW_TYPE_FLOAT, &ctx->BextDC[0], " ");
	TwAddVarRO(ctx->info_bar, " BextDC y", TW_TYPE_FLOAT, &ctx->BextDC[1], " ");
	TwAddVarRO(ctx->info_bar, " BextDC z", TW_TYPE_FLOAT, &ctx->BextDC[2], " ");

	TwAddSeparator(ctx->info_bar, "sep-0", NULL);
	TwAddVarRO(ctx->info_bar, "NPB", TW_TYPE_INT32,  &ctx->AtomsPerBlock, "help='number of atoms per block' ");
	TwAddVarRO(ctx->info_bar, "N_a", TW_TYPE_INT32,  &ctx->uABC[0],       "help='translations along a' ");
	TwAddVarRO(ctx->info_bar, "N_b", TW_TYPE_INT32,  &ctx->uABC[1],       "help='translations along b' ");
	TwAddVarRO(ctx->info_bar, "N_c", TW_TYPE_INT32,  &ctx->uABC[2],       "help='translations along c' ");	
	TwAddVarRO(ctx->info_bar, "NOS", TW_TYPE_INT32,  &ctx->NOS,           "help='Number of spins'      ");
	TwAddSeparator(ctx->info_bar, "sep", NULL);
	TwAddVarRO(ctx->info_bar, "ITR", TW_TYPE_INT32,  &ctx->currentIteration, "help='Total number of iterations' ");


	TwAddSeparator(ctx->info_bar, "sep0", NULL);

	TwAddVarRO(ctx->info_bar, "FPS", TW_TYPE_FLOAT,  &ctx->FPS, "help='Frame per second' precision=4");
	TwAddVarRO(ctx->info_bar, "IPS", TW_TYPE_FLOAT,  &ctx->IPS, "help='Iterations per secon' precision=4");	

	TwAddSeparator(ctx->info_bar, "sep01", NULL);


	TwAddVarRO(ctx->info_bar, "Etot", TW_TYPE_DOUBLE, &ctx->totalEnergy, " precision=10 help='Total energy' ");
	TwAddVarRO(ctx->info_bar, "etot", TW_TYPE_DOUBLE, &ctx->perSpEnergy, " precision=10 help='Energy density per spin'");
	TwAddVarRO(ctx->info_bar, "esat", TW_TYPE_DOUBLE, &ctx->totalEnergyFerro, " precision=10 help='Energy density for ferromagnetic state' ");
	TwAddVarRO(ctx->info_bar, "de", TW_TYPE_DOUBLE, &ctx->perSpEnergyMinusFerro, " precision=10 help='Energy density per spin wrt ferromagnetic state'");

	TwAddSeparator(ctx->info_bar, "sep1", NULL);	
	TwAddVarRO(ctx->info_bar, "M_x", TW_TYPE_DOUBLE, &ctx->Mtot[0], " help='x-component of total moment' precision=10");
	TwAddVarRO(ctx->info_bar, "M_y", TW_TYPE_DOUBLE, &ctx->Mtot[1], " help='y-component of total moment' precision=10");
	TwAddVarRO(ctx->info_bar, "M_z", TW_TYPE_DOUBLE, &ctx->Mtot[2], " help='z-component of total moment' precision=10");

	TwAddSeparator(ctx->info_bar, "sep2", NULL);
	TwAddVarRO(ctx->info_bar, "m_x", TW_TYPE_DOUBLE, &ctx->mtot[0], " help='x-component of average moment per spin' precision=10");
	TwAddVarRO(ctx->info_bar, "m_y", TW_TYPE_DOUBLE, &ctx->mtot[1], " help='y-component of average moment per spin' precision=10");
	TwAddVarRO(ctx->info_bar, "m_z", TW_TYPE_DOUBLE, &ctx->mtot[2], " help='z-component of average moment per spin' precision=10");
	TwAddSeparator(ctx->info_bar, "sep3", NULL);
	TwAddVarRO(ctx->info_bar, "Torque", TW_TYPE_DOUBLE, &ctx->MAX_TORQUE, " help='maximum torque acting on the spin' precision=10");
}

static int GLFWSpecialToWindowKey(int key)
{
	switch (key) {
		case GLFW_KEY_UP: return WINDOW_KEY_UP;
		case GLFW_KEY_DOWN: return WINDOW_KEY_DOWN;
		case GLFW_KEY_F1: return WINDOW_KEY_F1;
		case GLFW_KEY_F2: return WINDOW_KEY_F2;
		case GLFW_KEY_F3: return WINDOW_KEY_F3;
		case GLFW_KEY_F4: return WINDOW_KEY_F4;
		case GLFW_KEY_F5: return WINDOW_KEY_F5;
		case GLFW_KEY_F6: return WINDOW_KEY_F6;
		case GLFW_KEY_F12: return WINDOW_KEY_F12;
		default: return 0;
	}
}

// GLFW2's native key callback reports no modifier bitmask (unlike GLFW3),
// so Shift state is read directly via glfwGetKey, same as this project's
// former shim's own internal tw_glfw2_mods() helper did.
static int WindowShiftIsDown(void)
{
	return glfwGetKey(GLFW_KEY_LSHIFT) == GLFW_PRESS || glfwGetKey(GLFW_KEY_RSHIFT) == GLFW_PRESS;
}

void GLFWKeyCallback(int key, int action)
{
	magnoom_ctx *ctx = g_ctx;
	if (ctx->modal_open_requested) {
		if (action == GLFW_RELEASE) (void)TwEventKeyGLFW(key, action);
		return;
	}
	if (TwEventKeyGLFW(key, action)) return;
	if (ctx->modal_bar != NULL) return;

	int special = GLFWSpecialToWindowKey(key);
	if (special && action != GLFW_RELEASE) {
		HandleSpecialKey(ctx, special);
		return;
	}
	if (key == GLFW_KEY_ESC && action != GLFW_RELEASE) {
		HandleKeyDown(ctx, TW_KEY_ESCAPE);
		return;
	}
	if (action == GLFW_RELEASE && key >= 'A' && key <= 'Z') {
		unsigned char character = (unsigned char)(WindowShiftIsDown() ? key : key - 'A' + 'a');
		HandleKeyUp(ctx, character);
	}
}

// GLFW2's native char callback fires on both press and release (like the
// key callback), so this filters to press only, matching GLFW3-style
// char-input semantics (and this project's former shim, which did the
// same filtering internally).
void GLFWCharCallback(int character, int action)
{
	magnoom_ctx *ctx = g_ctx;
	if (action != GLFW_PRESS) return;
	if (ctx->modal_open_requested) return;
	if (character > 255 || TwEventCharGLFW(character, GLFW_PRESS)) return;
	if (ctx->modal_bar != NULL) return;
	HandleKeyDown(ctx, (unsigned char)character);
}

void GLFWMouseButtonCallback(int button, int action)
{
	magnoom_ctx *ctx = g_ctx;
	if (ctx->modal_open_requested) {
		if (action == GLFW_RELEASE) (void)TwEventMouseButtonGLFW(button, action);
		return;
	}
	if (TwEventMouseButtonGLFW(button, action)) return;
	if (ctx->modal_bar != NULL) return;
	int x = 0, y = 0;
	glfwGetMousePos(&x, &y);
	int windowButton = button == GLFW_MOUSE_BUTTON_LEFT ? WINDOW_MOUSE_LEFT
	                 : button == GLFW_MOUSE_BUTTON_RIGHT ? WINDOW_MOUSE_RIGHT
	                 : WINDOW_MOUSE_MIDDLE;
	HandleMouseButton(ctx, windowButton, action == GLFW_PRESS ? WINDOW_BUTTON_DOWN : WINDOW_BUTTON_UP,
	            (int)(x * ctx->ContentScaleX), (int)(y * ctx->ContentScaleY));
}

void GLFWCursorPosCallback(int x, int y)
{
	magnoom_ctx *ctx = g_ctx;
	int pixelX = (int)(x * ctx->ContentScaleX);
	int pixelY = (int)(y * ctx->ContentScaleY);
	if (ctx->modal_open_requested) return;
	if (TwEventMousePosGLFW(pixelX, pixelY)) return;
	if (ctx->modal_bar != NULL) return;
	HandleMouseDrag(ctx, pixelX, pixelY);
}

// GLFW2's native wheel callback reports an absolute position (unlike
// GLFW3's per-event delta), so the delta used below is the change since
// the last call; ctx->MouseWheelPosition itself is now just that raw
// absolute position (previously a value reconstructed from shim-faked
// deltas).
void GLFWScrollCallback(int position)
{
	magnoom_ctx *ctx = g_ctx;
	int delta = position - ctx->MouseWheelPosition;
	ctx->MouseWheelPosition = position;
	if (ctx->modal_open_requested) return;
	if (TwEventMouseWheelGLFW(ctx->MouseWheelPosition)) return;
	if (ctx->modal_bar != NULL) return;
	ctx->TransXYZ[2] += (float)delta * 0.5f;
}

// Registered as GLFW2's window-size callback (there is no separate
// framebuffer-size concept to register for); width/height arrive in
// window/logical units and are converted to pixel units via
// ContentScaleX/Y. The MeasureContentScale() call is a no-op here once
// contentScaleLocked is set (see setupOpenGL()): by the time the user can
// actually trigger a real resize, ContentScaleX/Y is already the fixed
// ratio established during startup, and must stay fixed rather than be
// re-derived from the window size this same call is about to change.
void GLFWWindowSizeCallback(int width, int height)
{
	magnoom_ctx *ctx = g_ctx;
	MeasureContentScale(ctx);
	int pixelWidth = (int)(width * ctx->ContentScaleX + 0.5);
	int pixelHeight = (int)(height * ctx->ContentScaleY + 0.5);
	if (pixelWidth < 1) pixelWidth = 1;
	if (pixelHeight < 1) pixelHeight = 1;
	ApplyFramebufferSize(ctx, pixelWidth, pixelHeight);
}


void HandleMouseButton(magnoom_ctx *ctx, int button, int state, int x, int y)
{
	int button_mask;
	switch (button) {
		case WINDOW_MOUSE_LEFT:   button_mask = LEFT;   break;
		case WINDOW_MOUSE_MIDDLE: button_mask = MIDDLE; break;
		case WINDOW_MOUSE_RIGHT:  button_mask = RIGHT;  break;
		default: return;
	}

	if (state == WINDOW_BUTTON_DOWN) {
		ctx->Xmouse = x;
		ctx->Ymouse = y;
		ctx->ActiveButton |= button_mask;
	} else {
		ctx->ActiveButton &= ~button_mask;
	}
}


void HandleMouseDrag(magnoom_ctx *ctx, int x, int y)
{
	const int dx = x - ctx->Xmouse;
	const int dy = y - ctx->Ymouse;

	if (ctx->ActiveButton & LEFT) {
		if (ctx->angle != 0) {
			float x_rotation[4];
			float y_rotation[4];
			SetQuaternionFromAxisAngle(ctx->axisY, dx*0.01, y_rotation);
			SetQuaternionFromAxisAngle(ctx->axisX, dy*0.01, x_rotation);
			MultiplyQuaternions(y_rotation, x_rotation, x_rotation);
			MultiplyQuaternions(x_rotation, ctx->q_Rotation, ctx->q_Rotation);
		}
	}

	if (ctx->ActiveButton & MIDDLE) {
		ctx->TransXYZ[2] += (dx-dy)*0.05;
	}

	if (ctx->ActiveButton & RIGHT) {
		ctx->TransXYZ[0] += dx*0.1;
		ctx->TransXYZ[1] -= dy*0.1;
	}

	ctx->Xmouse = x;
	ctx->Ymouse = y;
}

// the keyboard callback:
void HandleKeyDown( magnoom_ctx *ctx, unsigned char key )
{   
    {
      // your code here to handle the event	
	  // if( false ) fprintf( stderr, "Keyboard: '%c' (0x%0x)\n", c, c );
		switch( key )
		{
			case '<':
				switch(ctx->WhichSliceMode)
                {case A_AXIS:
					if (ctx->A_layer_min>1) {ctx->A_layer_min-=1; ChangeVectorMode(ctx, 0);}
				 break;
				 case B_AXIS:
				    if (ctx->B_layer_min>1) {ctx->B_layer_min-=1; ChangeVectorMode(ctx, 0);}
				 break;
				 case C_AXIS:
				 	if (ctx->C_layer_min>1) {ctx->C_layer_min-=1; ChangeVectorMode(ctx, 0);}
				 break;	
				 default:
				 break;		 
				}
			break;

			case '>':
				switch(ctx->WhichSliceMode)
                {case A_AXIS:
					if (ctx->A_layer_min<ctx->A_layer_max) {ctx->A_layer_min+=1; ChangeVectorMode(ctx, 0);}
				 break;
				 case B_AXIS:
				    if (ctx->B_layer_min<ctx->B_layer_max) {ctx->B_layer_min+=1; ChangeVectorMode(ctx, 0);}
				 break;
				 case C_AXIS:
				 	if (ctx->C_layer_min<ctx->C_layer_max) {ctx->C_layer_min+=1; ChangeVectorMode(ctx, 0);}
				 break;	
				 default:
				 break;		 
				}
			break;

			case ',':
				switch(ctx->WhichSliceMode)
                {case A_AXIS:
					if (ctx->A_layer_max>ctx->A_layer_min) {ctx->A_layer_max-=1; ChangeVectorMode (ctx, 0);}
				 break;
				 case B_AXIS:
				    if (ctx->B_layer_max>ctx->B_layer_min) {ctx->B_layer_max-=1; ChangeVectorMode (ctx, 0);}
				 break;
				 case C_AXIS:
				 	if (ctx->C_layer_max>ctx->C_layer_min) {ctx->C_layer_max-=1; ChangeVectorMode (ctx, 0);}
				 break;
				 default:
				 break;			 
				}
			break;

			case '.':
				switch(ctx->WhichSliceMode)
                {case A_AXIS:
					if (ctx->A_layer_max<ctx->uABC[0]) {ctx->A_layer_max+=1; ChangeVectorMode (ctx, 0);}
				 break;
				 case B_AXIS:
				    if (ctx->B_layer_max<ctx->uABC[1]) {ctx->B_layer_max+=1; ChangeVectorMode (ctx, 0);}
				 break;
				 case C_AXIS:
				 	if (ctx->C_layer_max<ctx->uABC[2]) {ctx->C_layer_max+=1; ChangeVectorMode (ctx, 0);}
				 break;	
				 default:
				 break;		 
				}
			break;

			case 'q':
			case 'Q':
				ctx->dRot[2] = - ctx->RotSpeed*0.5;
			break;

			case 'e':
			case 'E':
				ctx->dRot[2] = ctx->RotSpeed*0.5;
			break;

			case 'w':
			case 'W':
			    ctx->dRot[0] = - ctx->RotSpeed*0.5;
				break;
			case 's':
			case 'S':
				ctx->dRot[0] = ctx->RotSpeed*0.5;
				//TransXYZ[2]+=0.5;
				break;
			case 'a':
			case 'A':
				ctx->dRot[1] = - ctx->RotSpeed*0.5;
				//TransXYZ[0]-=1;	
				break;
			case 'd':
			case 'D':
				ctx->dRot[1] = ctx->RotSpeed*0.5;
				//TransXYZ[0]+=1;	
				break;

			case 'G':
			case 'g':
				ctx->dTransXYZ[2] = -ctx->TransSpeed*0.1;	
			break;

			case 'T':
			case 't':
				ctx->dTransXYZ[2] = ctx->TransSpeed*0.1;
			break;

			case 'x':
			case 'X':
				ExecuteCommand(ctx,  XUP );
			break;

			case 'y':
			case 'Y':
				ExecuteCommand(ctx,  YUP );
			break;

			case 'z':
			case 'Z':
				ExecuteCommand(ctx,  ZUP );
			break;

			case TW_KEY_ESCAPE:
				ExecuteCommand(ctx,  QUIT );
				break;

			//default:
				//fprintf( stderr, "Don't know what to do with keyboard hit: '%c' (0x%0x)\n", c, c );
		}

	}
}


void HandleKeyUp( magnoom_ctx *ctx, unsigned char key )
{
		switch( key )
		{			

			case 'q':
			case 'Q':
				ctx->dRot[2] = 0.0f;
			break;

			case 'e':
			case 'E':
				ctx->dRot[2] = 0.0f;
			break;

			case 'w':
			case 'W':
				//ctx->Rot[0] -= 0.25;//NSK
			    ctx->dRot[0] = 0.0f;
				//TransXYZ[2]-=0.5;
				break;
			case 's':
			case 'S':
				ctx->dRot[0] = 0.0f;
				//TransXYZ[2]+=0.5;
				break;
			case 'a':
			case 'A':
				ctx->dRot[1] = 0.0f;
				//TransXYZ[0]-=1;	
				break;
			case 'd':
			case 'D':
				ctx->dRot[1] = 0.0f;
				//TransXYZ[0]+=1;	
				break;

			case 'G':
			case 'g':
				ctx->dTransXYZ[2] = 0.0f;
			break;


			case 'T':
			case 't':
				ctx->dTransXYZ[2] = 0.0f;	
			break;
		}
}


void HandleSpecialKey( magnoom_ctx *ctx, int key ){
int isiconified;
    {
      // your code here to handle the event	
		switch( key )
		{
			case WINDOW_KEY_UP:
				switch (ctx->WhichSliceMode)
				{case A_AXIS:
					if (ctx->A_layer_max < ctx->uABC[0]) {ctx->A_layer_max+=1; ctx->A_layer_min+=1; ChangeVectorMode(ctx, 1);}
				 break;
				 case B_AXIS:
				    if (ctx->B_layer_max < ctx->uABC[1]) {ctx->B_layer_max+=1; ctx->B_layer_min+=1; ChangeVectorMode(ctx, 1);}
				 break;
				 case C_AXIS:
				 	if (ctx->C_layer_max < ctx->uABC[2]) {ctx->C_layer_max+=1; ctx->C_layer_min+=1; ChangeVectorMode(ctx, 1);}
				 break;	
				 default:
				 break;		 
				}
				break;
			case WINDOW_KEY_DOWN:
				switch (ctx->WhichSliceMode)
				{case A_AXIS:
					if (ctx->A_layer_min>1) {ctx->A_layer_min-=1; ctx->A_layer_max-=1; ChangeVectorMode(ctx, 1);}
				 break;
				 case B_AXIS:
				    if (ctx->B_layer_min>1) {ctx->B_layer_min-=1; ctx->B_layer_max-=1; ChangeVectorMode(ctx, 1);}
				 break;
				 case C_AXIS:
				 	if (ctx->C_layer_min>1) {ctx->C_layer_min-=1; ctx->C_layer_max-=1; ChangeVectorMode(ctx, 1);} 
				 break;
				 default:
				 break;
				}
				break;
			case  WINDOW_KEY_F1:
				TwGetParam(ctx->help_bar, NULL, "iconified", TW_PARAM_INT32, 1, &isiconified);
				if (isiconified){
			    	TwDefine(" TW_HELP iconified=false ");				
				}else{
					TwDefine(" TW_HELP iconified=true ");
				}
				break;
			case  WINDOW_KEY_F2:
				TwGetParam(ctx->view_bar, NULL, "iconified", TW_PARAM_INT32, 1, &isiconified);
				if (isiconified){
					TwDefine(" View iconified=false ");				
				}else{
					TwDefine(" View iconified=true ");
				}
				break;
			case  WINDOW_KEY_F3:
				TwGetParam(ctx->control_bar, NULL, "iconified", TW_PARAM_INT32, 1, &isiconified);
				if (isiconified){
					TwDefine(" Parameters&Controls iconified=false ");				
				}else{
					TwDefine(" Parameters&Controls iconified=true ");
				}
				break;
			case  WINDOW_KEY_F4:
				TwGetParam(ctx->initial_bar, NULL, "iconified", TW_PARAM_INT32, 1, &isiconified);
				if (isiconified){
					TwDefine(" Initial_State iconified=false ");				
				}else{
					TwDefine(" Initial_State iconified=true ");
				}
				break;
			case  WINDOW_KEY_F5:
				TwGetParam(ctx->BextAC_bar, NULL, "iconified", TW_PARAM_INT32, 1, &isiconified);
				if (isiconified){
					TwDefine(" BextAC iconified=false ");
				}else{
					TwDefine(" BextAC iconified=true ");
				}
				break;
			case  WINDOW_KEY_F6:
				TwGetParam(ctx->anisotropy_bar, NULL, "iconified", TW_PARAM_INT32, 1, &isiconified);
				if (isiconified){
					TwDefine(" Anisotropy iconified=false ");
				}else{
					TwDefine(" Anisotropy iconified=true ");
				}
				break;
			case  WINDOW_KEY_F12:
				TwGetParam(ctx->info_bar, NULL, "iconified", TW_PARAM_INT32, 1, &isiconified);
				if (isiconified){
					TwDefine(" Info iconified=false ");				
				}else{
					TwDefine(" Info iconified=true ");
				}
				break;
			//default:
				//fprintf( stderr, "Don't know what to do with special key: '%c' (0x%0x)\n", key, key );
		}
	}
}

void ExecuteCommand( magnoom_ctx *ctx, int id )
{
	switch( id )
	{
		case XUP:
			Xup(ctx);
			break;

		case YUP:
			Yup(ctx);
			break;

		case ZUP:
			Zup(ctx);
			break;

		case QUIT:
			ctx->WindowShouldClose = true;
			break;

		// default:
		// 	fprintf( stderr, "Don't know what to do with Button ID %d\n", id );
	}
}

void ChangeInitialState(magnoom_ctx *ctx)
{
	InitSpinComponents(ctx, ctx->PosX, ctx->PosY, ctx->PosZ, ctx->S, ctx->WhichInitialState);
	magnoom_reset_solver_state(ctx);
	ChangeVectorMode (ctx,  1 );
	ctx->SpecialEvent=1;
}

void ChangeVectorMode( magnoom_ctx *ctx, int id )
{
	switch ( id )
	{
		case 0: // change of mode e.g. from point to arrow
		ReallocateArrayDrawing(ctx);
		// Fill array for prototype (arrow or cane) array 
		UpdateProtoVerNorInd_Spins(ctx);
		// Fill big array for indecies for all arrows, cans, cones or boxes 
		UpdateIndices(ctx);
		case 1:	// change only the layer(s) for drawing
		UpdateVerticesNormalsColors(ctx);
		// FILTER slice mode uses N_filter subset
		if (ctx->WhichSliceMode == FILTER) {
			ctx->spin_mesh.index_count  = ctx->N_filter * ctx->IdNumProto;
			ctx->spin_mesh.component_count = ctx->N_filter * ctx->VCNumProto;
		} else {
			ctx->spin_mesh.index_count  = ctx->IdNum;
			ctx->spin_mesh.component_count = ctx->VCNum;
		}
		ctx->spin_mesh.component_capacity = ctx->spin_mesh.component_count;
		ctx->spin_mesh.index_capacity = ctx->spin_mesh.index_count;
		ctx->spin_mesh.uses_normals = (ctx->WhichVectorMode == BOX1 || ctx->WhichVectorMode == ARROW1 || ctx->WhichVectorMode == CONE1) ? 1 : 0;
		unsigned int upload_mask = ctx->spin_mesh.uses_normals
			? VBO_UPLOAD_ALL
			: VBO_UPLOAD_VERTICES | VBO_UPLOAD_COLORS | VBO_UPLOAD_INDICES;
		UploadVBOMesh(&ctx->spin_mesh, ctx->vertices, ctx->normals, ctx->colors, ctx->indices, upload_mask);
		break;
	}
}

void ChangeColorMap ( magnoom_ctx *ctx, int id )
{
	switch (id)
	{
		case 0:
		ctx->HueMap[0]=ctx->HueMapRGB[0]+ctx->ColorShift;
		ctx->HueMap[1]=ctx->HueMapRGB[1]+ctx->ColorShift;
		ctx->HueMap[2]=ctx->HueMapRGB[2]+ctx->ColorShift;
		ctx->HueMap[3]=ctx->HueMapRGB[3]+ctx->ColorShift;
		ctx->HueMap[4]=ctx->HueMapRGB[4]+ctx->ColorShift;
		ctx->HueMap[5]=ctx->HueMapRGB[5]+ctx->ColorShift;
		// InitRGB(RHue, GHue, BHue, HueMap);
		case 1:
		ChangeVectorMode(ctx, 1);
	}
}

void ReallocateArrayDrawing(magnoom_ctx *ctx)
{
	free(ctx->vertexProto); free(ctx->normalProto); free(ctx->indicesProto);
	free(ctx->vertices); free(ctx->normals); free(ctx->colors); free(ctx->indices);

	int NOS_L=0;
	int NOB_L=0;
  
	switch( ctx->WhichSliceMode)  // see also UpdateIndices(ctx)
	{
		case A_AXIS:
		NOS_L=ctx->NOS_AL * (1+ctx->A_layer_max - ctx->A_layer_min);
		if (ctx->A_layer_max - ctx->A_layer_min<2){
			NOB_L=ctx->NOB_AL * (1+ctx->A_layer_max - ctx->A_layer_min);
		}else{
			NOB_L=2*ctx->NOB_AL+2*(ctx->A_layer_max - ctx->A_layer_min-1)*(ctx->uABC[1]-1+ctx->uABC[2]-1);
            if (ctx->uABC[2]==1 || ctx->uABC[1]==1){
                NOB_L=ctx->NOB_AL * (1+ctx->A_layer_max - ctx->A_layer_min);
            } 
		}
		break;
		case B_AXIS:
		NOS_L=ctx->NOS_BL * (1+ctx->B_layer_max - ctx->B_layer_min);
		if (ctx->B_layer_max - ctx->B_layer_min<2){
			NOB_L=ctx->NOB_BL * (1+ctx->B_layer_max - ctx->B_layer_min);
		}else{
			NOB_L=2*ctx->NOB_BL + 2*(ctx->B_layer_max - ctx->B_layer_min-1)*(ctx->uABC[0]-1+ctx->uABC[2]-1);
            if (ctx->uABC[0]==1 || ctx->uABC[2]==1){
                NOB_L=ctx->NOB_BL * (1+ctx->B_layer_max - ctx->B_layer_min);
            } 
		}
		break;
		case C_AXIS:
		NOS_L=ctx->NOS_CL * (1+ctx->C_layer_max - ctx->C_layer_min);
		if (ctx->C_layer_max - ctx->C_layer_min<2){
			NOB_L=ctx->NOB_CL * (1+ctx->C_layer_max - ctx->C_layer_min);
		}else{
			NOB_L=2*ctx->NOB_CL + 2*(ctx->C_layer_max - ctx->C_layer_min-1)*(ctx->uABC[0]-1+ctx->uABC[1]-1);
            if (ctx->uABC[0]==1 || ctx->uABC[1]==1){
                NOB_L=ctx->NOB_CL * (1+ctx->C_layer_max - ctx->C_layer_min);
            } 
		}
		break;
		case FILTER:
		NOS_L=ctx->NOS;
		NOB_L=ctx->NOS/ctx->AtomsPerBlock;
		break;
	}

	switch( ctx->WhichVectorMode )
	{
		case ARROW1:		
			// arrowFaces - number of arrow faces
			ctx->ElNumProto = 5*ctx->arrowFaces-4; // number of triangles per arrow
			ctx->IdNumProto = 3*ctx->ElNumProto; // number of indixes per arrow
			ctx->VCNumProto = 3*(2*(1+ctx->arrowFaces)-2+4*ctx->arrowFaces+3*ctx->arrowFaces);
			// Allocate memory for arrow prototype  
			ctx->vertexProto  = (float  *)malloc(ctx->VCNumProto * sizeof( float  ));
			ctx->normalProto  = (float  *)malloc(ctx->VCNumProto * sizeof( float  ));
			ctx->indicesProto = (GLuint *)malloc(ctx->IdNumProto * sizeof( GLuint ));	
			ctx->ElNum = NOS_L * ctx->ElNumProto;
			ctx->IdNum = NOS_L * ctx->IdNumProto;
			ctx->VCNum = NOS_L * ctx->VCNumProto;
			// Allocate memory for all ARROW1 (spins) 
			ctx->vertices	= (float  *)malloc(ctx->VCNum * sizeof( float  ));
			ctx->normals 	= (float  *)malloc(ctx->VCNum * sizeof( float  ));
			ctx->colors 		= (float  *)malloc(ctx->VCNum * sizeof( float  ));
			ctx->indices		= (GLuint *)malloc(ctx->IdNum * sizeof( GLuint ));	
			break;

		case CONE1:		
			// arrowFaces - number of cone faces
			ctx->ElNumProto = 2*ctx->arrowFaces-2; // number of triangles per cone
			ctx->IdNumProto = 3*ctx->ElNumProto; // number of indixes per cone
			ctx->VCNumProto = 3*((ctx->arrowFaces)+3*ctx->arrowFaces);// number of vertecies per cone			
			// Allocate memory for arrow prototype  
			ctx->vertexProto  = (float  *)malloc(ctx->VCNumProto * sizeof( float  ));
			ctx->normalProto  = (float  *)malloc(ctx->VCNumProto * sizeof( float  ));
			ctx->indicesProto = (GLuint *)malloc(ctx->IdNumProto * sizeof( GLuint ));	
			ctx->ElNum = NOS_L * ctx->ElNumProto;
			ctx->IdNum = NOS_L * ctx->IdNumProto;
			ctx->VCNum = NOS_L * ctx->VCNumProto;	
			// Allocate memory for all CONES (spins) 
			ctx->vertices	= (float  *)malloc(ctx->VCNum * sizeof( float  ));
			ctx->normals 	= (float  *)malloc(ctx->VCNum * sizeof( float  ));
			ctx->colors 		= (float  *)malloc(ctx->VCNum * sizeof( float  ));
			ctx->indices		= (GLuint *)malloc(ctx->IdNum * sizeof( GLuint ));	
			break;

		case BOX1:		
			ctx->ElNumProto = 6*2; // number of triangles per box
			ctx->IdNumProto = 3*ctx->ElNumProto; // number of indixes per box
			ctx->VCNumProto = 6*4*3; // number of vertecies per box			
			// Allocate memory for box prototype  
			ctx->vertexProto  = (float  *)malloc(ctx->VCNumProto* sizeof( float  ));
			ctx->normalProto  = (float  *)malloc(ctx->VCNumProto* sizeof( float  ));
			ctx->indicesProto = (GLuint *)malloc(ctx->IdNumProto* sizeof( GLuint ));	
			ctx->ElNum = NOB_L * ctx->ElNumProto;
			ctx->IdNum = NOB_L * ctx->IdNumProto;
			ctx->VCNum = NOB_L * ctx->VCNumProto;	
			// Allocate memory for all boxes (spins averaged over Block) 
			ctx->vertices	= (float  *)malloc(ctx->VCNum * sizeof( float  ));
			ctx->normals 	= (float  *)malloc(ctx->VCNum * sizeof( float  ));
			ctx->colors 		= (float  *)malloc(ctx->VCNum * sizeof( float  ));
			ctx->indices		= (GLuint *)malloc(ctx->IdNum * sizeof( GLuint ));	
			break;

		case uPOINT:		
			// arrowFaces - number of arrow faces
			ctx->IdNumProto = 1; // number of indixes per point
			ctx->VCNumProto = 3; // number component of ctx->vertices per point
			// Allocate memory for point prototype
			ctx->vertexProto  = (float  *)malloc(ctx->VCNumProto * sizeof( float  ));
			ctx->normalProto  = NULL;
			ctx->indicesProto = (GLuint *)malloc(ctx->IdNumProto * sizeof( GLuint ));	
			ctx->IdNum = NOS_L * ctx->IdNumProto;
			ctx->VCNum = NOS_L * ctx->VCNumProto; // total number of all components of all ctx->vertices 
			// Allocate memory for all points (spins) 
			ctx->vertices	= (float  *)malloc(ctx->VCNum * sizeof( float  ));
			ctx->normals 	= NULL;
			ctx->colors 		= (float  *)malloc(ctx->VCNum * sizeof( float  ));
			ctx->indices		= (GLuint *)malloc(ctx->IdNum * sizeof( GLuint ));	
			break;

		case CANE:	// is default vector mode

		default:
			ctx->IdNumProto = 2; // number of indixes per cane
			ctx->VCNumProto = 3*2; // number component of ctx->vertices per cane
			// Allocate memory for cane prototype
			ctx->vertexProto  = (float  *)malloc(ctx->VCNumProto * sizeof( float  ));
			ctx->normalProto  = NULL;
			ctx->indicesProto = (GLuint *)malloc(ctx->IdNumProto * sizeof( GLuint ));	
			ctx->IdNum = NOS_L * ctx->IdNumProto;
			ctx->VCNum = NOS_L * ctx->VCNumProto; // total number of all components of all ctx->vertices 
			// Allocate memory for all CANES (spins) 
			ctx->vertices	= (float  *)malloc(ctx->VCNum * sizeof( float  ));
			ctx->normals 	= NULL;
			ctx->colors 		= (float  *)malloc(ctx->VCNum * sizeof( float  ));
			ctx->indices		= (GLuint *)malloc(ctx->IdNum * sizeof( GLuint ));
	}
}

void ReallocateArrayDrawing_BextDC(magnoom_ctx *ctx)
{
	free(ctx->vertexProto_BextDC); free(ctx->normalProto_BextDC); free(ctx->indicesProto_BextDC);
	free(ctx->vertices_BextDC); free(ctx->normals_BextDC); free(ctx->colors_BextDC); free(ctx->indices_BextDC);
	// arrowFaces - number of arrow faces
	int ElNum_BextDC = 5*ctx->arrowFaces_BextDC-4; // number of triangles per arrow
	ctx->IdNum_BextDC = 3*ElNum_BextDC; // number of indixes per arrow
	ctx->VCNum_BextDC = 3*(2*(1+ctx->arrowFaces_BextDC)-2+4*ctx->arrowFaces_BextDC+3*ctx->arrowFaces_BextDC);
	// Allocate memory for the Bext DC arrow prototype
	ctx->vertexProto_BextDC  	= (float  *)malloc(ctx->VCNum_BextDC * sizeof( float  ));
	ctx->normalProto_BextDC  	= (float  *)malloc(ctx->VCNum_BextDC * sizeof( float  ));
	ctx->indicesProto_BextDC 	= (GLuint *)malloc(ctx->IdNum_BextDC * sizeof( GLuint ));
	// Allocate memory for the Bext DC arrow
	ctx->vertices_BextDC		= (float  *)malloc(ctx->VCNum_BextDC * sizeof( float  ));
	ctx->normals_BextDC 		= (float  *)malloc(ctx->VCNum_BextDC * sizeof( float  ));
	ctx->colors_BextDC 		= (float  *)malloc(ctx->VCNum_BextDC * sizeof( float  ));
	ctx->indices_BextDC		= (GLuint *)malloc(ctx->IdNum_BextDC * sizeof( GLuint ));
}

void ReallocateArrayDrawing_BOX(magnoom_ctx *ctx)
{
	free(ctx->vertices_BOX); free(ctx->normals_BOX); free(ctx->colors_BOX); free(ctx->indices_BOX);			
	int ElNum_BOX = 6*2*12; // number of triangles 
	ctx->IdNum_BOX = 3*ElNum_BOX; // number of indixes
	ctx->VCNum_BOX = 6*4*3*12;
	ctx->vertices_BOX	= (float  *)malloc(ctx->VCNum_BOX * sizeof( float  ));
	ctx->normals_BOX 	= (float  *)malloc(ctx->VCNum_BOX * sizeof( float  ));
	ctx->colors_BOX 		= (float  *)malloc(ctx->VCNum_BOX * sizeof( float  ));
	ctx->indices_BOX		= (GLuint *)malloc(ctx->IdNum_BOX * sizeof( GLuint ));				
}

void ReallocateArrayDrawing_BASIS(magnoom_ctx *ctx)
{
	free(ctx->vertices_BASIS); free(ctx->normals_BASIS); free(ctx->colors_BASIS); free(ctx->indices_BASIS);			
	int ElNum_BASIS = 6*2*4; // number of triangles 
	ctx->IdNum_BASIS = 3*ElNum_BASIS; // number of indixes
	ctx->VCNum_BASIS = 6*4*3*4;
	ctx->vertices_BASIS	= (float  *)malloc(ctx->VCNum_BASIS * sizeof( float  ));
	ctx->normals_BASIS	= (float  *)malloc(ctx->VCNum_BASIS * sizeof( float  ));
	ctx->colors_BASIS	= (float  *)malloc(ctx->VCNum_BASIS * sizeof( float  ));
	ctx->indices_BASIS	= (GLuint *)malloc(ctx->IdNum_BASIS * sizeof( GLuint ));				
}

void ReallocateArrayDrawing_PBC(magnoom_ctx *ctx)
{
	free(ctx->vertices_PBC_A); free(ctx->normals_PBC_A); free(ctx->colors_PBC_A); free(ctx->indices_PBC_A);		
	free(ctx->vertices_PBC_B); free(ctx->normals_PBC_B); free(ctx->colors_PBC_B); free(ctx->indices_PBC_B);		
	free(ctx->vertices_PBC_C); free(ctx->normals_PBC_C); free(ctx->colors_PBC_C); free(ctx->indices_PBC_C);			
	int ElNum_PBC = 6*2*16; // number of triangles 
	ctx->IdNum_PBC = 3*ElNum_PBC; // number of indixes 
	ctx->VCNum_PBC = 6*4*3*16;
	ctx->vertices_PBC_A	= (float  *)malloc(ctx->VCNum_PBC * sizeof( float  ));
	ctx->vertices_PBC_B	= (float  *)malloc(ctx->VCNum_PBC * sizeof( float  ));
	ctx->vertices_PBC_C	= (float  *)malloc(ctx->VCNum_PBC * sizeof( float  ));

	ctx->normals_PBC_A 	= (float  *)malloc(ctx->VCNum_PBC * sizeof( float  ));
	ctx->normals_PBC_B 	= (float  *)malloc(ctx->VCNum_PBC * sizeof( float  ));
	ctx->normals_PBC_C 	= (float  *)malloc(ctx->VCNum_PBC * sizeof( float  ));

	ctx->colors_PBC_A 	= (float  *)malloc(ctx->VCNum_PBC * sizeof( float  ));
	ctx->colors_PBC_B 	= (float  *)malloc(ctx->VCNum_PBC * sizeof( float  ));
	ctx->colors_PBC_C 	= (float  *)malloc(ctx->VCNum_PBC * sizeof( float  ));

	ctx->indices_PBC_A	= (GLuint *)malloc(ctx->VCNum_PBC * sizeof( GLuint ));
	ctx->indices_PBC_B	= (GLuint *)malloc(ctx->VCNum_PBC * sizeof( GLuint ));
	ctx->indices_PBC_C	= (GLuint *)malloc(ctx->VCNum_PBC * sizeof( GLuint ));					
}

static void FillProtoVerNorInd(magnoom_ctx *ctx, float * V, float * N, GLuint * I, int faces, int mode, int style)//faces = arrowFaces
{
	int   i, j;
	float H = 1.00f;	// total length
	float dF=D2R*360.f/(float)faces;
	float R = H*0.20f;	// big radius
	float r = H*0.06f;	// small radius 
	float h = H*0.40f;	// H - head
	float P = H*ctx->Pivot;	// pivot - central point on which the arrow turns (shift up in Z from 0 to 1)
	float cosF0, cosF1, cosF2, sinF0, sinF1, sinF2;
	float tmp0[3], tmp1[3], tmp2[3], tmp3[3];
	float v[8][3]={	{0,0,0},
					{ctx->abc[0][0],ctx->abc[0][1],ctx->abc[0][2]},
					{ctx->abc[1][0],ctx->abc[1][1],ctx->abc[1][2]},
					{ctx->abc[0][0]+ctx->abc[1][0],ctx->abc[0][1]+ctx->abc[1][1],ctx->abc[0][2]+ctx->abc[1][2]},
					{ctx->abc[2][0],ctx->abc[2][1],ctx->abc[2][2]},
					{ctx->abc[2][0]+ctx->abc[0][0],ctx->abc[2][1]+ctx->abc[0][1],ctx->abc[2][2]+ctx->abc[0][2]},
					{ctx->abc[2][0]+ctx->abc[1][0],ctx->abc[2][1]+ctx->abc[1][1],ctx->abc[2][2]+ctx->abc[1][2]},
					{ctx->abc[2][0]+ctx->abc[0][0]+ctx->abc[1][0],
						ctx->abc[2][1]+ctx->abc[0][1]+ctx->abc[1][1],
							ctx->abc[2][2]+ctx->abc[0][2]+ctx->abc[1][2]},
				};
	if (style==1)
	{
		R = H*0.1f;	// big radius
		r = H*0.03f;	// small radius 
		h = H*0.65f;	// H - head
	}

	// Arrow //////////////////////////////////// Cane//////////////////////////////////////////////
	//          Top____________                //                                                 //
	//     H    / \           ^                //             O v1(x1,y1,z1)                      // 
	//     e   /   \          |                //            /                                    //
	//     a  /     \         |                //           /                                     //
	//     d /__R*2__\____    H    Z           //          / GL_LINE                              //
	//       T|     |    ^    |    ^           //         /                                       //
	//	     a|     |    h    |    |           //        /                                        //
	//       i|     |    |    |    |           //       /                                         //
	//       l|_____|____v____v_   o --- > X   //      o v0(x0,y0,z0)                             //
	//          r*2                            //                                                 //
	////////////////////////////////////////////////////////////////////////////////////////////////

	i = -1; // vertex component index
	j = -1; // current index in array of vertex

	switch(mode)
	{
	case ARROW1:
		//head of arrow composed of N=faces triangles 
		for (int n=0; n<faces; n++) // n side index
		{
			cosF1 = R*cos(dF*n); cosF2 = R*cos(dF*(n+1)); 
			sinF1 = R*sin(dF*n); sinF2 = R*sin(dF*(n+1));
			tmp0[0] = R*R/sqrt((H-h)*(H-h) + R*R);

			// i++; V[i] = 0.f;   N[i] = 0;	
			// i++; V[i] = 0.f;   N[i] = 0;	
			// i++; V[i] = H-P;   N[i] = 1;	

			// i++; V[i] = cosF2; N[i] = cosF2;	
			// i++; V[i] = sinF2; N[i] = sinF2;	
			// i++; V[i] = h-P;   N[i] = tmp0[0];

			// i++; V[i] = cosF1; N[i] = cosF1;
			// i++; V[i] = sinF1; N[i] = sinF1;	 
			// i++; V[i] = h-P;   N[i] = tmp0[0];	

			i++; V[i] = tmp0[0] = 0.f;	
			i++; V[i] = tmp0[1] = 0.f;	
			i++; V[i] = tmp0[2] = H-P;	

			i++; V[i] = tmp2[0] = cosF2;	
			i++; V[i] = tmp2[1] = sinF2;
			i++; V[i] = tmp2[2] = h-P;

			i++; V[i] = tmp1[0] = cosF1; N[i] = cosF1;
			i++; V[i] = tmp1[1] = sinF1; N[i] = sinF1;	 
			i++; V[i] = tmp1[2] = h-P;	 N[i] = 0;	

     		Enorm( tmp0, tmp1, tmp2, tmp3);

			N[i-8] = N[i-5] = N[i-2] = tmp3[0] ; // nx
			N[i-7] = N[i-4] = N[i-1] = tmp3[1] ; // ny
			N[i-6] = N[i-3] = N[i-0] = tmp3[2] ; // nz

			j++; I[j] = n*3 + 0 ; 
			j++; I[j] = n*3 + 1 ;
			j++; I[j] = n*3 + 2 ;
		}

		//bottom side of head
		for (int n=0; n<faces; n++) // n face index
		{
			cosF1 = R*cos(dF*n); sinF1 = R*sin(dF*n);
			i++; V[i] = cosF1;	N[i] = 0.f;	// x_0, x_1, ... // nx_0, nx_1, ...
			i++; V[i] = sinF1;	N[i] = 0.f;	// y_0, y_1, ... // ny_0, ny_1, ...  		
			i++; V[i] = h - P;	N[i] =-1.f;	// z_0, z_1, ... // nz_0, nz_1, ...   		
		}

		for (int n=0; n<faces-2; n++) // n face index
		{
			j++; I[j] = 3*faces;
			j++; I[j] = (n+1)%(faces) 	+ 3*faces;
			j++; I[j] = (n+2)%(faces) 	+ 3*faces; 
		}

		//tail of arrow composed of rectangles:
		// v2 o--o v4
		//    |\ |
		//    | \|
		// v1 o--o v3 

		for (int n=0; n<faces; n++) // n face index
		{
			cosF0 = 1*cos(dF*n+dF/2.f); sinF0 = 1*sin(dF*n+dF/2.f);	// for ctx->normals
			//cosF0 = 1*cos(dF*n); sinF0 = 1*sin(dF*n);	// for ctx->normals
			cosF1 = r*cos(dF*n); 		sinF1 = r*sin(dF*n);		// for v1,v2
			cosF2 = r*cos(dF*(n+1)); 	sinF2 = r*sin(dF*(n+1));	// for v3,v4
			//v1
			i++; V[i] = cosF1;	N[i] = cosF0;
			i++; V[i] = sinF1;	N[i] = sinF0;
			i++; V[i] = 0.f-P; 	N[i] = 0.f;
			//v2
			i++; V[i] = cosF1;	N[i] = cosF0;
			i++; V[i] = sinF1;	N[i] = sinF0;
			i++; V[i] = h - P;	N[i] = 0.f;
			//v3
			i++; V[i] = cosF2;	N[i] = cosF0;
			i++; V[i] = sinF2;	N[i] = sinF0;
			i++; V[i] = 0.f-P; 	N[i] = 0.f;
			//v4
			i++; V[i] = cosF2;	N[i] = cosF0;
			i++; V[i] = sinF2;	N[i] = sinF0;
			i++; V[i] = h-P; 	N[i] = 0.f;
			// //v1
			// i++; V[i] = cosF1;	N[i] = cosF1;
			// i++; V[i] = sinF1;	N[i] = sinF1;
			// i++; V[i] = 0.f-P; 	N[i] = 0.f;
			// //v2
			// i++; V[i] = cosF1;	N[i] = cosF1;
			// i++; V[i] = sinF1;	N[i] = sinF1;
			// i++; V[i] = h - P;	N[i] = 0.f;
			// //v3
			// i++; V[i] = cosF2;	N[i] = cosF2;
			// i++; V[i] = sinF2;	N[i] = sinF2;
			// i++; V[i] = 0.f-P; 	N[i] = 0.f;
			// //v4
			// i++; V[i] = cosF2;	N[i] = cosF2;
			// i++; V[i] = sinF2;	N[i] = sinF2;
			// i++; V[i] = h-P; 	N[i] = 0.f;
		}

		for (int n=0; n<faces; n++) // n face index+faces+1
		{
			j++; I[j] = (4*n+0) + 4*faces; //v1
			j++; I[j] = (4*n+1) + 4*faces; //v2
			j++; I[j] = (4*n+2) + 4*faces; //v3

			j++; I[j] = (4*n+2) + 4*faces; //v3		
			j++; I[j] = (4*n+1) + 4*faces; //v2
			j++; I[j] = (4*n+3) + 4*faces; //v4
		}
	
		//bottom side of tail
		for (int n=0; n<faces; n++) // n face index
		{
			cosF1 = r*cos(dF*n); sinF1 = r*sin(dF*n);
			i++; V[i] = cosF1;	N[i] = 0.f;	// x_0, x_1, ... // nx_0, nx_1, ...
			i++; V[i] = sinF1;	N[i] = 0.f;	// y_0, y_1, ... // ny_0, ny_1, ...  		
			i++; V[i] = 0.f-P; 	N[i] =-1.f;	// z_0, z_1, ... // nz_0, nz_1, ...  		
		}

		for (int n=0; n<faces-2; n++) // n face index
		{ 
			j++; I[j] = 8*faces;
			j++; I[j] = (n+1)%(faces) 	+ 8*faces;
			j++; I[j] = (n+2)%(faces) 	+ 8*faces;
		}
			
	break;
				
	case CONE1:
		//head of arrow composed of N=faces triangles 
		for (int n=0; n<faces; n++) // n side index
		{
			cosF1 = R*cos(dF*n); cosF2 = R*cos(dF*(n+1)); 
			sinF1 = R*sin(dF*n); sinF2 = R*sin(dF*(n+1));
			tmp0[0] = R*R/sqrt(H*H + R*R);

			i++; V[i] = 0.f;   N[i] = 0;	
			i++; V[i] = 0.f;   N[i] = 0;	
			i++; V[i] = H-P;   N[i] = 1;	

			i++; V[i] = cosF2; N[i] = cosF2;	
			i++; V[i] = sinF2; N[i] = sinF2;	
			i++; V[i] = -P;   N[i] = tmp0[0];

			i++; V[i] = cosF1; N[i] = cosF1;
			i++; V[i] = sinF1; N[i] = sinF1;	 
			i++; V[i] = -P;   N[i] = tmp0[0];

			// i++; V[i] = tmp0[0] = 0.f;		
			// i++; V[i] = tmp0[1] = 0.f;		
			// i++; V[i] = tmp0[2] = H-P;		

			// i++; V[i] = tmp2[0] = cosF2;	
			// i++; V[i] = tmp2[1] = sinF2;	
			// i++; V[i] = tmp2[2] = - P;	

			// i++; V[i] = tmp1[0] = cosF1;	
			// i++; V[i] = tmp1[1] = sinF1;	 
			// i++; V[i] = tmp1[2] = - P;		

			// Enorm( tmp0, tmp1, tmp2, tmp3);

			// N[i-8] = N[i-5] = N[i-2] = tmp3[0] ; // nx
			// N[i-7] = N[i-4] = N[i-1] = tmp3[1] ; // ny
			// N[i-6] = N[i-3] = N[i-0] = tmp3[2] ; // nz

            N[i-8] = N[i-5] = N[i-2]; // nx
            N[i-7] = N[i-4] = N[i-1]; // ny
            N[i-6] = N[i-3] = N[i-0]; // nz

			j++; I[j] = n*3 + 0 ; 
			j++; I[j] = n*3 + 1 ;
			j++; I[j] = n*3 + 2 ;
		}

		//bottom side of head
		for (int n=0; n<faces; n++) // n face index
		{
			cosF1 = R*cos(dF*n); sinF1 = R*sin(dF*n);
			i++; V[i] = cosF1;	N[i] = 0.f;	// x_0, x_1, ... // nx_0, nx_1, ...
			i++; V[i] = sinF1;	N[i] = 0.f;	// y_0, y_1, ... // ny_0, ny_1, ...  		
			i++; V[i] = - P;	N[i] =-1.f;	// z_0, z_1, ... // nz_0, nz_1, ...   		
		}

		for (int n=0; n<faces-2; n++) // n-2=num of triangles on the bottom side
		{ 
			j++; I[j] = 3*faces;
			j++; I[j] = (n+1)%(faces) 	+ 3*faces;
			j++; I[j] = (n+2)%(faces) 	+ 3*faces;
		}
	break;

	case BOX1:

			Enorm2( v[0], v[1], v[2], tmp1);
			Enorm2( v[0], v[4], v[1], tmp2);
			Enorm2( v[0], v[2], v[4], tmp3);
			float handedness = tmp1[0]*v[4][0] + tmp1[1]*v[4][1] + tmp1[2]*v[4][2] < 0.0f ? -1.0f : 1.0f;
			for (int component = 0; component < 3; ++component) {
				tmp1[component] *= handedness;
				tmp2[component] *= handedness;
				tmp3[component] *= handedness;
			}

			V[++i] = v[0][0]; N[i] = -tmp1[0]; V[i+12] = v[4][0]; N[i+12] = tmp1[0];
			V[++i] = v[0][0]; N[i] = -tmp1[1]; V[i+12] = v[4][1]; N[i+12] = tmp1[1];
			V[++i] = v[0][0]; N[i] = -tmp1[2]; V[i+12] = v[4][2]; N[i+12] = tmp1[2];

			V[++i] = v[1][0]; N[i] = -tmp1[0]; V[i+12] = v[6][0]; N[i+12] = tmp1[0];
			V[++i] = v[1][1]; N[i] = -tmp1[1]; V[i+12] = v[6][1]; N[i+12] = tmp1[1];
			V[++i] = v[1][2]; N[i] = -tmp1[2]; V[i+12] = v[6][2]; N[i+12] = tmp1[2];

			V[++i] = v[2][0]; N[i] = -tmp1[0]; V[i+12] = v[5][0]; N[i+12] = tmp1[0];
			V[++i] = v[2][1]; N[i] = -tmp1[1]; V[i+12] = v[5][1]; N[i+12] = tmp1[1];
			V[++i] = v[2][2]; N[i] = -tmp1[2]; V[i+12] = v[5][2]; N[i+12] = tmp1[2];

			V[++i] = v[3][0]; N[i] = -tmp1[0]; V[i+12] = v[7][0]; N[i+12] = tmp1[0];
			V[++i] = v[3][1]; N[i] = -tmp1[1]; V[i+12] = v[7][1]; N[i+12] = tmp1[1];
			V[++i] = v[3][2]; N[i] = -tmp1[2]; V[i+12] = v[7][2]; N[i+12] = tmp1[2];
            i+=12;
			I[++j] = 0; I[++j] = 1; I[++j] = 2; // first  triangle bottom
			I[++j] = 2; I[++j] = 1; I[++j] = 3; // second triangle bottom
			I[++j] = 4; I[++j] = 5; I[++j] = 6; // third  triangle top
			I[++j] = 6; I[++j] = 5; I[++j] = 7; // fourth triangle top

            // right/left
			V[++i] = v[0][0]; N[i] = -tmp2[0]; V[i+12] = v[3][0]; N[i+12] = tmp2[0];
			V[++i] = v[0][0]; N[i] = -tmp2[1]; V[i+12] = v[3][1]; N[i+12] = tmp2[1];
			V[++i] = v[0][0]; N[i] = -tmp2[2]; V[i+12] = v[3][2]; N[i+12] = tmp2[2];

			V[++i] = v[4][0]; N[i] = -tmp2[0]; V[i+12] = v[7][0]; N[i+12] = tmp2[0];
			V[++i] = v[4][1]; N[i] = -tmp2[1]; V[i+12] = v[7][1]; N[i+12] = tmp2[1];
			V[++i] = v[4][2]; N[i] = -tmp2[2]; V[i+12] = v[7][2]; N[i+12] = tmp2[2];

			V[++i] = v[1][0]; N[i] = -tmp2[0]; V[i+12] = v[2][0]; N[i+12] = tmp2[0];
			V[++i] = v[1][1]; N[i] = -tmp2[1]; V[i+12] = v[2][1]; N[i+12] = tmp2[1];
			V[++i] = v[1][2]; N[i] = -tmp2[2]; V[i+12] = v[2][2]; N[i+12] = tmp2[2];

			V[++i] = v[5][0]; N[i] = -tmp2[0]; V[i+12] = v[6][0]; N[i+12] = tmp2[0];
			V[++i] = v[5][1]; N[i] = -tmp2[1]; V[i+12] = v[6][1]; N[i+12] = tmp2[1];
			V[++i] = v[5][2]; N[i] = -tmp2[2]; V[i+12] = v[6][2]; N[i+12] = tmp2[2];
            i+=12;
			I[++j] = 8+0; I[++j] = 8+1; I[++j] = 8+2; // first  triangle 
			I[++j] = 8+2; I[++j] = 8+1; I[++j] = 8+3; // second triangle 
			I[++j] = 8+4; I[++j] = 8+5; I[++j] = 8+6; // third  triangle 
			I[++j] = 8+6; I[++j] = 8+5; I[++j] = 8+7; // fourth triangle 

            // front/back
			V[++i] = v[0][0]; N[i] = -tmp3[0]; V[i+12] = v[1][0]; N[i+12] = tmp3[0];
			V[++i] = v[0][0]; N[i] = -tmp3[1]; V[i+12] = v[1][1]; N[i+12] = tmp3[1];
			V[++i] = v[0][0]; N[i] = -tmp3[2]; V[i+12] = v[1][2]; N[i+12] = tmp3[2];

			V[++i] = v[2][0]; N[i] = -tmp3[0]; V[i+12] = v[5][0]; N[i+12] = tmp3[0];
			V[++i] = v[2][1]; N[i] = -tmp3[1]; V[i+12] = v[5][1]; N[i+12] = tmp3[1];
			V[++i] = v[2][2]; N[i] = -tmp3[2]; V[i+12] = v[5][2]; N[i+12] = tmp3[2];

			V[++i] = v[4][0]; N[i] = -tmp3[0]; V[i+12] = v[3][0]; N[i+12] = tmp3[0];
			V[++i] = v[4][1]; N[i] = -tmp3[1]; V[i+12] = v[3][1]; N[i+12] = tmp3[1];
			V[++i] = v[4][2]; N[i] = -tmp3[2]; V[i+12] = v[3][2]; N[i+12] = tmp3[2];

			V[++i] = v[6][0]; N[i] = -tmp3[0]; V[i+12] = v[7][0]; N[i+12] = tmp3[0];
			V[++i] = v[6][1]; N[i] = -tmp3[1]; V[i+12] = v[7][1]; N[i+12] = tmp3[1];
			V[++i] = v[6][2]; N[i] = -tmp3[2]; V[i+12] = v[7][2]; N[i+12] = tmp3[2];

			I[++j] = 2*8+0; I[++j] = 2*8+1; I[++j] = 2*8+2; // first  triangle 
			I[++j] = 2*8+2; I[++j] = 2*8+1; I[++j] = 2*8+3; // second triangle 
			I[++j] = 2*8+4; I[++j] = 2*8+5; I[++j] = 2*8+6; // third  triangle 
			I[++j] = 2*8+6; I[++j] = 2*8+5; I[++j] = 2*8+7; // fourth triangle 
			if (handedness < 0.0f) {
				for (int triangle = 0; triangle < 12; ++triangle) {
					GLuint swap = I[3*triangle + 1];
					I[3*triangle + 1] = I[3*triangle + 2];
					I[3*triangle + 2] = swap;
				}
			}
	break;

	case uPOINT:
			i++; V[i] = 0.f;		
			i++; V[i] = 0.f;		
			i++; V[i] = 0.f;
			j++; I[j] = j; // 0
	break;

	case CANE:	// canes is default visual mode 
	default:
			i++; V[i] = 0.f;		
			i++; V[i] = 0.f;		
			i++; V[i] = H-P;
			j++; I[j] = j; // 0

			i++; V[i] = 0.f;		
			i++; V[i] = 0.f;		
			i++; V[i] = -P;
			j++; I[j] = j; // 1	
	}	
}

void UpdateProtoVerNorInd_Spins(magnoom_ctx *ctx)
{
	FillProtoVerNorInd(ctx, ctx->vertexProto, ctx->normalProto, ctx->indicesProto,
		ctx->arrowFaces, ctx->WhichVectorMode, 0);
}

void UpdateProtoVerNorInd_BextDC(magnoom_ctx *ctx)
{
	FillProtoVerNorInd(ctx, ctx->vertexProto_BextDC, ctx->normalProto_BextDC, ctx->indices_BextDC,
		ctx->arrowFaces_BextDC, ARROW1, 1);
}

void UpdateIndices(magnoom_ctx *ctx)
{	// VerN number of vertesiec components per one prototype arrow
	const GLuint *Iinp = ctx->indicesProto;
	const int Kinp = ctx->IdNumProto;
	GLuint *Iout = ctx->indices;
	const int VerN = ctx->VCNumProto;
	int i,j;
	int NOS_L=0;
	int NOB_L=0;
    // NumberElementInSlice(NOS_L, NOB_L);
    switch( ctx->WhichSliceMode)
    {
        case A_AXIS:
        NOS_L=ctx->NOS_AL * (1+ctx->A_layer_max - ctx->A_layer_min);
        if (ctx->A_layer_max - ctx->A_layer_min<2){
            NOB_L=ctx->NOB_AL * (1+ctx->A_layer_max - ctx->A_layer_min);
        }else{
            NOB_L=2*ctx->NOB_AL+2*(ctx->A_layer_max - ctx->A_layer_min-1)*(ctx->uABC[1]-1+ctx->uABC[2]-1);
            if (ctx->uABC[2]==1 || ctx->uABC[1]==1){
                NOB_L=ctx->NOB_AL * (1+ctx->A_layer_max - ctx->A_layer_min);
            } 
        }
        break;
        case B_AXIS:
        NOS_L=ctx->NOS_BL * (1+ctx->B_layer_max - ctx->B_layer_min);
        if (ctx->B_layer_max - ctx->B_layer_min<2){
            NOB_L=ctx->NOB_BL * (1+ctx->B_layer_max - ctx->B_layer_min);
        }else{
            NOB_L=2*ctx->NOB_BL + 2*(ctx->B_layer_max - ctx->B_layer_min-1)*(ctx->uABC[0]-1+ctx->uABC[2]-1);
            if (ctx->uABC[0]==1 || ctx->uABC[2]==1){
                NOB_L=ctx->NOB_BL * (1+ctx->B_layer_max - ctx->B_layer_min);
            } 
        }
        break;
        case C_AXIS:
        NOS_L=ctx->NOS_CL * (1+ctx->C_layer_max - ctx->C_layer_min);
        if (ctx->C_layer_max - ctx->C_layer_min<2){
            NOB_L=ctx->NOB_CL * (1+ctx->C_layer_max - ctx->C_layer_min);
        }else{
            NOB_L=2*ctx->NOB_CL + 2*(ctx->C_layer_max - ctx->C_layer_min-1)*(ctx->uABC[0]-1+ctx->uABC[1]-1);
            if (ctx->uABC[0]==1 || ctx->uABC[1]==1){
                NOB_L=ctx->NOB_CL * (1+ctx->C_layer_max - ctx->C_layer_min);
            } 
        }
        break;
        case FILTER:
        NOS_L=ctx->NOS;
        NOB_L=ctx->NOS/ctx->AtomsPerBlock;
        break;
    }

	if (ctx->WhichVectorMode==BOX1){
		for (int n = 0; n<NOB_L; n++)//NOB_L number of blocks per layer(s)
		{
			i = n*(VerN/3); // shift in vertex array 
			j = n*Kinp; // shift in index array (for arrow and cone vertex num<index num!)
			for (int k=0; k<Kinp; k++) Iout[j+k] = i+Iinp[k];		
		}
	}else{
		for (int n = 0; n<NOS_L; n++)//NOS_L number of spins per layer(s)
		{
			i = n*(VerN/3); // shift in vertex array 
			j = n*Kinp; // shift in index array (for arrow and cone vertex num<index num!)
			for (int k=0; k<Kinp; k++) Iout[j+k] = i+Iinp[k];		
		}		
	}		
}



void UpdateVerticesNormalsColors(magnoom_ctx *ctx)
{
	const float *Vinp = ctx->vertexProto;
	const float *Ninp = ctx->normalProto;
	const int Kinp = ctx->VCNumProto;
	float *Vout = ctx->vertices;
	float *Nout = ctx->normals;
	float *Cout = ctx->colors;
	const float *Px = ctx->PosX;
	const float *Py = ctx->PosY;
	const float *Pz = ctx->PosZ;
	const double *spinArr = ctx->bS;
	//float tmpV1[3], tmpV2[3], tmpV3[3], RGB[3], U, A;
	float S[3], RGB[3], U, A;
	float vlength;
	int i,j,n,N;
	int anini=0;
	int anfin=ctx->uABC[0];
	int bnini=0;
	int bnfin=ctx->uABC[1];
	int cnini=0;
	int cnfin=ctx->uABC[2];
	bool F;//factor

	switch( ctx->WhichSliceMode)//see also ReallocateArrayDrawing(ctx)
	{
		case A_AXIS:
			anini=ctx->A_layer_min-1;
	        anfin=ctx->A_layer_max;
		break;
		case B_AXIS:
			bnini=ctx->B_layer_min-1;
	        bnfin=ctx->B_layer_max;
		break;
		case C_AXIS:
			cnini=ctx->C_layer_min-1;
	        cnfin=ctx->C_layer_max;
		break;
		default:
		break;
	}
	j=-1;
	switch (ctx->WhichVectorMode)
	{
	case ARROW1:
	case CONE1:
		if (ctx->WhichSliceMode==FILTER)
			{
			ctx->N_filter=0;
			for (int an = 0; an<ctx->uABC[0]; an++) {
			for (int bn = 0; bn<ctx->uABC[1]; bn++) {
			for (int cn = 0; cn<ctx->uABC[2]; cn++) {
				n = an+bn*ctx->uABC[0]+cn*ctx->uABC[0]*ctx->uABC[1];// index of the block
				n = n*ctx->AtomsPerBlock;//index of the first spin in the block
				if (ctx->GreedFilter){
					F=(an>=ctx->GreedFilterMinA && an<=ctx->GreedFilterMaxA) &&
						  (bn>=ctx->GreedFilterMinB && bn<=ctx->GreedFilterMaxB) &&
						  (cn>=ctx->GreedFilterMinC && cn<=ctx->GreedFilterMaxC) ;
					if (ctx->GreedFilterInvert) F=!F;					
				}else{F=true;}


				if (F)
				{
					for (int atom=0; atom<ctx->AtomsPerBlock; atom++)//atom is the index of the atom in block
					{
					    N = n + atom;
						//slow version is commented but easy to read: 
						S[0] = VEC_X(spinArr,N);
						S[1] = VEC_Y(spinArr,N);
						S[2] = VEC_Z(spinArr,N);
						int phi = atan2int( S[1], S[0] );// return integer angle phi 0 - 360
						// if (SpinFilter1){
						// 	if (!PhiInvert1 && !PhiInvert2){
						// 		F=((S[2]>=Sz_min1 && S[2]<=Sz_max1) &&  (phi>=phi_min1 && phi<=phi_max1)) ||
						// 	      ((S[2]>=Sz_min2 && S[2]<=Sz_max2) &&  (phi>=phi_min2 && phi<=phi_max2))  ;	
						// 	}else if (!PhiInvert1 && PhiInvert2){
						// 		F=((S[2]>=Sz_min1 && S[2]<=Sz_max1) &&  (phi>=phi_min1 && phi<=phi_max1)) ||
						// 	      ((S[2]>=Sz_min2 && S[2]<=Sz_max2) && !(phi>=phi_min2 && phi<=phi_max2))  ;	
						// 	}else if ( PhiInvert1 && PhiInvert2){
						// 		F=((S[2]>=Sz_min1 && S[2]<=Sz_max1) && !(phi>=phi_min1 && phi<=phi_max1)) ||
						// 	      ((S[2]>=Sz_min2 && S[2]<=Sz_max2) && !(phi>=phi_min2 && phi<=phi_max2))  ;									
						// 	}				
						// }else{F=true;}

						// if ( (S[2]>=Sz_min1 && S[2]<=Sz_max1) && (phi>=phi_min1 && phi<=phi_max1) ) 
						// {
						// 	F=true;}else{
						// 		F=false;
						// 	}
						F =	(ctx->SpinFilter1 && ((S[2]>=ctx->Sz_min1 && S[2]<=ctx->Sz_max1) && (phi>=ctx->phi_min1 && phi<=ctx->phi_max1)) ) ||
							(ctx->SpinFilter2 && ((S[2]>=ctx->Sz_min2 && S[2]<=ctx->Sz_max2) && (phi>=ctx->phi_min2 && phi<=ctx->phi_max2)) ) ||
							(ctx->SpinFilter3 && ((S[2]>=ctx->Sz_min3 && S[2]<=ctx->Sz_max3) && (phi>=ctx->phi_min3 && phi<=ctx->phi_max3)) ) ;													

						if (F)
						{
							HSVtoRGB(ctx, S, RGB, ctx->InvertValue, ctx->InvertHue);//metka
                            vlength = Unitf(S,S);
							ctx->N_filter++;
					        // HSVtoRGB(S, RGB, InvertValue, InvertHue);
							j++;
							for (int k=0; k<Kinp/3; k++) // k runs over ctx->vertices of the arrow/cone 
							{
								i = j*Kinp + 3*k;	// vertex index
								//slow version is commented but easy to read: 
								// tmpV1[0] = Vinp[3*k+0];
								// tmpV1[1] = Vinp[3*k+1];
								// tmpV1[2] = Vinp[3*k+2];
								// NewBasisCartesian(tmpV1, S, tmpV3); // to find arrow vector components w.r.t basis based on Sx,Sy,Sz
								// Vout[i+0] = tmpV3[0]*Scale + Px[n];	// new x-component of vertex + translation
								// Vout[i+1] = tmpV3[1]*Scale + Py[n];	// new y-component of vertex + translation
								// Vout[i+2] = tmpV3[2]*Scale + Pz[n];	// new z-component of vertex + translation
								U = 1.0f/(S[0]*S[0]+S[1]*S[1]+(1e-37f)); 
								
								A = (-S[1]*Vinp[3*k+0] + S[0]*Vinp[3*k+1])*(1. - S[2])*U; 

								Vout[i+0] =ctx->Kind[N]*( (-S[1]*A + Vinp[3*k+0]*S[2] + S[0]*Vinp[3*k+2]			  )*ctx->Scale*vlength + Px[N]);
								Vout[i+1] =ctx->Kind[N]*( ( S[0]*A + Vinp[3*k+1]*S[2] + S[1]*Vinp[3*k+2]			  )*ctx->Scale*vlength + Py[N]);
								Vout[i+2] =ctx->Kind[N]*( ( Vinp[3*k+2]*S[2] - (S[0]*Vinp[3*k+0]+S[1]*Vinp[3*k+1]) )*ctx->Scale*vlength + Pz[N]);	

								//slow version is commented but easy to read:
								// tmpV1[0] = Ninp[3*k+0];
								// tmpV1[1] = Ninp[3*k+1];
								// tmpV1[2] = Ninp[3*k+2];
								//NewBasisCartesian(tmpV1, S, tmpV3);	// to find ctx->vertices ctx->normals w.r.t basis based on Sx,Sy,Sz
								// Nout[i+0] = tmpV3[0];		// x-component of vertex normal
								// Nout[i+1] = tmpV3[1];		// y-component of vertex normal
								// Nout[i+2] = tmpV3[2];		// z-component of vertex normal

								A = (-S[1]*Ninp[3*k+0] + S[0]*Ninp[3*k+1])*(1. - S[2])*U;

								Nout[i+0] =-S[1]*A + Ninp[3*k+0]*S[2] + S[0]*Ninp[3*k+2];
								Nout[i+1] = S[0]*A + Ninp[3*k+1]*S[2] + S[1]*Ninp[3*k+2];
								Nout[i+2] = Ninp[3*k+2]*S[2] - (S[0]*Ninp[3*k+0]+S[1]*Ninp[3*k+1]);

								Cout[i+0] = RGB[0];			// x-component of vertex normal
								Cout[i+1] = RGB[1];			// y-component of vertex normal
								Cout[i+2] = RGB[2];			// z-component of vertex normal
							}
						}//if
					}
				}
			}
			}
			}
		}else{//IF SLICING MODE = A-AXIS, B-AXIS, C-AXIS 
			for (int an = anini; an<anfin; an++) 
			{
			for (int bn = bnini; bn<bnfin; bn++) 
			{
			for (int cn = cnini; cn<cnfin; cn++) 
			{
				n = an+bn*ctx->uABC[0]+cn*ctx->uABC[0]*ctx->uABC[1];// index of the block
				n = n*ctx->AtomsPerBlock;//index of the first spin in the block
				for (int atom=0; atom<ctx->AtomsPerBlock; atom++)//atom is the index of the atom in block
				{
				    N = n + atom;
					//slow version is commented but easy to read: 
					S[0] = VEC_X(spinArr,N);
					S[1] = VEC_Y(spinArr,N);
					S[2] = VEC_Z(spinArr,N);
                    HSVtoRGB(ctx, S, RGB, ctx->InvertValue, ctx->InvertHue);//metka
					vlength = Unitf(S,S);
			        // HSVtoRGB( S, RGB, InvertValue, InvertHue);
					j++;
					if (S[2]==-1){
					for (int k=0; k<Kinp/3; k++){// k runs over ctx->vertices 
							i = j*Kinp + 3*k;
							Vout[i+0] = ctx->Kind[N]*((-Vinp[3*k+0])*ctx->Scale*vlength + Px[N]);
							Vout[i+1] = ctx->Kind[N]*(( Vinp[3*k+1])*ctx->Scale*vlength + Py[N]);
							Vout[i+2] = ctx->Kind[N]*((-Vinp[3*k+2])*ctx->Scale*vlength + Pz[N]);	

							Nout[i+0] = -Ninp[3*k+0];
							Nout[i+1] =  Ninp[3*k+1];
							Nout[i+2] = -Ninp[3*k+2];

							Cout[i+0] = RGB[0];			// x-component of vertex normal
							Cout[i+1] = RGB[1];			// y-component of vertex normal
							Cout[i+2] = RGB[2];			// z-component of vertex normal
						}
					}else{
					for (int k=0; k<Kinp/3; k++) // k runs over ctx->vertices of the arrow/cone 
						{
							i = j*Kinp + 3*k;	// vertex index
							//slow version is commented but easy to read: 
							// tmpV1[0] = Vinp[3*k+0];
							// tmpV1[1] = Vinp[3*k+1];
							// tmpV1[2] = Vinp[3*k+2];
							// NewBasisCartesian(tmpV1, S, tmpV3); // to find arrow vector components w.r.t basis based on Sx,Sy,Sz
							// Vout[i+0] = tmpV3[0]*Scale + Px[n];	// new x-component of vertex + translation
							// Vout[i+1] = tmpV3[1]*Scale + Py[n];	// new y-component of vertex + translation
							// Vout[i+2] = tmpV3[2]*Scale + Pz[n];	// new z-component of vertex + translation
							U = 1.0f/(S[0]*S[0]+S[1]*S[1]+(1e-37f)); 
							
							A = (-S[1]*Vinp[3*k+0] + S[0]*Vinp[3*k+1])*(1. - S[2])*U; 

							Vout[i+0] =ctx->Kind[N]*((-S[1]*A + Vinp[3*k+0]*S[2] + S[0]*Vinp[3*k+2]			 )*ctx->Scale*vlength + Px[N]);
							Vout[i+1] =ctx->Kind[N]*(( S[0]*A + Vinp[3*k+1]*S[2] + S[1]*Vinp[3*k+2]			 )*ctx->Scale*vlength + Py[N]);
							Vout[i+2] =ctx->Kind[N]*(( Vinp[3*k+2]*S[2] - (S[0]*Vinp[3*k+0]+S[1]*Vinp[3*k+1]) )*ctx->Scale*vlength + Pz[N]);

							//slow version is commented but easy to read:
							// tmpV1[0] = Ninp[3*k+0];
							// tmpV1[1] = Ninp[3*k+1];
							// tmpV1[2] = Ninp[3*k+2];
							//NewBasisCartesian(tmpV1, S, tmpV3);	// to find ctx->vertices ctx->normals w.r.t basis based on Sx,Sy,Sz
							// Nout[i+0] = tmpV3[0];		// x-component of vertex normal
							// Nout[i+1] = tmpV3[1];		// y-component of vertex normal
							// Nout[i+2] = tmpV3[2];		// z-component of vertex normal

							A = (-S[1]*Ninp[3*k+0] + S[0]*Ninp[3*k+1])*(1. - S[2])*U;

							Nout[i+0] =-S[1]*A + Ninp[3*k+0]*S[2] + S[0]*Ninp[3*k+2];
							Nout[i+1] = S[0]*A + Ninp[3*k+1]*S[2] + S[1]*Ninp[3*k+2];
							Nout[i+2] = Ninp[3*k+2]*S[2] - (S[0]*Ninp[3*k+0]+S[1]*Ninp[3*k+1]);

							Cout[i+0] = RGB[0];			// x-component of vertex normal
							Cout[i+1] = RGB[1];			// y-component of vertex normal
							Cout[i+2] = RGB[2];			// z-component of vertex normal
						}
					}
				}
			}
			}
			}
		}//if
	break;

	case BOX1:
		if (ctx->WhichSliceMode==FILTER)
		{
			ctx->N_filter=0;
			for (int an = 0; an<ctx->uABC[0]; an++) {
			for (int bn = 0; bn<ctx->uABC[1]; bn++) {
			for (int cn = 0; cn<ctx->uABC[2]; cn++) {	
				if (ctx->GreedFilter){
						F =	(an>=ctx->GreedFilterMinA && an<=ctx->GreedFilterMaxA) &&
							(bn>=ctx->GreedFilterMinB && bn<=ctx->GreedFilterMaxB) &&
							(cn>=ctx->GreedFilterMinC && cn<=ctx->GreedFilterMaxC) ;
						if (ctx->GreedFilterInvert) F=!F;					
					}else{F=true;}

				if (F)
				{
					S[0] = 0;
					S[1] = 0;
					S[2] = 0;
					n = an+bn*ctx->uABC[0]+cn*ctx->uABC[0]*ctx->uABC[1];// index of the block
					for (int atom=0; atom<ctx->AtomsPerBlock; atom++)//atom is the index of the atom in block
					{
					    N = atom + n*ctx->AtomsPerBlock;//n*AtomsPerBlock=index of the first spin in the block;
		 
						S[0]+= VEC_X(spinArr,N);
						S[1]+= VEC_Y(spinArr,N);
						S[2]+= VEC_Z(spinArr,N);
				    }
					int phi = atan2int( S[1], S[0] );// return integer angle phi 0 - 360
					//if (S[2]>=Sz_min1 && S[2]<=Sz_max1)
					//if ( (S[2]>=Sz_min1 && S[2]<=Sz_max1) && (phi>=phi_min1 && phi<=phi_max1) )
					F =	(ctx->SpinFilter1 && ((S[2]>=ctx->Sz_min1 && S[2]<=ctx->Sz_max1) && (phi>=ctx->phi_min1 && phi<=ctx->phi_max1)) ) ||
						(ctx->SpinFilter2 && ((S[2]>=ctx->Sz_min2 && S[2]<=ctx->Sz_max2) && (phi>=ctx->phi_min2 && phi<=ctx->phi_max2)) ) ||
						(ctx->SpinFilter3 && ((S[2]>=ctx->Sz_min3 && S[2]<=ctx->Sz_max3) && (phi>=ctx->phi_min3 && phi<=ctx->phi_max3)) ) ;

					if (F)
					{
						ctx->N_filter++;
                        HSVtoRGB(ctx, S, RGB, ctx->InvertValue, ctx->InvertHue);//metka
					    (void)Unitf(S,S);
					    // HSVtoRGB( S, RGB, InvertValue, InvertHue);

						j++;
						for (int k=0; k<Kinp/3; k++) // k runs over ctx->vertices of the box 
						{
							i = j*Kinp + 3*k;	// vertex index
                            N = n*ctx->AtomsPerBlock;//first index of the atom in the blok defines visible/invisible
							Vout[i+0] = (Vinp[3*k+0] + ctx->BlockPosX[n])*ctx->Kind[N];
							Vout[i+1] = (Vinp[3*k+1] + ctx->BlockPosY[n])*ctx->Kind[N];
							Vout[i+2] = (Vinp[3*k+2] + ctx->BlockPosZ[n])*ctx->Kind[N];	

							Nout[i+0] = Ninp[3*k+0];
							Nout[i+1] = Ninp[3*k+1];
							Nout[i+2] = Ninp[3*k+2];

							Cout[i+0] = RGB[0];			
							Cout[i+1] = RGB[1];			
							Cout[i+2] = RGB[2];	
						}
					}//if ( (S[2]>=...
				}//if (F)	
			}
			}
			}
			//IdNum=j*Kinp;
			// for (i=i+1;i<Kout/Kinp;i++){
			// 	Vout[i+0] = Vout[i+1] = Vout[i+2] = 0;
			// 	Nout[i+0] = Nout[i+1] = Nout[i+2] = 0;
			// 	Cout[i+0] = Cout[i+1] = Cout[i+2] = 0;
			// }
		}else{
			for (int an = anini; an<anfin; an++) {
			for (int bn = bnini; bn<bnfin; bn++) {
			for (int cn = cnini; cn<cnfin; cn++) {	
				if (an == anini || an == anfin-1 || bn == bnini || bn == bnfin-1 || cn == cnini || cn == cnfin-1){	
					S[0] = 0;
					S[1] = 0;
					S[2] = 0;

					n = an+bn*ctx->uABC[0]+cn*ctx->uABC[0]*ctx->uABC[1];// index of the block
					for (int atom=0; atom<ctx->AtomsPerBlock; atom++)//atom is the index of the atom in block
					{
					    N = atom + n*ctx->AtomsPerBlock;//n*AtomsPerBlock=index of the first spin in the block;
		 
						S[0]+= VEC_X(spinArr,N);
						S[1]+= VEC_Y(spinArr,N);
						S[2]+= VEC_Z(spinArr,N);
				    }
                    HSVtoRGB(ctx, S, RGB, ctx->InvertValue, ctx->InvertHue);//metka
				    (void)Unitf(S,S);
				    // HSVtoRGB( S, RGB, InvertValue, InvertHue);

					j++;
					for (int k=0; k<Kinp/3; k++) // k runs over ctx->vertices of the box 
					{
						i = j*Kinp + 3*k;	// vertex index
                        N=n*ctx->AtomsPerBlock;//first index of the atom in the blok defines visible/invisible
						Vout[i+0] = (Vinp[3*k+0] + ctx->BlockPosX[n])*ctx->Kind[N];
						Vout[i+1] = (Vinp[3*k+1] + ctx->BlockPosY[n])*ctx->Kind[N];
						Vout[i+2] = (Vinp[3*k+2] + ctx->BlockPosZ[n])*ctx->Kind[N];	

						Nout[i+0] = Ninp[3*k+0];
						Nout[i+1] = Ninp[3*k+1];
						Nout[i+2] = Ninp[3*k+2];

						Cout[i+0] = RGB[0];			
						Cout[i+1] = RGB[1];			
						Cout[i+2] = RGB[2];	
					}	
				}	
			}
			}
			}	
		}
	break;

	case uPOINT:
		if (ctx->WhichSliceMode==FILTER)
		{
			ctx->N_filter=0;
			for (int an = 0; an<ctx->uABC[0]; an++) // n runs over spins 
			{
			for (int bn = 0; bn<ctx->uABC[1]; bn++) // n runs over spins 
			{
			for (int cn = 0; cn<ctx->uABC[2]; cn++) // n runs over spins 
			{
				if (ctx->GreedFilter){
					F =	(an>=ctx->GreedFilterMinA && an<=ctx->GreedFilterMaxA) &&
						(bn>=ctx->GreedFilterMinB && bn<=ctx->GreedFilterMaxB) &&
						(cn>=ctx->GreedFilterMinC && cn<=ctx->GreedFilterMaxC) ;
					if (ctx->GreedFilterInvert) F=!F;					
				}else{F=true;}

				if (F)
				{
					n = an+bn*ctx->uABC[0]+cn*ctx->uABC[0]*ctx->uABC[1];// index of the block
					n = n*ctx->AtomsPerBlock;//index of the first spin in the block
					for (int atom=0; atom<ctx->AtomsPerBlock; atom++)//atom is the index of the atom in block
					{
						S[0] = VEC_X(spinArr,n+atom);
						S[1] = VEC_Y(spinArr,n+atom);
						S[2] = VEC_Z(spinArr,n+atom);
                        HSVtoRGB(ctx, S, RGB, ctx->InvertValue, ctx->InvertHue);
						vlength = Unitf(S,S);//metka
						int phi = atan2int( S[1], S[0] );// return integer angle phi 0 - 360
						//if (S[2]>=Sz_min1 && S[2]<=Sz_max1)
						//if ( (S[2]>=Sz_min1 && S[2]<=Sz_max1) && (phi>=phi_min1 && phi<=phi_max1) )
						F =	(ctx->SpinFilter1 && ((S[2]>=ctx->Sz_min1 && S[2]<=ctx->Sz_max1) && (phi>=ctx->phi_min1 && phi<=ctx->phi_max1)) ) ||
							(ctx->SpinFilter2 && ((S[2]>=ctx->Sz_min2 && S[2]<=ctx->Sz_max2) && (phi>=ctx->phi_min2 && phi<=ctx->phi_max2)) ) ||
							(ctx->SpinFilter3 && ((S[2]>=ctx->Sz_min3 && S[2]<=ctx->Sz_max3) && (phi>=ctx->phi_min3 && phi<=ctx->phi_max3)) ) ;
						if (F)
						{
							ctx->N_filter++;
					        j++;
							i = j*Kinp;			// index of first cane vertex 
							Vout[i+0] = Px[n+atom];	// new x-component of vertex + translation
							Vout[i+1] = Py[n+atom];	// new y-component of vertex + translation
							Vout[i+2] = Pz[n+atom];	// new z-component of vertex + translation
							Cout[i+0] = RGB[0]*ctx->Kind[n+atom];	// x-component of vertex normal
							Cout[i+1] = RGB[1]*ctx->Kind[n+atom];	// y-component of vertex normal
							Cout[i+2] = RGB[2]*ctx->Kind[n+atom];	// z-component of vertex normal
					    }
					}
				}//if (F)
			}
			}
			}
		}else{// slice mode
			for (int an = anini; an<anfin; an++) // n runs over spins 
			{
			for (int bn = bnini; bn<bnfin; bn++) // n runs over spins 
			{
			for (int cn = cnini; cn<cnfin; cn++) // n runs over spins 
			{
				n = an+bn*ctx->uABC[0]+cn*ctx->uABC[0]*ctx->uABC[1];// index of the block
				n = n*ctx->AtomsPerBlock;//index of the first spin in the block
				for (int atom=0; atom<ctx->AtomsPerBlock; atom++)//atom is the index of the atom in block
				{
					S[0] = VEC_X(spinArr,n+atom);
					S[1] = VEC_Y(spinArr,n+atom);
					S[2] = VEC_Z(spinArr,n+atom);
                    HSVtoRGB(ctx, S, RGB, ctx->InvertValue, ctx->InvertHue);//metka
					vlength = Unitf(S,S);
			        // HSVtoRGB( S, RGB, InvertValue, InvertHue);
			        j++;
					i = j*Kinp;			// index of first cane vertex 
					int Factor = ctx->Kind[n+atom];
					if (Factor==0) {
					Vout[i+0] = 1000000;	// new x-component of vertex + translation
					Vout[i+1] = 1000000;	// new y-component of vertex + translation
					Vout[i+2] = 1000000;	// new z-component of vertex + translation
					}else{
					Vout[i+0] = Px[n+atom];	// new x-component of vertex + translation
					Vout[i+1] = Py[n+atom];	// new y-component of vertex + translation
					Vout[i+2] = Pz[n+atom];	// new z-component of vertex + translation
					}
					Cout[i+0] = RGB[0];	// x-component of vertex normal
					Cout[i+1] = RGB[1];	// y-component of vertex normal
					Cout[i+2] = RGB[2];	// z-component of vertex normal
				}	
			}
			}
			}
		}
	break;

	case CANE:
	default:
		if (ctx->WhichSliceMode==FILTER)
		{
			ctx->N_filter=0;
			for (int an = 0; an<ctx->uABC[0]; an++) // an runs over slice along A_AXIS
			{
			for (int bn = 0; bn<ctx->uABC[1]; bn++) // bn runs over slice along B_AXIS
			{
			for (int cn = cnini; cn<ctx->uABC[2]; cn++) // cn runs over slice along C_AXIS
			{
				if (ctx->GreedFilter){
					F =	(an>=ctx->GreedFilterMinA && an<=ctx->GreedFilterMaxA) &&
						(bn>=ctx->GreedFilterMinB && bn<=ctx->GreedFilterMaxB) &&
						(cn>=ctx->GreedFilterMinC && cn<=ctx->GreedFilterMaxC) ;
					if (ctx->GreedFilterInvert) F=!F;					
				}else{F=true;}

				if (F)
				{
					n = an+bn*ctx->uABC[0]+cn*ctx->uABC[0]*ctx->uABC[1];// index of the block
					n = n*ctx->AtomsPerBlock;//index of the first spin in the block
					for (int atom=0; atom<ctx->AtomsPerBlock; atom++)//atom is the index of the atom in block
					{
						S[0] = VEC_X(spinArr,n+atom);
						S[1] = VEC_Y(spinArr,n+atom);
						S[2] = VEC_Z(spinArr,n+atom);
						int phi = atan2int( S[1], S[0] );// return integer angle phi 0 - 360
						//if (S[2]>=Sz_min1 && S[2]<=Sz_max1)
						//if ( (S[2]>=Sz_min1 && S[2]<=Sz_max1) && (phi>=phi_min1 && phi<=phi_max1) )
						F =	(ctx->SpinFilter1 && ((S[2]>=ctx->Sz_min1 && S[2]<=ctx->Sz_max1) && (phi>=ctx->phi_min1 && phi<=ctx->phi_max1)) ) ||
							(ctx->SpinFilter2 && ((S[2]>=ctx->Sz_min2 && S[2]<=ctx->Sz_max2) && (phi>=ctx->phi_min2 && phi<=ctx->phi_max2)) ) ||
							(ctx->SpinFilter3 && ((S[2]>=ctx->Sz_min3 && S[2]<=ctx->Sz_max3) && (phi>=ctx->phi_min3 && phi<=ctx->phi_max3)) ) ;	
						if (F)
						{
                            HSVtoRGB(ctx, S, RGB, ctx->InvertValue, ctx->InvertHue);
                            vlength = Unitf(S,S);//metka
							ctx->N_filter++;
					        j++;
							//i = (n-nini)*Kinp;							// index of ferst cane vertex 
							i = j*Kinp;
							Vout[i+0] = ctx->Kind[n+atom]*( S[0]*(1-ctx->Pivot)*ctx->Scale*vlength + Px[n+atom]);	// new x-component of vertex + translation
							Vout[i+1] = ctx->Kind[n+atom]*( S[1]*(1-ctx->Pivot)*ctx->Scale*vlength + Py[n+atom]);	// new y-component of vertex + translation
							Vout[i+2] = ctx->Kind[n+atom]*( S[2]*(1-ctx->Pivot)*ctx->Scale*vlength + Pz[n+atom]);	// new z-component of vertex + translation
							//i = n*Kinp/3*4;		// ctx->colors contains 4 floats
							Cout[i+0] = RGB[0];					// x-component of vertex normal
							Cout[i+1] = RGB[1];					// y-component of vertex normal
							Cout[i+2] = RGB[2];					// z-component of vertex normal
							//Cout[i+3] = 1.f;
							//printf( "|V1=%f,%f,%f \n",Vout[i+0],Vout[i+1],Vout[i+2]);
							//i = (n-nini)*Kinp + 3*1;					// index of ferst cane vertex 
							i = j*Kinp+ 3*1;
							Vout[i+0] = ctx->Kind[n+atom]*(-S[0]*(ctx->Pivot)*ctx->Scale*vlength + Px[n+atom]);		// new x-component of vertex + translation
							Vout[i+1] = ctx->Kind[n+atom]*(-S[1]*(ctx->Pivot)*ctx->Scale*vlength + Py[n+atom]);		// new y-component of vertex + translation
							Vout[i+2] = ctx->Kind[n+atom]*(-S[2]*(ctx->Pivot)*ctx->Scale*vlength + Pz[n+atom]);		// new z-component of vertex + translation
							//i = n*Kinp/3*4+4;			// ctx->colors contains 4 floats
							Cout[i+0] = RGB[0];//*ctx->Kind[n+atom];					// x-component of vertex normal
							Cout[i+1] = RGB[1];//*ctx->Kind[n+atom];					// y-component of vertex normal
							Cout[i+2] = RGB[2];//*ctx->Kind[n+atom];					// z-component of vertex normal
							//Cout[i+3] = 1.f;
							//printf( "|V1=%f,%f,%f \n",Vout[i+0],Vout[i+1],Vout[i+2]);
						}
					}
				}//if (F)
			}
			}
			}
		}else{//slice mode
			for (int an = anini; an<anfin; an++) // an runs over slice along A_AXIS
			{
			for (int bn = bnini; bn<bnfin; bn++) // bn runs over slice along B_AXIS
			{
			for (int cn = cnini; cn<cnfin; cn++) // cn runs over slice along C_AXIS
			{
				n = an+bn*ctx->uABC[0]+cn*ctx->uABC[0]*ctx->uABC[1];// index of the block
				n = n*ctx->AtomsPerBlock;//index of the first spin in the block
				for (int atom=0; atom<ctx->AtomsPerBlock; atom++)//atom is the index of the atom in block
				{
					S[0] = VEC_X(spinArr,n+atom);
					S[1] = VEC_Y(spinArr,n+atom);
					S[2] = VEC_Z(spinArr,n+atom);
                    HSVtoRGB(ctx, S, RGB, ctx->InvertValue, ctx->InvertHue);//metka
					vlength = Unitf(S,S);
			        // HSVtoRGB( S, RGB, InvertValue, InvertHue);
			        j++;
					//i = (n-nini)*Kinp;							// index of ferst cane vertex 
					i = j*Kinp;
					int Factor = ctx->Kind[n+atom];
					if (Factor==0) Factor=HIDDEN_VECTOR_SCALE;
					Vout[i+0] = Factor*( S[0]*(1-ctx->Pivot)*ctx->Scale*vlength + Px[n+atom]);	// new x-component of vertex + translation
					Vout[i+1] = Factor*( S[1]*(1-ctx->Pivot)*ctx->Scale*vlength + Py[n+atom]);	// new y-component of vertex + translation
					Vout[i+2] = Factor*( S[2]*(1-ctx->Pivot)*ctx->Scale*vlength + Pz[n+atom]);	// new z-component of vertex + translation
					//i = n*Kinp/3*4;		// ctx->colors contains 4 floats
					Cout[i+0] = RGB[0];					// x-component of vertex normal
					Cout[i+1] = RGB[1];					// y-component of vertex normal
					Cout[i+2] = RGB[2];					// z-component of vertex normal
					//Cout[i+3] = 1.f;
					//printf( "|V1=%f,%f,%f \n",Vout[i+0],Vout[i+1],Vout[i+2]);
					//i = (n-nini)*Kinp + 3*1;					// index of ferst cane vertex 
					i = j*Kinp+ 3*1;
					Vout[i+0] = Factor*(-S[0]*(ctx->Pivot)*ctx->Scale*vlength + Px[n+atom]);		// new x-component of vertex + translation
					Vout[i+1] = Factor*(-S[1]*(ctx->Pivot)*ctx->Scale*vlength + Py[n+atom]);		// new y-component of vertex + translation
					Vout[i+2] = Factor*(-S[2]*(ctx->Pivot)*ctx->Scale*vlength + Pz[n+atom]);		// new z-component of vertex + translation
					//i = n*Kinp/3*4+4;			// ctx->colors contains 4 floats
					Cout[i+0] = ctx->Kind[n+atom]*RGB[0];					// x-component of vertex normal
					Cout[i+1] = ctx->Kind[n+atom]*RGB[1];					// y-component of vertex normal
					Cout[i+2] = ctx->Kind[n+atom]*RGB[2];					// z-component of vertex normal
					//Cout[i+3] = 1.f;
					//printf( "|V1=%f,%f,%f \n",Vout[i+0],Vout[i+1],Vout[i+2]);
				}
			}
			}
			}
		}
	}
}

void UpdateVerticesNormalsColors_BextDC(magnoom_ctx *ctx)
{
	const float *Vinp = ctx->vertexProto_BextDC;
	const float *Ninp = ctx->normalProto_BextDC;
	const int Kinp = ctx->VCNum_BextDC;
	float *Vout = ctx->vertices_BextDC;
	float *Nout = ctx->normals_BextDC;
	float *Cout = ctx->colors_BextDC;
	const float Px = ctx->Box[0][0]*0.6;
	const float Py = ctx->Box[1][1]*0.6;
	const float Pz = ctx->Box[2][2]*0.6;
	const float Vx = ctx->BextDCDirection[0];
	const float Vy = ctx->BextDCDirection[1];
	const float Vz = ctx->BextDCDirection[2];
	int i;
    float Sx=Vx;
    float Sy=Vy;
    float Sz=Vz;
	float U,A;
	if (Sz==-1){
		for (int k=0; k<Kinp/3; k++){// k runs over ctx->vertices 
			i = 3*k;	// vertex index
			Vout[i+0] = (-Vinp[i+0])*ctx->BextDCMagnitude*ctx->Scale_BextDC + Px;
			Vout[i+1] = ( Vinp[i+1])*ctx->BextDCMagnitude*ctx->Scale_BextDC + Py;
			Vout[i+2] = (-Vinp[i+2])*ctx->BextDCMagnitude*ctx->Scale_BextDC + Pz;
		
			Nout[i+0] = -Ninp[i+0];
			Nout[i+1] =  Ninp[i+1];
			Nout[i+2] = -Ninp[i+2];
			if (ctx->WhichBackgroundColor == BLACK){
				Cout[i+0] = 0.9;
				Cout[i+1] = 0.9;
				Cout[i+2] = 0.9;
			}else{
				Cout[i+0] = 0.3;
				Cout[i+1] = 0.3;
				Cout[i+2] = 0.3;
			}
		}
	}else{	
		for (int k=0; k<Kinp/3; k++){// k runs over ctx->vertices 
			i = 3*k;	// vertex index
			U = 1.0f/(Sx*Sx + Sy*Sy+(1e-37f)); 		
			A = (-Sy*Vinp[3*k+0] + Sx*Vinp[3*k+1])*(1. - Sz)*U; 
			Vout[i+0] = (-Sy*A + Vinp[3*k+0]*Sz + Sx*Vinp[3*k+2]			)*ctx->BextDCMagnitude*ctx->Scale_BextDC + Px;
			Vout[i+1] = ( Sx*A + Vinp[3*k+1]*Sz + Sy*Vinp[3*k+2]			)*ctx->BextDCMagnitude*ctx->Scale_BextDC + Py;
			Vout[i+2] = ( Vinp[3*k+2]*Sz - (Sx*Vinp[3*k+0]+Sy*Vinp[3*k+1])	)*ctx->BextDCMagnitude*ctx->Scale_BextDC + Pz;

			A = (-Sy*Ninp[3*k+0] + Sx*Ninp[3*k+1])*(1. - Sz)*U; 		
			Nout[i+0] =-Sy*A + Ninp[3*k+0]*Sz + Sx*Ninp[3*k+2];
			Nout[i+1] = Sx*A + Ninp[3*k+1]*Sz + Sy*Ninp[3*k+2];
			Nout[i+2] = Ninp[3*k+2]*Sz - (Sx*Ninp[3*k+0]+Sy*Ninp[3*k+1]);
			if (ctx->WhichBackgroundColor == BLACK){
				Cout[i+0] = 0.9;
				Cout[i+1] = 0.9;
				Cout[i+2] = 0.9;
			}else{
				Cout[i+0] = 0.3;			
				Cout[i+1] = 0.3;
				Cout[i+2] = 0.3;		
			}
		}
	}
}

void
parallelepiped(const float abc[][3], const float tr[3], float scale1, float scale2, float scale3,
	const float color[3], int offset_index, float* V, float* N, float* C, GLuint* I)
{
/*
  p23 o-----------o p123
     /|          /|
 p3 o-----------o p13
  ^ | |p2       | |
 c| | o---------|-o p12 
    |/b         |/     
 p0 o-----------o p1
      -->a  
*/
	float p0[3]={0,0,0};
	
	float p1[3]={abc[0][0],abc[0][1],abc[0][2]}; (void)Unitf(p1,p1);
	p1[0]*=scale1; p1[1]*=scale1; p1[2]*=scale1;
	
	float p2[3]={abc[1][0],abc[1][1],abc[1][2]}; (void)Unitf(p2,p2 );
	p2[0]*=scale2; p2[1]*=scale2; p2[2]*=scale2;
	
	float p3[3]={abc[2][0],abc[2][1],abc[2][2]}; (void)Unitf(p3,p3);
	p3[0]*=scale3; p3[1]*=scale3; p3[2]*=scale3;
	
	float p12[3]={p1[0]+p2[0],p1[1]+p2[1],p1[2]+p2[2]};
	float p13[3]={p1[0]+p3[0],p1[1]+p3[1],p1[2]+p3[2]};
	float p23[3]={p2[0]+p3[0],p2[1]+p3[1],p2[2]+p3[2]};
	float p123[3]={p1[0]+p23[0],p1[1]+p23[1],p1[2]+p23[2]};
	float normal1[3];
	float normal2[3];
	float normal3[3];
	Enorm( p0, p2, p3, normal1);// left and right ~a
	Enorm( p0, p1, p3, normal2);// front and back ~b
	Enorm( p0, p1, p2, normal3);// top and bottom ~c
	float handedness = normal1[0]*p1[0] + normal1[1]*p1[1] + normal1[2]*p1[2] < 0.0f ? -1.0f : 1.0f;
	for (int component = 0; component < 3; ++component) {
		normal1[component] *= handedness;
		normal2[component] *= handedness;
		normal3[component] *= handedness;
	}
	p0[0]  +=tr[0]; p0[1]  +=tr[1]; p0[2]  +=tr[2];
	p1[0]  +=tr[0]; p1[1]  +=tr[1]; p1[2]  +=tr[2];
	p2[0]  +=tr[0]; p2[1]  +=tr[1]; p2[2]  +=tr[2];
	p3[0]  +=tr[0]; p3[1]  +=tr[1]; p3[2]  +=tr[2];
	p12[0] +=tr[0]; p12[1] +=tr[1]; p12[2] +=tr[2];
	p13[0] +=tr[0]; p13[1] +=tr[1]; p13[2] +=tr[2];
	p23[0] +=tr[0]; p23[1] +=tr[1]; p23[2] +=tr[2];
	p123[0]+=tr[0]; p123[1]+=tr[1]; p123[2]+=tr[2];

	int i = offset_index*6*4*3-1;//vertex component counter
	int j = offset_index*6*2*3-1;//vertex index counter 6 sides, 2 triangles, 3 ctx->vertices per triangle
	int k = 0;//used with j

	//top ctx->vertices + ctx->normals, two triangles: p3-p123-p13, p3-p23-p123
	V[++i] = p3[0]; N[i] = normal3[0]; C[i]=color[0];
	V[++i] = p3[1]; N[i] = normal3[1]; C[i]=color[1];
	V[++i] = p3[2]; N[i] = normal3[2]; C[i]=color[2];

	V[++i] = p123[0]; N[i] = normal3[0]; C[i]=color[0];
	V[++i] = p123[1]; N[i] = normal3[1]; C[i]=color[1];
	V[++i] = p123[2]; N[i] = normal3[2]; C[i]=color[2];

	V[++i] = p13[0]; N[i] = normal3[0]; C[i]=color[0];
	V[++i] = p13[1]; N[i] = normal3[1]; C[i]=color[1];
	V[++i] = p13[2]; N[i] = normal3[2]; C[i]=color[2];

	V[++i] = p23[0]; N[i] = normal3[0]; C[i]=color[0];
	V[++i] = p23[1]; N[i] = normal3[1]; C[i]=color[1];
	V[++i] = p23[2]; N[i] = normal3[2]; C[i]=color[2];
	//top ctx->indices p3-p123-p13:
	k = 0;
	I[++j] = offset_index*6*4 + k * 4 + 0; //p3
	I[++j] = offset_index*6*4 + k * 4 + 1; //p123
	I[++j] = offset_index*6*4 + k * 4 + 2; //p13
    //p3-p23-p123
	I[++j] = offset_index*6*4 + k * 4 + 0; //p3
	I[++j] = offset_index*6*4 + k * 4 + 3; //p23
	I[++j] = offset_index*6*4 + k * 4 + 1; //p123

	//bottom ctx->vertices + ctx->normals, two triangles: p0-p1-p12, p0-p12-p2
	V[++i] = p0[0]; N[i] = -normal3[0]; C[i]=color[0];
	V[++i] = p0[1]; N[i] = -normal3[1]; C[i]=color[1];
	V[++i] = p0[2]; N[i] = -normal3[2]; C[i]=color[2];

	V[++i] = p12[0]; N[i] = -normal3[0]; C[i]=color[0];
	V[++i] = p12[1]; N[i] = -normal3[1]; C[i]=color[1];
	V[++i] = p12[2]; N[i] = -normal3[2]; C[i]=color[2];

	V[++i] = p1[0]; N[i] = -normal3[0]; C[i]=color[0];
	V[++i] = p1[1]; N[i] = -normal3[1]; C[i]=color[1];
	V[++i] = p1[2]; N[i] = -normal3[2]; C[i]=color[2];

	V[++i] = p2[0]; N[i] = -normal3[0]; C[i]=color[0];
	V[++i] = p2[1]; N[i] = -normal3[1]; C[i]=color[1];
	V[++i] = p2[2]; N[i] = -normal3[2]; C[i]=color[2];
	//bottom ctx->indices p0-p1-p12:
	k = 1;
	I[++j] = offset_index*6*4 + k * 4 + 2; //p0
	I[++j] = offset_index*6*4 + k * 4 + 1; //p1
	I[++j] = offset_index*6*4 + k * 4 + 0; //p12
    //p0-p12-p2
	I[++j] = offset_index*6*4 + k * 4 + 1; //p0
	I[++j] = offset_index*6*4 + k * 4 + 3; //p12
	I[++j] = offset_index*6*4 + k * 4 + 0; //p2

	//front ctx->vertices + ctx->normals, two triangles: p0-p3-p13, p0-p13-p1
	V[++i] = p0[0]; N[i] = normal2[0]; C[i]=color[0];
	V[++i] = p0[1]; N[i] = normal2[1]; C[i]=color[1];
	V[++i] = p0[2]; N[i] = normal2[2]; C[i]=color[2];
 
	V[++i] = p3[0]; N[i] = normal2[0]; C[i]=color[0];
	V[++i] = p3[1]; N[i] = normal2[1]; C[i]=color[1];
	V[++i] = p3[2]; N[i] = normal2[2]; C[i]=color[2];

	V[++i] = p13[0]; N[i] = normal2[0]; C[i]=color[0];
	V[++i] = p13[1]; N[i] = normal2[1]; C[i]=color[1];
	V[++i] = p13[2]; N[i] = normal2[2]; C[i]=color[2];

	V[++i] = p1[0]; N[i] = normal2[0]; C[i]=color[0];
	V[++i] = p1[1]; N[i] = normal2[1]; C[i]=color[1];
	V[++i] = p1[2]; N[i] = normal2[2]; C[i]=color[2];
	//front ctx->indices p0-p3-p13:
	k = 2;
	I[++j] = offset_index*6*4 + k * 4 + 0; //p0
	I[++j] = offset_index*6*4 + k * 4 + 1; //p3
	I[++j] = offset_index*6*4 + k * 4 + 2; //p13
    //p0-p13-p1
	I[++j] = offset_index*6*4 + k * 4 + 0; //p0
	I[++j] = offset_index*6*4 + k * 4 + 2; //p13
	I[++j] = offset_index*6*4 + k * 4 + 3; //p1

	//back ctx->vertices + ctx->normals, two triangles: p2-p12-p123, p2-p123-p23
	V[++i] = p2[0]; N[i] = -normal2[0]; C[i]=color[0];
	V[++i] = p2[1]; N[i] = -normal2[1]; C[i]=color[1];
	V[++i] = p2[2]; N[i] = -normal2[2]; C[i]=color[2];

	V[++i] = p12[0]; N[i] = -normal2[0]; C[i]=color[0];
	V[++i] = p12[1]; N[i] = -normal2[1]; C[i]=color[1];
	V[++i] = p12[2]; N[i] = -normal2[2]; C[i]=color[2];

	V[++i] = p123[0]; N[i] = -normal2[0]; C[i]=color[0];
	V[++i] = p123[1]; N[i] = -normal2[1]; C[i]=color[1];
	V[++i] = p123[2]; N[i] = -normal2[2]; C[i]=color[2];

	V[++i] = p23[0]; N[i] = -normal2[0]; C[i]=color[0];
	V[++i] = p23[1]; N[i] = -normal2[1]; C[i]=color[1];
	V[++i] = p23[2]; N[i] = -normal2[2]; C[i]=color[2];
	//back ctx->indices p2-p12-p123:
	k = 3;
	I[++j] = offset_index*6*4 + k * 4 + 0; //p2
	I[++j] = offset_index*6*4 + k * 4 + 1; //p12
	I[++j] = offset_index*6*4 + k * 4 + 2; //p123
    //p2-p123-p23
	I[++j] = offset_index*6*4 + k * 4 + 0; //p2
	I[++j] = offset_index*6*4 + k * 4 + 2; //p123
	I[++j] = offset_index*6*4 + k * 4 + 3; //p23
 
	//righ ctx->vertices + ctx->normals, two triangles: p1-p13-p123, p1-p123-p12
	V[++i] = p1[0]; N[i] = normal1[0]; C[i]=color[0];
	V[++i] = p1[1]; N[i] = normal1[1]; C[i]=color[1];
	V[++i] = p1[2]; N[i] = normal1[2]; C[i]=color[2];

	V[++i] = p13[0]; N[i] = normal1[0]; C[i]=color[0];
	V[++i] = p13[1]; N[i] = normal1[1]; C[i]=color[1];
	V[++i] = p13[2]; N[i] = normal1[2]; C[i]=color[2];

	V[++i] = p123[0]; N[i] = normal1[0]; C[i]=color[0];
	V[++i] = p123[1]; N[i] = normal1[1]; C[i]=color[1];
	V[++i] = p123[2]; N[i] = normal1[2]; C[i]=color[2];

	V[++i] = p12[0]; N[i] = normal1[0]; C[i]=color[0];
	V[++i] = p12[1]; N[i] = normal1[1]; C[i]=color[1];
	V[++i] = p12[2]; N[i] = normal1[2]; C[i]=color[2];
	//righ ctx->indices p1-p13-p123:
	k = 4;
	I[++j] = offset_index*6*4 + k * 4 + 0; //p1
	I[++j] = offset_index*6*4 + k * 4 + 1; //p13
	I[++j] = offset_index*6*4 + k * 4 + 2; //p123
    //p1-p123-p12
	I[++j] = offset_index*6*4 + k * 4 + 0; //p1
	I[++j] = offset_index*6*4 + k * 4 + 2; //p123
	I[++j] = offset_index*6*4 + k * 4 + 3; //p12	
 
	//left ctx->vertices + ctx->normals, two triangles: p0-p2-p23, p0-p23-p3
	V[++i] = p0[0]; N[i] = -normal1[0]; C[i]=color[0];
	V[++i] = p0[1]; N[i] = -normal1[1]; C[i]=color[1];
	V[++i] = p0[2]; N[i] = -normal1[2]; C[i]=color[2];

	V[++i] = p2[0]; N[i] = -normal1[0]; C[i]=color[0];
	V[++i] = p2[1]; N[i] = -normal1[1]; C[i]=color[1];
	V[++i] = p2[2]; N[i] = -normal1[2]; C[i]=color[2];

	V[++i] = p23[0]; N[i] = -normal1[0]; C[i]=color[0];
	V[++i] = p23[1]; N[i] = -normal1[1]; C[i]=color[1];
	V[++i] = p23[2]; N[i] = -normal1[2]; C[i]=color[2];

	V[++i] = p3[0]; N[i] = -normal1[0]; C[i]=color[0];
	V[++i] = p3[1]; N[i] = -normal1[1]; C[i]=color[1];
	V[++i] = p3[2]; N[i] = -normal1[2]; C[i]=color[2];
	//left ctx->indices p0-p2-p23:
	k = 5;
	I[++j] = offset_index*6*4 + k * 4 + 0; //p0
	I[++j] = offset_index*6*4 + k * 4 + 1; //p2
	I[++j] = offset_index*6*4 + k * 4 + 2; //p23
    //p0-p23-p3
	I[++j] = offset_index*6*4 + k * 4 + 0; //p0
	I[++j] = offset_index*6*4 + k * 4 + 2; //p23
	I[++j] = offset_index*6*4 + k * 4 + 3; //p3	

	if (handedness < 0.0f) {
		int first_index = offset_index*6*2*3;
		for (int triangle = 0; triangle < 12; ++triangle) {
			GLuint swap = I[first_index + 3*triangle + 1];
			I[first_index + 3*triangle + 1] = I[first_index + 3*triangle + 2];
			I[first_index + 3*triangle + 2] = swap;
		}
	}
}


void UpdateVerticesNormalsColors_BOX(magnoom_ctx *ctx)
{
	float *vertices = ctx->vertices_BOX;
	float *normals = ctx->normals_BOX;
	float *colors = ctx->colors_BOX;
	GLuint *indices = ctx->indices_BOX;
	const float (*box)[3] = ctx->Box;
	float 	d = ctx->WireWidth;
	float 	Tr[3] = {-(box[0][0]+box[1][0]+box[2][0])/2.f,
					 -(box[0][1]+box[1][1]+box[2][1])/2.f,
					 -(box[0][2]+box[1][2]+box[2][2])/2.f 
					};
	Tr[0] -= d/2;
	Tr[1] -= d/2;
	Tr[2] -= d/2;
	float tr[3] = {Tr[0], Tr[1], Tr[2]};
	float length_a = ctx->uABC[0] * sqrt(ctx->abc[0][0]*ctx->abc[0][0] + ctx->abc[0][1]*ctx->abc[0][1] + ctx->abc[0][2]*ctx->abc[0][2])+d;
	float length_b = ctx->uABC[1] * sqrt(ctx->abc[1][0]*ctx->abc[1][0] + ctx->abc[1][1]*ctx->abc[1][1] + ctx->abc[1][2]*ctx->abc[1][2])+d;
	float length_c = ctx->uABC[2] * sqrt(ctx->abc[2][0]*ctx->abc[2][0] + ctx->abc[2][1]*ctx->abc[2][1] + ctx->abc[2][2]*ctx->abc[2][2])+d;
	float color[3]={0.7,0.7,0.7};

	parallelepiped(ctx->abc, tr, length_a, d, d, color, 0, vertices, normals, colors, indices );//(0,0,0)-->(1,0,0)
	parallelepiped(ctx->abc, tr, d, length_b, d, color, 1, vertices, normals, colors, indices );//(0,0,0)-->(0,1,0)
	parallelepiped(ctx->abc, tr, d, d, length_c, color, 2, vertices, normals, colors, indices );//(0,0,0)-->(0,0,1)

	tr[0] = Tr[0]+ctx->abc[1][0]*ctx->uABC[1]; tr[1] = Tr[1]+ctx->abc[1][1]*ctx->uABC[1]; tr[2] = Tr[2]+ctx->abc[1][2]*ctx->uABC[1];
	parallelepiped(ctx->abc, tr, length_a, d, d, color, 3, vertices, normals, colors, indices );//(0,1,0)-->(1,1,0)

	tr[0]=Tr[0]+ctx->abc[2][0]*ctx->uABC[2]; tr[1]=Tr[1]+ctx->abc[2][1]*ctx->uABC[2]; tr[2]=Tr[2]+ctx->abc[2][2]*ctx->uABC[2];
	parallelepiped(ctx->abc, tr, length_a, d, d, color, 4, vertices, normals, colors, indices );//(0,0,1)-->(0,1,1)

	tr[0]+=ctx->abc[1][0]*ctx->uABC[1]; tr[1]+=ctx->abc[1][1]*ctx->uABC[1]; tr[2]+=ctx->abc[1][2]*ctx->uABC[1];
	parallelepiped(ctx->abc, tr, length_a, d, d, color, 5, vertices, normals, colors, indices );//(1,0,1)-->(1,1,1)

	tr[0]=Tr[0]+ctx->abc[0][0]*ctx->uABC[0]; tr[1]=Tr[1]+ctx->abc[0][1]*ctx->uABC[0]; tr[2]=Tr[2]+ctx->abc[0][2]*ctx->uABC[0];
	parallelepiped(ctx->abc, tr, d, length_b, d, color, 6, vertices, normals, colors, indices );//(1,0,0)-->(1,1,0)

	tr[0]=Tr[0]+ctx->abc[2][0]*ctx->uABC[2]; tr[1]=Tr[1]+ctx->abc[2][1]*ctx->uABC[2]; tr[2]=Tr[2]+ctx->abc[2][2]*ctx->uABC[2];
	parallelepiped(ctx->abc, tr, d, length_b, d, color, 7, vertices, normals, colors, indices );//(0,0,1)-->(0,1,1)

	tr[0]+=ctx->abc[0][0]*ctx->uABC[0]; tr[1]+=ctx->abc[0][1]*ctx->uABC[0]; tr[2]+=ctx->abc[0][2]*ctx->uABC[0];
	parallelepiped(ctx->abc, tr, d, length_b, d, color, 8, vertices, normals, colors, indices );//(0,1,1)-->(1,1,1)

	tr[0]=Tr[0]+ctx->abc[0][0]*ctx->uABC[0]; tr[1]=Tr[1]+ctx->abc[0][1]*ctx->uABC[0]; tr[2]=Tr[2]+ctx->abc[0][2]*ctx->uABC[0];
	parallelepiped(ctx->abc, tr, d, d, length_c, color, 9, vertices, normals, colors, indices );//(1,0,0)-->(1,0,1)

	tr[0]=Tr[0]+ctx->abc[1][0]*ctx->uABC[1]; tr[1]=Tr[1]+ctx->abc[1][1]*ctx->uABC[1]; tr[2]=Tr[2]+ctx->abc[1][2]*ctx->uABC[1];
	parallelepiped(ctx->abc, tr, d, d, length_c, color, 10, vertices, normals, colors, indices );//(0,1,0)-->(0,1,1)

	tr[0]+=ctx->abc[0][0]*ctx->uABC[0]; tr[1]+=ctx->abc[0][1]*ctx->uABC[0]; tr[2]+=ctx->abc[0][2]*ctx->uABC[0];
	parallelepiped(ctx->abc, tr, d, d, length_c, color, 11, vertices, normals, colors, indices );//(1,1,0)-->(1,1,1)
}

void UpdateVerticesNormalsColors_BASIS(magnoom_ctx *ctx)
{
	float *vertices = ctx->vertices_BASIS;
	float *normals = ctx->normals_BASIS;
	float *colors = ctx->colors_BASIS;
	GLuint *indices = ctx->indices_BASIS;
	float 	d = ctx->WireWidth;
	float	cube[3][3] = {
				{	1.0f, 0.0f, 0.0f }, // a
				{	0.0f, 1.0f, 0.0f }, // b
				{	0.0f, 0.0f, 1.0f }};// c
	float tr[3] = {-d/2, -d/2, -d/2};
	float length = 20*d;
	float colorR[3]={1, 0, 0};
	float colorG[3]={0, 1, 0};
	float colorB[3]={0, 0, 1};
	float color0[3]={0.1,0.1,0.1};
	
	parallelepiped(cube, tr, length, d, d, colorR, 0, vertices, normals, colors, indices );//X
	parallelepiped(cube, tr, d, length, d, colorG, 1, vertices, normals, colors, indices );//Y
	parallelepiped(cube, tr, d, d, length, colorB, 2, vertices, normals, colors, indices );//Z
	d *= 1.5; tr[0]=-d/2; tr[1]=-d/2; tr[2]=-d/2;

	parallelepiped(cube, tr, d, d, d, color0, 3, vertices, normals, colors, indices );//O

}

void UpdateVerticesNormalsColors_PBC(magnoom_ctx *ctx, int K0)
{
	float *vertices;
	float *normals;
	float *colors;
	GLuint *indices;
	switch (K0) {
		case 0:
			vertices = ctx->vertices_PBC_A;
			normals = ctx->normals_PBC_A;
			colors = ctx->colors_PBC_A;
			indices = ctx->indices_PBC_A;
			break;
		case 1:
			vertices = ctx->vertices_PBC_B;
			normals = ctx->normals_PBC_B;
			colors = ctx->colors_PBC_B;
			indices = ctx->indices_PBC_B;
			break;
		case 2:
			vertices = ctx->vertices_PBC_C;
			normals = ctx->normals_PBC_C;
			colors = ctx->colors_PBC_C;
			indices = ctx->indices_PBC_C;
			break;
		default:
			return;
	}
	const float (*box)[3] = ctx->Box;
	float 	d = ctx->WireWidth;
	float 	Tr[3] = {-(box[0][0]+box[1][0]+box[2][0]+d)/2.f,
					 -(box[0][1]+box[1][1]+box[2][1]+d)/2.f,
					 -(box[0][2]+box[1][2]+box[2][2]+d)/2.f 
					};
	float length_a = d;
	float length_b = d;
	float length_c = d;
	if (K0==0) length_a = 5*d;
	if (K0==1) length_b = 5*d;
	if (K0==2) length_c = 5*d;
	float tr[3];
	float color[3]={0.7,0.7,0.7};
	int K1=(K0+1)%3;
	int K2=(K0+2)%3;

	tr[0] = Tr[0]+ctx->abc[K0][0]*(ctx->uABC[K0]+6*d); 
	tr[1] = Tr[1]+ctx->abc[K0][1]*(ctx->uABC[K0]+6*d);
	tr[2] = Tr[2]+ctx->abc[K0][2]*(ctx->uABC[K0]+6*d);

	parallelepiped(ctx->abc, tr, length_a, length_b, length_c, color, 0, vertices, normals, colors, indices );

	tr[0]=Tr[0]+ctx->abc[K0][0]*(ctx->uABC[K0]+16*d); 
	tr[1]=Tr[1]+ctx->abc[K0][1]*(ctx->uABC[K0]+16*d);
	tr[2]=Tr[2]+ctx->abc[K0][2]*(ctx->uABC[K0]+16*d);
    parallelepiped(ctx->abc, tr, length_a, length_b, length_c, color, 1, vertices, normals, colors, indices );

	tr[0]=Tr[0]+ctx->abc[K0][0]*(-10*d); 
	tr[1]=Tr[1]+ctx->abc[K0][1]*(-10*d);
	tr[2]=Tr[2]+ctx->abc[K0][2]*(-10*d);
	parallelepiped(ctx->abc, tr, length_a, length_b, length_c, color, 2, vertices, normals, colors, indices );

	tr[0]=Tr[0]+ctx->abc[K0][0]*(-20*d); 
	tr[1]=Tr[1]+ctx->abc[K0][1]*(-20*d);
	tr[2]=Tr[2]+ctx->abc[K0][2]*(-20*d);
	parallelepiped(ctx->abc, tr, length_a, length_b, length_c, color, 3, vertices, normals, colors, indices );

	tr[0]=Tr[0]+ctx->abc[K1][0]*ctx->uABC[K1]+ctx->abc[K0][0]*(ctx->uABC[K0]+6*d);
	tr[1]=Tr[1]+ctx->abc[K1][1]*ctx->uABC[K1]+ctx->abc[K0][1]*(ctx->uABC[K0]+6*d);
	tr[2]=Tr[2]+ctx->abc[K1][2]*ctx->uABC[K1]+ctx->abc[K0][2]*(ctx->uABC[K0]+6*d);
	parallelepiped(ctx->abc, tr, length_a, length_b, length_c, color, 4, vertices, normals, colors, indices );

	tr[0]=Tr[0]+ctx->abc[K1][0]*ctx->uABC[K1]+ctx->abc[K0][0]*(ctx->uABC[K0]+16*d);
	tr[1]=Tr[1]+ctx->abc[K1][1]*ctx->uABC[K1]+ctx->abc[K0][1]*(ctx->uABC[K0]+16*d);
	tr[2]=Tr[2]+ctx->abc[K1][2]*ctx->uABC[K1]+ctx->abc[K0][2]*(ctx->uABC[K0]+16*d);
	parallelepiped(ctx->abc, tr, length_a, length_b, length_c, color, 5, vertices, normals, colors, indices );

	tr[0]=Tr[0]+ctx->abc[K1][0]*ctx->uABC[K1]+ctx->abc[K0][0]*(-10*d);
	tr[1]=Tr[1]+ctx->abc[K1][1]*ctx->uABC[K1]+ctx->abc[K0][1]*(-10*d);
	tr[2]=Tr[2]+ctx->abc[K1][2]*ctx->uABC[K1]+ctx->abc[K0][2]*(-10*d);	
	parallelepiped(ctx->abc, tr, length_a, length_b, length_c, color, 6, vertices, normals, colors, indices );

	tr[0]=Tr[0]+ctx->abc[K1][0]*ctx->uABC[K1]+ctx->abc[K0][0]*(-20*d);
	tr[1]=Tr[1]+ctx->abc[K1][1]*ctx->uABC[K1]+ctx->abc[K0][1]*(-20*d);
	tr[2]=Tr[2]+ctx->abc[K1][2]*ctx->uABC[K1]+ctx->abc[K0][2]*(-20*d);	
	parallelepiped(ctx->abc, tr, length_a, length_b, length_c, color, 7, vertices, normals, colors, indices );

	tr[0]=Tr[0]+ctx->abc[K2][0]*ctx->uABC[K2]+ctx->abc[K0][0]*(ctx->uABC[K0]+6*d);
	tr[1]=Tr[1]+ctx->abc[K2][1]*ctx->uABC[K2]+ctx->abc[K0][1]*(ctx->uABC[K0]+6*d);
	tr[2]=Tr[2]+ctx->abc[K2][2]*ctx->uABC[K2]+ctx->abc[K0][2]*(ctx->uABC[K0]+6*d);
	parallelepiped(ctx->abc, tr, length_a, length_b, length_c, color, 8, vertices, normals, colors, indices );

	tr[0]=Tr[0]+ctx->abc[K2][0]*ctx->uABC[K2]+ctx->abc[K0][0]*(ctx->uABC[K0]+16*d);
	tr[1]=Tr[1]+ctx->abc[K2][1]*ctx->uABC[K2]+ctx->abc[K0][1]*(ctx->uABC[K0]+16*d);
	tr[2]=Tr[2]+ctx->abc[K2][2]*ctx->uABC[K2]+ctx->abc[K0][2]*(ctx->uABC[K0]+16*d);

	parallelepiped(ctx->abc, tr, length_a, length_b, length_c, color, 9, vertices, normals, colors, indices );
	tr[0]=Tr[0]+ctx->abc[K2][0]*ctx->uABC[K2]+ctx->abc[K0][0]*(-10*d);
	tr[1]=Tr[1]+ctx->abc[K2][1]*ctx->uABC[K2]+ctx->abc[K0][1]*(-10*d);
	tr[2]=Tr[2]+ctx->abc[K2][2]*ctx->uABC[K2]+ctx->abc[K0][2]*(-10*d);
	parallelepiped(ctx->abc, tr, length_a, length_b, length_c, color, 10, vertices, normals, colors, indices );
	tr[0]=Tr[0]+ctx->abc[K2][0]*ctx->uABC[K2]+ctx->abc[K0][0]*(-20*d);
	tr[1]=Tr[1]+ctx->abc[K2][1]*ctx->uABC[K2]+ctx->abc[K0][1]*(-20*d);
	tr[2]=Tr[2]+ctx->abc[K2][2]*ctx->uABC[K2]+ctx->abc[K0][2]*(-20*d);
	parallelepiped(ctx->abc, tr, length_a, length_b, length_c, color, 11, vertices, normals, colors, indices );

	tr[0]=Tr[0]+ctx->abc[K1][0]*ctx->uABC[K1]+ctx->abc[K2][0]*ctx->uABC[K2]+ctx->abc[K0][0]*(ctx->uABC[K0]+6*d);
	tr[1]=Tr[1]+ctx->abc[K1][1]*ctx->uABC[K1]+ctx->abc[K2][1]*ctx->uABC[K2]+ctx->abc[K0][1]*(ctx->uABC[K0]+6*d);
	tr[2]=Tr[2]+ctx->abc[K1][2]*ctx->uABC[K1]+ctx->abc[K2][2]*ctx->uABC[K2]+ctx->abc[K0][2]*(ctx->uABC[K0]+6*d);
	parallelepiped(ctx->abc, tr, length_a, length_b, length_c, color, 12, vertices, normals, colors, indices );

	tr[0]=Tr[0]+ctx->abc[K1][0]*ctx->uABC[K1]+ctx->abc[K2][0]*ctx->uABC[K2]+ctx->abc[K0][0]*(ctx->uABC[K0]+16*d);
	tr[1]=Tr[1]+ctx->abc[K1][1]*ctx->uABC[K1]+ctx->abc[K2][1]*ctx->uABC[K2]+ctx->abc[K0][1]*(ctx->uABC[K0]+16*d);
	tr[2]=Tr[2]+ctx->abc[K1][2]*ctx->uABC[K1]+ctx->abc[K2][2]*ctx->uABC[K2]+ctx->abc[K0][2]*(ctx->uABC[K0]+16*d);
	parallelepiped(ctx->abc, tr, length_a, length_b, length_c, color, 13, vertices, normals, colors, indices );

	tr[0]=Tr[0]+ctx->abc[K1][0]*ctx->uABC[K1]+ctx->abc[K2][0]*ctx->uABC[K2]+ctx->abc[K0][0]*(-10*d);
	tr[1]=Tr[1]+ctx->abc[K1][1]*ctx->uABC[K1]+ctx->abc[K2][1]*ctx->uABC[K2]+ctx->abc[K0][1]*(-10*d);
	tr[2]=Tr[2]+ctx->abc[K1][2]*ctx->uABC[K1]+ctx->abc[K2][2]*ctx->uABC[K2]+ctx->abc[K0][2]*(-10*d);
	parallelepiped(ctx->abc, tr, length_a, length_b, length_c, color, 14, vertices, normals, colors, indices );

	tr[0]=Tr[0]+ctx->abc[K1][0]*ctx->uABC[K1]+ctx->abc[K2][0]*ctx->uABC[K2]+ctx->abc[K0][0]*(-20*d);
	tr[1]=Tr[1]+ctx->abc[K1][1]*ctx->uABC[K1]+ctx->abc[K2][1]*ctx->uABC[K2]+ctx->abc[K0][1]*(-20*d);
	tr[2]=Tr[2]+ctx->abc[K1][2]*ctx->uABC[K1]+ctx->abc[K2][2]*ctx->uABC[K2]+ctx->abc[K0][2]*(-20*d);
	parallelepiped(ctx->abc, tr, length_a, length_b, length_c, color, 15, vertices, normals, colors, indices );
}

void UpdateSpinPositions(magnoom_ctx *ctx)
{
	const float (*abc)[3] = ctx->abc;
	const int *uABC = ctx->uABC;
	const float (*BD)[3] = ctx->Block;
	const int NBD = ctx->AtomsPerBlock;
	const float (*box)[3] = ctx->Box;
	float *Px = ctx->PosX;
	float *Py = ctx->PosY;
	float *Pz = ctx->PosZ;
	float Tr[3] = {
					(box[0][0]+box[1][0]+box[2][0])/2.f,
					(box[0][1]+box[1][1]+box[2][1])/2.f,
					(box[0][2]+box[1][2]+box[2][2])/2.f
				  };
	// float Tr[3] = {	
	// 				(box[0][0]-abc[0][0]+box[1][0]-abc[1][0]+box[2][0]-abc[2][0])/2.f,
	// 				(box[0][1]-abc[0][1]+box[1][1]-abc[1][1]+box[2][1]-abc[2][1])/2.f,
	// 				(box[0][2]-abc[0][2]+box[1][2]-abc[1][2]+box[2][2]-abc[2][2])/2.f
	// 			  };
	int si=-1; // spin index
	int bi=-1; // block index
	for( int L=0;L<uABC[2];L++)// translation of basic domain along vector 'c' L times
	{
		for(int K=0;K<uABC[1];K++)// translation of basic domain along vector 'b' K times
		{
			for(int J=0;J<uABC[0];J++) // translation of basic domain along vector 'a' J times
			{	
				for(int I=0; I < NBD; I++) // runs over atoms in basic domain 
				{	
					si++; 
					Px[si] = BD[I][0] + abc[0][0]*J + abc[1][0]*K + abc[2][0]*L-Tr[0]; 
					Py[si] = BD[I][1] + abc[0][1]*J + abc[1][1]*K + abc[2][1]*L-Tr[1];
					Pz[si] = BD[I][2] + abc[0][2]*J + abc[1][2]*K + abc[2][2]*L-Tr[2];
				}
				bi++;
				ctx->BlockPosX[bi] = abc[0][0]*J + abc[1][0]*K + abc[2][0]*L-Tr[0]; 
				ctx->BlockPosY[bi] = abc[0][1]*J + abc[1][1]*K + abc[2][1]*L-Tr[1];
				ctx->BlockPosZ[bi] = abc[0][2]*J + abc[1][2]*K + abc[2][2]*L-Tr[2];				
			}
		}
	}	
}
///////////////////////////////////////////////////////////////////////////////////////////////////






void drawVBO(magnoom_ctx *ctx)
{
	GLsizei spin_count = (GLsizei)ctx->spin_mesh.index_count;
	int is_solid = (ctx->WhichVectorMode == BOX1 || ctx->WhichVectorMode == ARROW1 || ctx->WhichVectorMode == CONE1);

	if (is_solid) {
		ctx->spin_mesh.uses_normals = 1;
		DrawVBOMeshIndexed(&ctx->spin_mesh, GL_TRIANGLES, spin_count, 0);
		return;
	}

	GLboolean restore_lighting = glIsEnabled(GL_LIGHTING);
	glDisable(GL_LIGHTING);
	ctx->spin_mesh.uses_normals = 0;

	if (ctx->WhichVectorMode == uPOINT) {
		DrawVBOMeshIndexed(&ctx->spin_mesh, GL_POINTS, spin_count, 0);
		glPointSize(10.f * ctx->Scale);
		glDrawElements(GL_POINTS, spin_count, GL_UNSIGNED_INT, (void*)0);
		if (restore_lighting) glEnable(GL_LIGHTING);
		return;
	}

	/* CANE (default) mode — draw twice: first lines+small points (stride 0),
	   then large points (interleaved stride 6*sizeof(float)). */
	{
		vbo_mesh *m = &ctx->spin_mesh;
		GLsizei stride0  = 0;
		GLsizei stride6  = 6 * (GLsizei)sizeof(float);

		/* Pass 1 — stride 0 */
		glBindBuffer(GL_ARRAY_BUFFER, m->color_buffer);
		glColorPointer(3, GL_FLOAT, stride0, (void*)0);
		glBindBuffer(GL_ARRAY_BUFFER, m->vertex_buffer);
		glVertexPointer(3, GL_FLOAT, stride0, (void*)0);
		glEnableClientState(GL_COLOR_ARRAY);
		glEnableClientState(GL_VERTEX_ARRAY);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m->index_buffer);

		glLineWidth(5.0f * ctx->Scale);
		glDrawElements(GL_LINES, spin_count, GL_UNSIGNED_INT, (void*)0);
		glPointSize(4.0f * ctx->Scale);
		glDrawElements(GL_POINTS, spin_count, GL_UNSIGNED_INT, (void*)0);

		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_COLOR_ARRAY);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

		/* Pass 2 — interleaved stride 6 */
		glBindBuffer(GL_ARRAY_BUFFER, m->color_buffer);
		glColorPointer(3, GL_FLOAT, stride6, (void*)0);
		glBindBuffer(GL_ARRAY_BUFFER, m->vertex_buffer);
		glVertexPointer(3, GL_FLOAT, stride6, (void*)0);
		glEnableClientState(GL_COLOR_ARRAY);
		glEnableClientState(GL_VERTEX_ARRAY);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m->index_buffer);

		glPointSize(10.5f * ctx->Scale);
		glDrawElements(GL_POINTS, spin_count / 2, GL_UNSIGNED_INT, (void*)0);

		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_COLOR_ARRAY);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

		if (restore_lighting) glEnable(GL_LIGHTING);
	}
}




//////////////////////////////////////////////////////////////////////////////////////////
// vbo_mesh generic helpers                                                             //
//////////////////////////////////////////////////////////////////////////////////////////

void InitVBOMesh(vbo_mesh *mesh, GLenum usage)
{
	mesh->vertex_buffer = 0;
	mesh->normal_buffer = 0;
	mesh->color_buffer  = 0;
	mesh->index_buffer  = 0;
	mesh->component_count    = 0;
	mesh->component_capacity = 0;
	mesh->index_count    = 0;
	mesh->index_capacity = 0;
	mesh->usage = usage;
}

void CreateVBOMesh(vbo_mesh *mesh)
{
	glGenBuffers(1, &mesh->vertex_buffer);
	glGenBuffers(1, &mesh->normal_buffer);
	glGenBuffers(1, &mesh->color_buffer);
	glGenBuffers(1, &mesh->index_buffer);
}

void UploadVBOMesh(vbo_mesh *restrict mesh,
		   const float *restrict vertices,
		   const float *restrict normals,
		   const float *restrict colors,
		   const GLuint *restrict indices,
		   unsigned int upload_mask)
{
	size_t vert_bytes = mesh->component_capacity * sizeof(float);
	size_t idx_bytes  = mesh->index_capacity  * sizeof(GLuint);

	if (upload_mask & VBO_UPLOAD_VERTICES) {
		glBindBuffer(GL_ARRAY_BUFFER, mesh->vertex_buffer);
		glBufferData(GL_ARRAY_BUFFER, vert_bytes, NULL, mesh->usage);
		glBufferSubData(GL_ARRAY_BUFFER, 0, mesh->component_count * sizeof(float), vertices);
	}
	if (upload_mask & VBO_UPLOAD_NORMALS) {
		glBindBuffer(GL_ARRAY_BUFFER, mesh->normal_buffer);
		glBufferData(GL_ARRAY_BUFFER, vert_bytes, NULL, mesh->usage);
		glBufferSubData(GL_ARRAY_BUFFER, 0, mesh->component_count * sizeof(float), normals);
	}
	if (upload_mask & VBO_UPLOAD_COLORS) {
		glBindBuffer(GL_ARRAY_BUFFER, mesh->color_buffer);
		glBufferData(GL_ARRAY_BUFFER, vert_bytes, NULL, mesh->usage);
		glBufferSubData(GL_ARRAY_BUFFER, 0, mesh->component_count * sizeof(float), colors);
	}
	if (upload_mask & VBO_UPLOAD_INDICES) {
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh->index_buffer);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, idx_bytes, NULL, mesh->usage);
		glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, mesh->index_capacity * sizeof(GLuint), indices);
	}
}

void DrawVBOMeshIndexed(const vbo_mesh *mesh, GLenum primitive, GLsizei count, GLsizei stride)
{
	glBindBuffer(GL_ARRAY_BUFFER, mesh->color_buffer);
	glColorPointer(3, GL_FLOAT, (GLsizei)stride, (void*)0);

	if (mesh->uses_normals) {
		glBindBuffer(GL_ARRAY_BUFFER, mesh->normal_buffer);
		glNormalPointer(GL_FLOAT, 0, (void*)0);
	}

	glBindBuffer(GL_ARRAY_BUFFER, mesh->vertex_buffer);
	glVertexPointer(3, GL_FLOAT, (GLsizei)stride, (void*)0);

	glEnableClientState(GL_COLOR_ARRAY);
	if (mesh->uses_normals) glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_VERTEX_ARRAY);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh->index_buffer);
	glDrawElements(primitive, count, GL_UNSIGNED_INT, (void*)0);

	glDisableClientState(GL_VERTEX_ARRAY);
	if (mesh->uses_normals) glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void DestroyVBOMesh(vbo_mesh *mesh)
{
	if (mesh->vertex_buffer) glDeleteBuffers(1, &mesh->vertex_buffer);
	if (mesh->normal_buffer) glDeleteBuffers(1, &mesh->normal_buffer);
	if (mesh->color_buffer)  glDeleteBuffers(1, &mesh->color_buffer);
	if (mesh->index_buffer)  glDeleteBuffers(1, &mesh->index_buffer);
	InitVBOMesh(mesh, mesh->usage);
}
