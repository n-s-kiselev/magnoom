GLuint 			GLUT_window;
const char *	WINDOWTITLE = { "Magnoom v1.0" };
int				window_width	= 800;//1024;//400;//
int				window_height	= 600;//768;//300;//
float 			asp_rat			= (float)( ((double)window_width)/((double)window_height) );
float 			asp_rat_inv		= (float)( ((double)window_height)/((double)window_width) );

float			CameraEye[3]	= { 0.0, 0.0, 50.0}; // "camera position"
float			CameraC[3]		= { 0.0, 0.0,  0.0}; // "look at point"
float			CameraUp[3]		= { 0.0, 1.0,  0.0}; // "where is up direction"

// light parameters
GLfloat			light_ambient[]  = {0.1, 0.1, 0.1, 1.0};
GLfloat			light_diffuse[]  = {0.9, 0.9, 0.9, 1.0};
GLfloat			light_specular[] = {0.0, 0.0, 0.0, 1.0};

GLfloat			light_position[] = {0.0, 0.0, 1000.0, 0.0};
GLfloat			light_direction[] = {0.0, 0.0, 10.0, 0.0};

// material parameters
GLfloat			material_ambient[]  = {0.8, 0.8, 0.8, 1.0};
GLfloat			material_diffuse[]  = {0.2, 0.2, 0.2, 1.0};
GLfloat			material_specular[] = {0.0, 0.0, 0.0, 1.0};

float			PerspSet[4]		= {60.0, asp_rat, 0.01, 50000}; // {Setings: Field of view vertical, apect ratio, zNear, zFar}  

typedef enum	{ORTHO,	PERSP} enProjections;	// declare new enum type for projections
enProjections	WhichProjection = PERSP; // PERSP by default 	

typedef enum	{RND, HOMO, SKYRM1, SKYRM2, SKYRM3, BOBBER_T, BOBBER_B, BOBBER_L, BOBBER_L_T, BOBBER_L_B, HOPFION1, SPIRAL, SKYRMION_L} enIniState; // which mode
enIniState		WhichInitialState = RND;	// RND by default 

// what the glui package defines as true and false:
const int 		GLUITRUE  = { true  };
const int 		GLUIFALSE = { false };
// active mouse buttons (or them together):
const int 		LEFT   = { 4 };
const int 		MIDDLE = { 2 };
const int 		RIGHT  = { 1 };

// which button:
typedef enum 	{XUP, YUP, ZUP, ADD_SET, RESET, QUIT, PLAY, RECORD} enButton;

// the color numbers, this order must match the radio button order
typedef enum	{WHITE, BLACK, RED, GREEN, BLUE} enColors;
enColors		WhichColor = BLACK;	// index into Colors[ ]

typedef enum	{RGB, RYGB} enColorScheme;
enColorScheme	WhichColorScheme = RGB;

int temp_color[3] = { 55, 55, 155 };

// the color definitions, this order must match the radio button order
const GLfloat 	Colors[ ][3] = 
{
				{ 0.975, 0.975, 0.975 },		// white
				{ 0.1, 0.1, 0.1 },	// black
				{ 1., 0.8, 0.8 },	// red
				{ 0.8, 1., 0.8 },	// green
				{ 0.8, 0.8, 1. },	// blue
};

typedef enum	{ARROW1 
				,CONE1 
				,CANE
				,uPOINT 
				,BOX1
				} enVectorMode; // which mode
enVectorMode	WhichVectorMode	= BOX1;	// CANE by default 
char shortBufer[80];
char inputfilename[64] = "input.csv";
char outputfilename[64] = "output.csv";

#define ESCAPE  0x1b // the escape key:
#define SPACE   0x20 // the space  key:
#define PLUS    0x2b // "+" key
#define MINUS   0x2d // "-" key

//Color map control parameters
int				ColorShift=0; // by default the red collor corresponds to phi=0
int				InvertHue=0; // by default the sequence is red-green-blue for phi=0,120,240(-120) degrees respectively
int				InvertValue=0; //n_z=+1 (white), -1 (black). For InvertHue=1 vice versa 

// non-constant global variables:
int				ActiveButton;	// current mous button that is down
int				Xmouse, Ymouse;	// mouse values
float			Rot[3]={0,0,0};	// rotation angles in degrees
float			TransXYZ[3]={0,0,0};	// set by glui translation widgets
const int       NumCamPosSave=5;
int             CurrentCameraPositionBank=0;
float           CameraPosition[NumCamPosSave][7];// array which contains camera positions 
float			Scale = 1.f;	// scaling factors for arrows [0.1-2] 
float			Pivot = 0.55f;
float 			WireWidth = 0.2;

float			Scale_H = (float)(uABC[0]+uABC[1]+uABC[2]);	// scaling factors for arrows [0.1-2] 

// Slicing parameters
int				A_layer_min = 1;	// which layer (along a tr. vect. ) to show max=uABC[0]
int				B_layer_min = 1;	// which layer (along b tr. vect. ) to show max=uABC[1]
int				C_layer_min = 1;	// which layer (along c tr. vect. ) to show max=uABC[2]
int				A_layer_max = 1;	// which layer (along a tr. vect. ) to show max=uABC[0]
int				B_layer_max = 1;	// which layer (along b tr. vect. ) to show max=uABC[1]
int				C_layer_max = 1;	// which layer (along c tr. vect. ) to show max=uABC[2]

typedef enum	{A_AXIS, B_AXIS, C_AXIS, FILTER} enSliceMode; // which mode
enSliceMode	    WhichSliceMode	= C_AXIS;	// CANE by default 

int   N_filter=0;
float theta_max1=PI/2+0.13; //0.01;//
float Sz_min1=cos(theta_max1);
float theta_min1=PI/2-0.13; //0;//   
float Sz_max1=cos(theta_min1);

float theta_max2=PI/4+0.13; //0.01;//
float Sz_min2=cos(theta_max2);
float theta_min2=PI/4-0.13; //0;//   
float Sz_max2=cos(theta_min2);

int phi_max1=360;
int phi_min1=0;

// Parameters for initial state 
float			chSize = 35; // characteristic size of initial state in units of "a"
float			chDir[3] = {0,1,0}; // characteristic size of initial state in units of "a"
GLuint 			iStart;
GLuint 			iNum;

// spin textue action
float			RotateAllSpins = 0.0f;

// axes parameters:
int				AxesOn=1;		// != 0 means to draw the axes
int				BoxOn=1;		// != 0 means to draw the axes
GLuint			AxesList, BoxList; // OpenGL lists to hold the scenario for the drawing something
GLuint			BoundaryListA, BoundaryListB, BoundaryListC;
GLfloat			AXES_WIDTH	= 2.;
GLfloat			AXES_LENGTH	= 2.;

///////////////////////////////////////////////////////////////////////////////////////////////////
//tweak menu
TwBar *help_bar; // Pointer to the default help tweak bar
TwBar *view_bar; // Pointer to the tweak bar with widgets controling view options
TwBar *control_bar; // Pointer to the tweak bar with widgets controling calculations
TwBar *initial_bar; // Pointer to the tweak bar with widgets controling generation of initial state
TwBar *ac_field_bar; // Pointer to the tweak bar with widgets controling generation of initial state
TwBar *info_bar; // Pointer to the tweak bar with widgets controling generation of initial state
float axis[] = { 0.7f, 0.7f, 0.0f }; // initial model rotation
float angle = 0.8f;

// Shapes material
// float g_MatAmbient[] = { 0.5f, 0.0f, 0.0f, 1.0f };
// float g_MatDiffuse[] = { 1.0f, 1.0f, 0.0f, 1.0f };
// Light parameter
float g_LightMultiplier = 1.0f;
float g_LightDirection[] = { -0.2f, 1.0f, 0.0f };
int   Light_On=1;



////////////////////////////////// Parameters for visualisation ///////////////////////////////////

GLfloat*	vertexProto   	= NULL; // array of vertexes for prototipe arrow or cane
GLfloat*	normalProto   	= NULL; // array of normals for prototipe arrow (not used in vector "cane mode")
GLuint*		indicesProto  	= NULL; // array of indices for prototipe arrow or cane

GLfloat*	vertexProto_H 	= NULL; // array of vertexes for prototipe arrow or cane
GLfloat*	normalProto_H 	= NULL; // array of normals for prototipe arrow (not used in vector "cane mode")
GLuint*		indicesProto_H	= NULL; // array of indices for prototipe arrow or cane

GLfloat*	vertices		= NULL; // array of vertexes for tatal vector field 
GLfloat*	normals			= NULL; // array of normals for tatal vector field
GLfloat*	colors			= NULL; // array of colors 
GLuint*		indices			= NULL; // array of indices for tatal vector field

GLfloat*	vertices_H		= NULL; // array of vertexes for tatal vector field 
GLfloat*	normals_H		= NULL; // array of normals for tatal vector field
GLfloat*	colors_H		= NULL; // array of colors 
GLuint*		indices_H		= NULL; // array of indices for tatal vector field

GLfloat*	vertices_BOX	= NULL; // array of vertexes for tatal vector field 
GLfloat*	normals_BOX		= NULL; // array of normals for tatal vector field
GLfloat*	colors_BOX		= NULL; // array of colors 
GLuint*		indices_BOX		= NULL; // array of indices for tatal vector field

GLfloat*	vertices_BASIS	= NULL; // array of vertexes for tatal vector field 
GLfloat*	normals_BASIS	= NULL; // array of normals for tatal vector field
GLfloat*	colors_BASIS	= NULL; // array of colors 
GLuint*		indices_BASIS	= NULL; // array of indices for tatal vector field

GLfloat*	vertices_PBC_A	= NULL; // array of vertexes for tatal vector field 
GLfloat*	normals_PBC_A	= NULL; // array of normals for tatal vector field
GLfloat*	colors_PBC_A	= NULL; // array of colors 
GLuint*		indices_PBC_A	= NULL; // array of indices for tatal vector field

GLfloat*	vertices_PBC_B	= NULL; // array of vertexes for tatal vector field 
GLfloat*	normals_PBC_B	= NULL; // array of normals for tatal vector field
GLfloat*	colors_PBC_B	= NULL; // array of colors 
GLuint*		indices_PBC_B	= NULL; // array of indices for tatal vector field

GLfloat*	vertices_PBC_C	= NULL; // array of vertexes for tatal vector field 
GLfloat*	normals_PBC_C	= NULL; // array of normals for tatal vector field
GLfloat*	colors_PBC_C	= NULL; // array of colors 
GLuint*		indices_PBC_C	= NULL; // array of indices for tatal vector field

int			arrowFaces		= 6; // number of arrow faces, default number
int			arrowFaces_H	= 30; // number of arrow faces for applied field vector

GLuint		vboIdV;   // ID of VBO for vertex arrays
GLuint		vboIdN;   // ID of VBO for normal arrays
GLuint		vboIdC;   // ID of VBO for color arrays
GLuint		iboIdI;   // ID of IBO for index arrays

GLuint		vboIdV_H;   // ID of VBO for vertex arrays
GLuint		vboIdN_H;   // ID of VBO for normal arrays
GLuint		vboIdC_H;   // ID of VBO for color arrays
GLuint		iboIdI_H;   // ID of IBO for index arrays

GLuint		vboIdV_BOX;   // ID of VBO for vertex arrays
GLuint		vboIdN_BOX;   // ID of VBO for normal arrays
GLuint		vboIdC_BOX;   // ID of VBO for color arrays
GLuint		iboIdI_BOX;   // ID of IBO for index arrays

GLuint		vboIdV_BASIS;   // ID of VBO for vertex arrays
GLuint		vboIdN_BASIS;   // ID of VBO for normal arrays
GLuint		vboIdC_BASIS;   // ID of VBO for color arrays
GLuint		iboIdI_BASIS;   // ID of IBO for index arrays

GLuint		vboIdV_PBC_A, vboIdV_PBC_B, vboIdV_PBC_C;   // ID of VBO for vertex arrays
GLuint		vboIdN_PBC_A, vboIdN_PBC_B, vboIdN_PBC_C;   // ID of VBO for normal arrays
GLuint		vboIdC_PBC_A, vboIdC_PBC_B, vboIdC_PBC_C;   // ID of VBO for color arrays
GLuint		iboIdI_PBC_A, iboIdI_PBC_B, iboIdI_PBC_C;   // ID of IBO for index arrays

int			ElNumProto;   // number of triangles per arrow
int			IdNumProto;   // number of indixes per arrow
int			VCNumProto;   // number of Vertex and normals Component per arrow

int			ElNum; // total number of triangles for the whole vector field = ElNumProto * Number of spins
int			IdNum; // total number of indixes for the whole vector field = IdNumProto * Number of spins
int			VCNum; // total number of component of vertices for whole vector field = VCNumProto * Number of spins

int			IdNum_H;
int			VCNum_H;

int			IdNum_BOX;
int			VCNum_BOX;

int			IdNum_BASIS;
int			VCNum_BASIS;

int			IdNum_PBC;
int			VCNum_PBC;

int			Play=0;

int			DataTransfer=1;

void			ChangeVectorMode( int );
void			ChangeColorMap( int );
void			ChangeInitialState( int );
void			Buttons( int );
void			Keyboard( unsigned char, int, int );
void			KeyboardAdd( int, int, int );
void			MouseButton( int, int, int, int );
void			MouseMotion( int, int );
float			ElapsedSeconds( );
// color control functions
void			HSVtoRGB(float[3], float [3] , int, int);
void			InitRGB(float* , float* , float* , int*);
// VBO array preparing functions
void			ReallocateArrayDrawing();
void			UpdatePrototypeVerNorInd(float*,float*,GLuint*,int,int,int);
void			CreateNewVBO();
void			UpdateVBO(GLuint * , GLuint * , GLuint * , GLuint * , float * , float * , float * , GLuint * );
void			UpdateVBO_H(GLuint * , GLuint * , GLuint * , GLuint * , float * , float * , float * , GLuint * );
void			UpdateSpinComponents(float * , float * , float * , int);
void			UpdateSpinPositions(float[][3], int[3], float[][3], int, float[3][3], float*, float*, float*);
void			InitSpinComponents(float * , float * , float * , double * , double * , double * , int N);
void			UpdateIndices(GLuint * , int, GLuint *, int, int);
void			UpdateVerticesNormalsColors(float *, float *, int, float *, float *, float *, int, float * , float * , float *, double * , double * , double *, int);
void 			UpdateVerticesNormalsColors_H(float *, float *, int Kinp, float *, float *, float *, float, float, float, float, float, float);
// drawing functions
void			GetBox(float[][3], int[3], float[4][3]);
void			drawVBO();
void			drawVBO_H();
void			drawVBO_BOX();
void 			drawVBO_BASIS();
void			drawVBO_PBC_A();
void			drawVBO_PBC_B();
void			drawVBO_PBC_C();
void			idle();
void			setupTweakBar();
// return the number of seconds since the start of the program:
float ElapsedSeconds( )	{
	int ms = glutGet( GLUT_ELAPSED_TIME );	// get # of milliseconds since the start of the program
	return (float)ms / 1000.;				// convert it to seconds:
}

void Resize( int window_width, int window_height) // called when user resizes the window
{
	asp_rat     = (float)((   ((double)window_width)/((double)window_height)   ));
	asp_rat_inv = (float)((   ((double)window_height)/((double)window_width)   ));

	glutSetWindow( GLUT_window );
	glViewport(0, 0, window_width, window_height);
	    // Send the new window size to AntTweakBar
    TwWindowSize(window_width, window_height);
}

void Xup( )
{
	Rot[1] = 0.;
	Rot[0] = Rot[2] = 270.;
	TransXYZ[0] = TransXYZ[1] = 0.;
}

void Yup( )
{
	Rot[2] = 180;
	Rot[0] = 270.;
	TransXYZ[0] = TransXYZ[1] = 0.;
}

void Zup( )
{
	Rot[0] = Rot[1] = Rot[2] = TransXYZ[0] = TransXYZ[1] = 0.;
}

// use glut to display a string of characters using a raster font:
void
DoRasterString( float x, float y, float z, char const *s )
{
	char const *c;			// one character to print

	glRasterPos3f( (GLfloat)x, (GLfloat)y, (GLfloat)z );
	for( c=s; *c; c++ )
	{
		glutBitmapCharacter( GLUT_BITMAP_8_BY_13, *c );
		//GLUT_BITMAP_8_BY_13
		//GLUT_BITMAP_9_BY_15
		//GLUT_BITMAP_TIMES_ROMAN_10
		//GLUT_BITMAP_TIMES_ROMAN_24
		//GLUT_BITMAP_HELVETICA_10
		//GLUT_BITMAP_HELVETICA_12
	}
}

// use glut to display a string of characters using a stroke font:
void
DoStrokeString( float x, float y, float z, float ht, char const *s )
{
	char const *c;			// one character to print
	glLineWidth(1.0f);
	glPushMatrix( );
		glTranslatef( (GLfloat)x, (GLfloat)y, (GLfloat)z );
		float sf = ht /( 119.05 + 33.33 );
		glScalef( (GLfloat)sf, (GLfloat)sf, (GLfloat)sf );
		for( c=s; *c; c++ )
		{
			glutStrokeCharacter( GLUT_STROKE_ROMAN, *c );
		}
	glPopMatrix( );
}

void GetBox(float abc[][3], int uABC[3], float box[3][3])
{
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

void ChangeBoxSize(int Na, int Nb, int Nc){
	if (Play==1) Play=0;
}


void Display (void)
{
	GLdouble Hight;
	float Vtemp[3];

	glutSetWindow( GLUT_window );// set which window we want to do the graphics into
	glClearColor( Colors[WhichColor][0], Colors[WhichColor][1], Colors[WhichColor][2], 0. );// setup the clear values
	glDrawBuffer( GL_BACK );// erase the background
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_ACCUM_BUFFER_BIT);
	glMatrixMode( GL_PROJECTION ); glLoadIdentity( );

	if( WhichProjection == PERSP )
	gluPerspective(PerspSet[0], asp_rat, PerspSet[2], PerspSet[3]);	
	else{
		Vtemp[0] = CameraEye[0] - CameraC[0];
		Vtemp[1] = CameraEye[1] - CameraC[0];
		Vtemp[2] = CameraEye[2] - CameraC[0]+TransXYZ[2];
		Hight = Unitf(Vtemp,Vtemp); //distance to the point which camera looks at
		Hight = -Hight*tan(D2R*PerspSet[0]/2.f); // hight of the view frame at the plane perp to cam view line
		glOrtho( Hight*asp_rat, -Hight*asp_rat, Hight,  -Hight,  PerspSet[2], PerspSet[3] );
	//         (     left,          right,      bottom,   top,       near,        far     )
	}
	// place the objects into the scene:
	glMatrixMode( GL_MODELVIEW ); glLoadIdentity( );
	// set the eye position, look-at position, and up-vector:
	gluLookAt(	CameraEye[0],	CameraEye[1],	CameraEye[2],   // position  
				CameraC[0],		CameraC[1],		CameraC[2],     // look at
				CameraUp[0],	CameraUp[1],	CameraUp[2]);   // up

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
	float v[4]; 
    v[0] = v[1] = v[2] = g_LightMultiplier*0.4f; v[3] = 1.0f;
    glLightfv(GL_LIGHT0, GL_AMBIENT, v);
    v[0] = v[1] = v[2] = g_LightMultiplier*0.8f; v[3] = 1.0f;
    glLightfv(GL_LIGHT0, GL_DIFFUSE, v);
    v[0] = -g_LightDirection[0]; v[1] = -g_LightDirection[1]; v[2] = -g_LightDirection[2]; v[3] = 0.0f;
    glLightfv(GL_LIGHT0, GL_POSITION, v);

	// translate the scene:
	glTranslatef( (GLfloat)TransXYZ[0], (GLfloat)TransXYZ[1], -(GLfloat)TransXYZ[2] );
	// rotate the scene:
	glRotatef( (GLfloat)Rot[0], 1., 0., 0. ); //@X
	glRotatef( (GLfloat)Rot[1], 0., 1., 0. ); //@Y
	glRotatef( (GLfloat)Rot[2], 0., 0., 1. ); //@Z

	drawVBO(); // Draw VBO for spins
	drawVBO_H(); // Draw VBO for vector representing the firld direction 
	
	// possibly draw the box and periodic boundary condition :
	if( BoxOn != 0 ) {	
		drawVBO_BOX();
		if(Boundary[0]!=0) drawVBO_PBC_A();;
		if(Boundary[1]!=0) drawVBO_PBC_B();
		if(Boundary[2]!=0) drawVBO_PBC_C();
	}
	// possibly draw the axes:
	if( AxesOn != 0 ) drawVBO_BASIS();//glCallList( AxesList );//
    // Draw tweak bars
    TwDraw();
    // Present frame buffer
    glutSwapBuffers();
    // Recall Display at next frame
    glutPostRedisplay();
}

void setupOpenGL ()
{
	//  Connect to the windowing system + create a window
	//  with the specified dimensions and position
	//  + set the display mode + specify the window title.
	int glutArgc = 0;
	glutInit(&glutArgc, NULL); //glutInit(&argc, argv);
	glutInitWindowSize (window_width, window_height);
	glutInitWindowPosition (0, 0);
    glutInitDisplayMode( GLUT_RGB|GLUT_DOUBLE|GLUT_DEPTH|GLUT_MULTISAMPLE|GLUT_ALPHA);
	GLUT_window = glutCreateWindow (WINDOWTITLE);
	#if !defined(__APPLE__)
	glewInit();
	#endif

	for (int i=0;i<NumCamPosSave;i++){
		CameraPosition[i][0]=0;
		CameraPosition[i][1]=0;
		CameraPosition[i][2]=0;
		CameraPosition[i][3]=0;
		CameraPosition[i][4]=0;
		CameraPosition[i][5]=0;
		CameraPosition[i][6]=PerspSet[0];
	}

	InitRGB(RHue, GHue, BHue, HueMapRGB);
	// Set the GLUT callback functions
	glutDisplayFunc( Display );// DisplayFunc -- redraws the OpenGl main window
	glutReshapeFunc( Resize );// ReshapeFunc -- handles the user resizing the window
	glutKeyboardFunc( Keyboard );// KeyboardFunc -- handles a keyboard input
	glutSpecialFunc( KeyboardAdd );// SpecialFunc -- handle special keys on the keyboard e.g. F1, Num+...
	glutMouseFunc( MouseButton );// MouseFunc -- handle the mouse button going down or up
	glutMotionFunc( MouseMotion );// MotionFunc -- handle the mouse moving with a button down
	glutPassiveMotionFunc((GLUTmousemotionfun)TwEventMouseMotionGLUT);
	// glutKeyboardFunc((GLUTkeyboardfun)TwEventKeyboardGLUT);
	// glutMouseFunc((GLUTmousebuttonfun)TwEventMouseButtonGLUT);
	// glutMotionFunc((GLUTmousemotionfun)TwEventMouseMotionGLUT);
	// glutSpecialFunc((GLUTspecialfun)TwEventSpecialGLUT);
	// glutVisibilityFunc( Visibility );// VisibilityFunc -- handle a change in window visibility
	// glutPassiveMotionFunc( NULL );// PassiveMotionFunc -- handle the mouse moving with a button up
	// glutEntryFunc( NULL );// EntryFunc	-- handle the cursor entering or leaving the window
	// glutSpaceballMotionFunc( NULL );// SpaceballMotionFunc -- handle spaceball translation
	// glutSpaceballRotateFunc( NULL );// SpaceballRotateFunc -- handle spaceball rotation
	// glutSpaceballButtonFunc( NULL );// SpaceballButtonFunc -- handle spaceball button hits
	// glutButtonBoxFunc( NULL );// ButtonBoxFunc -- handle button box hits
	// glutDialsFunc( NULL );// DialsFunc -- handle dial rotations
	// glutTabletMotionFunc( NULL );// TabletMotionFunc -- handle digitizing tablet motion
	// glutTabletButtonFunc( NULL );// TabletButtonFunc -- handle digitizing tablet button hits
	// glutMenuStateFunc( NULL );// MenuStateFunc -- declare when a pop-up menu is in use
	// glutTimerFunc( 0, NULL, 0 );// TimerFunc -- trigger something to happen a certain time from now
	glutIdleFunc( idle );// IdleFunc -- what to do when nothing else is going on

	// Initialize AntTweakBar
    TwInit(TW_OPENGL, NULL);
    // - Send 'glutGetModifers' function pointer to AntTweakBar;
    //   required because the GLUT key event functions do not report key modifiers states.
    TwGLUTModifiersFunc(glutGetModifiers);
	setupTweakBar();

	glEnable(GL_DEPTH_TEST);

	glShadeModel (GL_SMOOTH); //set the shader to smooth shader GL_SMOOTH/GL_FLAT

	glEnable(GL_COLOR_MATERIAL);
	glCullFace(GL_FRONT);//GL_FRONT//GL_FRONT_AND_BACK
	glEnable(GL_CULL_FACE);

	glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);

	glEnable(GL_LINE_SMOOTH);               
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

	glEnable(GL_POLYGON_SMOOTH);
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
}

void idle ()
{  
	currentTime = glutGet(GLUT_ELAPSED_TIME);
	timeInterval = currentTime - previousTime;

	if(timeInterval > 40 && Play==1)//40ms gives approximately 25 FPS +/-1 if the engine works faster then 25 IPS
	{
		if( DATA_TRANSFER_MUTEX==TAKE_DATA )
		{
			if (DataTransfer==1)
			//if (Play==1)
			{
			ChangeVectorMode(1);
			totalEnergy = GetTotalEnergy( 	bSx, bSy, bSz, 
						NeighborPairs, AIdxBlock, NIdxBlock, NIdxGridA, NIdxGridB, NIdxGridC, SIdx,
						Jij, Bij, Dij, VDMx, VDMy, VDMz, VKu, Ku, Kc, VHf, Hf, Etot, Mtot, NOS );
			mtot[0] = Mtot[0]/NOS;
			mtot[1] = Mtot[1]/NOS;
			mtot[2] = Mtot[2]/NOS;
			perSpEnergy = totalEnergy/NOS;
			totalEnergyFerro = GetTotalEnergyFerro( VHf[0], VHf[1], VHf[2], 
						NeighborPairs, AIdxBlock, NIdxBlock, NIdxGridA, NIdxGridB, NIdxGridC, SIdx,
						Jij, Bij, Dij, VDMx, VDMy, VDMz, VKu, Ku, Kc, VHf, Hf, Etot, NOS );
			totalEnergyFerro = totalEnergyFerro/NOS;	
			perSpEnergyMinusFerro = perSpEnergy - totalEnergyFerro;
			}
			if (timeInterval > 1000)//~2 seconds
			{
				FPS = frameCount / (timeInterval * 0.002f);
				previousTime = currentTime;
				frameCount = 0;	
				IPS = (currentIteration - previousIteration)/ (timeInterval * 0.002f);;
				previousIteration = currentIteration;
			}
			pthread_mutex_lock( &show_mutex);
				DATA_TRANSFER_MUTEX = WAIT_DATA; // meaning that OpenGL is waiting for new data from engine
			pthread_mutex_unlock( &show_mutex);
		} 

		glutSetWindow (GLUT_window);
		glutPostRedisplay ();
		frameCount++;
	}
}

void TW_CALL CB_Run( void *clientData )
{
	if (Play==0)
	{
		Play=1;
		TwDefine(" Parameters&Controls/Run  label='STOP simulation' ");
		pthread_mutex_lock(&culc_mutex);
			ENGINE_MUTEX=DO_IT;
		pthread_mutex_unlock(&culc_mutex);
	}else{
		Play=0; 
		TwDefine(" Parameters&Controls/Run  label='RUN simulation' ");
		pthread_mutex_lock(&culc_mutex);
			ENGINE_MUTEX=WAIT;
		pthread_mutex_unlock(&culc_mutex);	
	}
}

void TW_CALL CB_SetScale(const void *value, void *clientData )
{
	(void)clientData; // unused
    Scale = *( float *)value; // copy value to Scale
// if(WhichVectorMode==ARROW1||WhichVectorMode==CONE1) 
//    {UpdatePrototypeVerNorInd(vertexProto, normalProto, indicesProto, arrowFaces, WhichVectorMode,0);}
	ChangeVectorMode (0);
}

//  Callback function called by the tweak bar to get the 'AutoRotate' value
void TW_CALL CB_GetScale(void *value, void *clientData)
{
    (void)clientData; // unused
    *(float *)value = Scale; // just copy Scale to value
}


void TW_CALL CB_SetVectorMode(const void *value, void *clientData )
{
	(void)clientData; // unused
    WhichVectorMode = *( enVectorMode *)value; // copy value to WhichVectorMode
    ChangeVectorMode (0);
}


void TW_CALL CB_GetVectorMode(void *value, void *clientData)
{
    (void)clientData; // unused
    *(int *)value = WhichVectorMode; // just copy WhichVectorMode to value
}


void TW_CALL CB_SetFaces(const void *value, void *clientData )
{
	(void)clientData; // unused
    arrowFaces = *( int *)value; // copy value to arrowFaces
	ChangeVectorMode (0);
}


void TW_CALL CB_GetFaces(void *value, void *clientData)
{
    (void)clientData; // unused
    *(int *)value = arrowFaces; // just copy arrowFaces to value
}

void TW_CALL CB_SetPivot(const void *value, void *clientData )
{
	(void)clientData; // unused
    Pivot = *( float *)value; // copy value to Pivot
	ChangeVectorMode (0);
}


void TW_CALL CB_GetPivot(void *value, void *clientData)
{
    (void)clientData; // unused
    *(float *)value = Pivot; // just copy Pivot to value
}

void TW_CALL CB_SaveCameraPosition ( void *clientData )
{
	CameraPosition[CurrentCameraPositionBank][0]=Rot[0];
	CameraPosition[CurrentCameraPositionBank][1]=Rot[1];
	CameraPosition[CurrentCameraPositionBank][2]=Rot[2];
	CameraPosition[CurrentCameraPositionBank][3]=TransXYZ[0];
	CameraPosition[CurrentCameraPositionBank][4]=TransXYZ[1];
	CameraPosition[CurrentCameraPositionBank][5]=TransXYZ[2];
	CameraPosition[CurrentCameraPositionBank][6]=PerspSet[0];
}

void TW_CALL CB_GetCameraPosition ( void *clientData )
{
	Rot[0]=CameraPosition[CurrentCameraPositionBank][0];
	Rot[1]=CameraPosition[CurrentCameraPositionBank][1];
	Rot[2]=CameraPosition[CurrentCameraPositionBank][2];
	TransXYZ[0]=CameraPosition[CurrentCameraPositionBank][3];
	TransXYZ[1]=CameraPosition[CurrentCameraPositionBank][4];
	TransXYZ[2]=CameraPosition[CurrentCameraPositionBank][5];
	PerspSet[0]=CameraPosition[CurrentCameraPositionBank][6];
}

void TW_CALL CB_SetColorShift(const void *value, void *clientData )
{
	(void)clientData; // unused
    ColorShift = *( int *)value; // copy value to ColorShift
		if(WhichColorScheme==RGB)
		{
		HueMap[0]=HueMapRGB[0]+ColorShift;;
		HueMap[1]=HueMapRGB[1]+ColorShift;;
		HueMap[2]=HueMapRGB[2]+ColorShift;;
		HueMap[3]=HueMapRGB[3]+ColorShift;;
		HueMap[4]=HueMapRGB[4]+ColorShift;;
		HueMap[5]=HueMapRGB[5]+ColorShift;;
		} 	
		else
		{
		HueMap[0]=HueMapRYGB[0]+ColorShift;;
		HueMap[1]=HueMapRYGB[1]+ColorShift;;
		HueMap[2]=HueMapRYGB[2]+ColorShift;;
		HueMap[3]=HueMapRYGB[3]+ColorShift;;
		HueMap[4]=HueMapRYGB[4]+ColorShift;;
		HueMap[5]=HueMapRYGB[5]+ColorShift;;
		}
		InitRGB(RHue, GHue, BHue, HueMap);
		ChangeVectorMode (1);
}

void TW_CALL CB_GetColorShift(void *value, void *clientData)
{
    (void)clientData; // unused
    *(int *)value = ColorShift; // just copy ColorShift to value
}


void TW_CALL CB_SetColorScheme(const void *value, void *clientData )
{
	(void)clientData; // unused
    WhichColorScheme = *( enColorScheme *)value; // copy value to WhichColorScheme
		if(WhichColorScheme==RGB)
		{
		HueMap[0]=HueMapRGB[0]+ColorShift;;
		HueMap[1]=HueMapRGB[1]+ColorShift;;
		HueMap[2]=HueMapRGB[2]+ColorShift;;
		HueMap[3]=HueMapRGB[3]+ColorShift;;
		HueMap[4]=HueMapRGB[4]+ColorShift;;
		HueMap[5]=HueMapRGB[5]+ColorShift;;
		} 	
		else
		{
		HueMap[0]=HueMapRYGB[0]+ColorShift;;
		HueMap[1]=HueMapRYGB[1]+ColorShift;;
		HueMap[2]=HueMapRYGB[2]+ColorShift;;
		HueMap[3]=HueMapRYGB[3]+ColorShift;;
		HueMap[4]=HueMapRYGB[4]+ColorShift;;
		HueMap[5]=HueMapRYGB[5]+ColorShift;;
		}
		InitRGB(RHue, GHue, BHue, HueMap);
		ChangeVectorMode (1);
}


void TW_CALL CB_GetColorScheme(void *value, void *clientData)
{
    (void)clientData; // unused
    *(enColorScheme *)value = WhichColorScheme; // just copy WhichColorScheme to value
}


void TW_CALL CB_SetInvHue(const void *value, void *clientData )
{
	(void)clientData; // unused
    InvertHue = *( int *)value; // copy value to InvertHue
	ChangeVectorMode (1);
}

void TW_CALL CB_GetInvHue(void *value, void *clientData)
{
    (void)clientData; // unused
    *(int *)value = InvertHue; // just copy InvertHue to value
}

void TW_CALL CB_SetInvVal(const void *value, void *clientData )
{
	(void)clientData; // unused
    InvertValue = *( int *)value; // copy value to InvertValue
	ChangeVectorMode (1);
}

void TW_CALL CB_GetInvVal(void *value, void *clientData)
{
    (void)clientData; // unused
    *(int *)value = InvertValue; // just copy InvertValue to value
}


void TW_CALL CB_SetHfield(const void *value, void *clientData )
{
	(void)clientData; // unused
    Hf = *( float *)value; // copy value to InvertValue
    if (Hf>1.1*fabs(Jij[0])) Hf=1.1*fabs(Jij[0]);
	UpdateVerticesNormalsColors_H(vertexProto_H, normalProto_H, VCNum_H, vertices_H, normals_H, colors_H, Box[0][0]*0.6, Box[1][1]*0.6, Box[2][2]*0.6, VHf[0], VHf[1], VHf[2]);
	UpdateVBO_H(&vboIdV_H, &vboIdN_H, &vboIdC_H, &iboIdI_H, vertices_H, normals_H, colors_H, indices_H);
}

void TW_CALL CB_GetHfield(void *value, void *clientData)
{
    (void)clientData; // unused
    *(float *)value = Hf; // just copy InvertValue to value
}

// void TW_CALL SetCallback(const void *value, void *clientData)
// { 
//     myVariable = *(const MyVariableType *)value;  // for instance
// }

void TW_CALL CB_SetHfieldTheta(const void *value, void *clientData )
{
	(void)clientData; // unused
    VHtheta = *(float*)value;
    VHf[0]=sin(PI*VHtheta/180)*cos(PI*VHphi/180);
	VHf[1]=sin(PI*VHtheta/180)*sin(PI*VHphi/180);
	VHf[2]=cos(PI*VHtheta/180);
	UpdateVerticesNormalsColors_H(vertexProto_H, normalProto_H, VCNum_H, vertices_H, normals_H, colors_H, Box[0][0]*0.6, Box[1][1]*0.6, Box[2][2]*0.6, VHf[0], VHf[1], VHf[2]);
	UpdateVBO_H(&vboIdV_H, &vboIdN_H, &vboIdC_H, &iboIdI_H, vertices_H, normals_H, colors_H, indices_H);
}

void TW_CALL CB_GetHfieldTheta(void *value, void *clientData)
{
    (void)clientData; 
    *(float*)value = VHtheta; 
}

void TW_CALL CB_SetHfieldPhi(const void *value, void *clientData )
{
	(void)clientData; // unused
    VHphi = *(float*)value;
    VHf[0]=sin(PI*VHtheta/180)*cos(PI*VHphi/180);
	VHf[1]=sin(PI*VHtheta/180)*sin(PI*VHphi/180);
	VHf[2]=cos(PI*VHtheta/180);
	UpdateVerticesNormalsColors_H(vertexProto_H, normalProto_H, VCNum_H, vertices_H, normals_H, colors_H, Box[0][0]*0.6, Box[1][1]*0.6, Box[2][2]*0.6, VHf[0], VHf[1], VHf[2]);
	UpdateVBO_H(&vboIdV_H, &vboIdN_H, &vboIdC_H, &iboIdI_H, vertices_H, normals_H, colors_H, indices_H);
}

void TW_CALL CB_GetHfieldPhi(void *value, void *clientData)
{
    (void)clientData; 
    *(float*)value = VHphi; 
}

void TW_CALL CB_SetHfieldXYZ(const void *value, void *clientData )
{
	(void)clientData; // unused
 //    VHtheta = *(float*)value;
 //    VHf[0]=sin(PI*VHtheta/180)*cos(PI*VHphi/180);
	// VHf[1]=sin(PI*VHtheta/180)*sin(PI*VHphi/180);
	// VHf[2]=cos(PI*VHtheta/180);
	UpdateVerticesNormalsColors_H(vertexProto_H, normalProto_H, VCNum_H, vertices_H, normals_H, colors_H, Box[0][0]*0.6, Box[1][1]*0.6, Box[2][2]*0.6, VHf[0], VHf[1], VHf[2]);
	UpdateVBO_H(&vboIdV_H, &vboIdN_H, &vboIdC_H, &iboIdI_H, vertices_H, normals_H, colors_H, indices_H);
}

void TW_CALL CB_GetHfieldXYZ(void *value, void *clientData)
{
    (void)clientData; 
    // *(float*)value = VHtheta; 
}



void TW_CALL CB_SetACPeriod(const void *value, void *clientData )
{
    (void)clientData; // unused
    Period_dc = *( float *)value; // copy value to Period_dc
    Omega_dc = TPI/Period_dc;
}

void TW_CALL CB_GetACPeriod(void *value, void *clientData)
{
    (void)clientData; // unused
    *(float *)value = Period_dc; // just copy Period_dc to value
}

void TW_CALL CB_SetOmega(const void *value, void *clientData )
{
    (void)clientData; // unused
    Omega_dc = *( float *)value; // copy value to Period_dc
    Period_dc = TPI/Omega_dc;
}

void TW_CALL CB_GetOmega(void *value, void *clientData)
{
    (void)clientData; // unused
    *(float *)value = Omega_dc; // just copy Omega_dc to value
}

void TW_CALL CB_SetInitial( void *clientData )
{
	ChangeInitialState( WhichInitialState );
}

void TW_CALL CB_RotateAllSpins( void *clientData )
{	
	if(fabs(chDir[0])+fabs(chDir[1])+fabs(chDir[2])!=0)
	{
		double tmp[3];
		for (int i=0; i<NOS; i++) {
			RotateVector(Sx[i], Sy[i], Sz[i], chDir[0], chDir[1], chDir[2], RotateAllSpins, tmp);
			bSx[i] = Sx[i] = tmp[0];
			bSy[i] = Sy[i] = tmp[1];
			bSz[i] = Sz[i] = tmp[2];
		}
		ChangeVectorMode(1);		
	}
}

void TW_CALL CB_InvertX( void *clientData )
{
	for (int i=0; i<NOS; i++) {Sx[i] = -Sx[i]; bSx[i] = -bSx[i];}
	ChangeVectorMode(1);
}

void TW_CALL CB_InvertY( void *clientData )
{
	for (int i=0; i<NOS; i++) {Sy[i] = -Sy[i]; bSy[i] = -bSy[i];}
	ChangeVectorMode(1);
}

void TW_CALL CB_InvertZ( void *clientData ){
	for (int i=0; i<NOS; i++) {Sz[i] = -Sz[i]; bSz[i] = -bSz[i];}
	ChangeVectorMode(1);
}

void TW_CALL CB_CleanSxSySzFile( void *clientData ){
  fclose (outFile);//outFile is a global variable - pointer FILE* see also CALC_THREAD in ENGINE.cpp
	outFile = fopen ("table.csv","w");
	if (outFile!=NULL) {fputs ("iter,time,Mx,My,Mz,E_tot,\n",outFile);}
 }


void TW_CALL CB_SetSliceMode(const void *value, void *clientData ){
	(void)clientData; // unused
    WhichSliceMode = *( enSliceMode *)value; // copy value to WhichSliceMode
    ChangeVectorMode(0);
}


void TW_CALL CB_GetSliceMode(void *value, void *clientData){
    (void)clientData; // unused
    *(int *)value = WhichSliceMode; // just copy WhichSliceMode to value
}


void TW_CALL CB_GetThetaMax(void *value, void *clientData){
    (void)clientData; // unused
    *(float *)value = theta_max1; 
}


void TW_CALL CB_SetThetaMax(const void *value, void *clientData ){
	(void)clientData; // unused
	float test= *( float *)value; 
	if (test>=theta_min1 ){
        theta_max1 = test;
        Sz_min1=cos(theta_max1);
        ChangeVectorMode(0);
	}
}


void TW_CALL CB_GetThetaMin(void *value, void *clientData){
    (void)clientData; // unused
    *(float *)value = theta_min1; 
}

void TW_CALL CB_SetThetaMin(const void *value, void *clientData ){
	(void)clientData; // unused
	float test= *( float *)value; 
	if (test<=theta_max1 ){
        theta_min1 = test; 
        Sz_max1=cos(theta_min1);
        ChangeVectorMode(0);
	}
}

void TW_CALL CB_GetPhiMax(void *value, void *clientData){
    (void)clientData; // unused
    *(int *)value = phi_max1; 
}


void TW_CALL CB_SetPhiMax(const void *value, void *clientData ){
	(void)clientData; // unused
	int test= *( int *)value; 
	if (test>=phi_min1 ){
        phi_max1 = test;
        ChangeVectorMode(0);
	}
}

void TW_CALL CB_GetPhiMin(void *value, void *clientData){
    (void)clientData; // unused
    *(int *)value = phi_min1; 
}


void TW_CALL CB_SetPhiMin(const void *value, void *clientData ){
	(void)clientData; // unused
	int test= *( int *)value; 
	if (test<=phi_max1 ){
        phi_min1 = test;
        ChangeVectorMode(0);
	}
}

void TW_CALL CB_ResetIterations( void *clientData ){
  ITERATION=0;
  //glutPostRedisplay();
}


void TW_CALL CB_SaveCSV( void *clientData )
{
FILE * pFile;
  pFile = fopen (outputfilename,"w");
  if (pFile!=NULL)
  {
  	
	 //  	fputs ("px,py,pz,nx,ny,nz,\n",pFile);
	 //  	for (int i=0;i<NOS;i++)
	 //  	{
		// snprintf(shortBufer,80,"%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,\n",Px[i],Py[i],Pz[i],bSx[i],bSy[i],bSz[i]);
	 //    fputs (shortBufer,pFile);  		
	 //  	}
	 //    fclose (pFile);

		int an,bn,cn,atom,n,N;
		int anini=0;
		int anfin=uABC[0];
		int bnini=0;
		int bnfin=uABC[1];
		int cnini=0;
		int cnfin=uABC[2];
		if (save_slice==1){
			switch( WhichSliceMode){
				case A_AXIS:
					anini=A_layer_min-1;
			        anfin=A_layer_max;
				break;
				case B_AXIS:
					bnini=B_layer_min-1;
			        bnfin=B_layer_max;
				break;
				case C_AXIS:
					cnini=C_layer_min-1;
			        cnfin=C_layer_max;
				break;
				default:
				break;
			}
		}

		if (WhichAverageMode==ALONG_0){
			printf("Write to csv-file!\n");
		  	fputs ("px,py,pz,nx,ny,nz,\n",pFile);
		  	for (cn = cnini; cn<cnfin; cn++) {
				for (bn = bnini; bn<bnfin; bn++) {
					for (an = anini; an<anfin; an++) {
						n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
						n = n*AtomsPerBlock;//index of the first spin in the block
						for (atom=0; atom<AtomsPerBlock; atom++){
						    N = n + atom;
							snprintf(shortBufer,80,"%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,\n",Px[N],Py[N],Pz[N],bSx[N],bSy[N],bSz[N]);
							fputs (shortBufer,pFile); 
						}
					}
				}
			}
		    fclose (pFile);
		}else{// WhichAverageMode==ALONG_A or ALONG_B or ALONG_C
			double tSpin[3];
			float tPositin[3]; 
			fputs ("px,py,pz,<nx>,<ny>,<nz>,|<n>|\n",pFile);
			if (WhichAverageMode==ALONG_C){
				for (bn = bnini; bn<bnfin; bn++){
					for (an = anini; an<anfin; an++){
						tSpin[0]=0;
						tSpin[1]=0;
						tSpin[2]=0;
						tPositin[0]=abc[0][0]*0.5+abc[1][0]*0.5+abc[2][0]*0.5;
						tPositin[1]=abc[0][1]*0.5+abc[1][1]*0.5+abc[2][1]*0.5;
						tPositin[2]=abc[0][2]*0.5+abc[1][2]*0.5+abc[2][2]*0.5;
						for (cn = cnini; cn<cnfin; cn++){
							n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
							n = n*AtomsPerBlock;//index of the first spin in the block
							for (atom=0; atom<AtomsPerBlock; atom++){
							    N = n + atom;
							    tSpin[0]+=bSx[N];
							    tSpin[1]+=bSy[N];
							    tSpin[2]+=bSz[N];
							}
						}
						tPositin[0]+=abc[0][0]*an+abc[1][0]*bn;
						tPositin[1]+=abc[0][1]*an+abc[1][1]*bn;
						tPositin[2]+=abc[0][2]*an+abc[1][2]*bn;
						double modulus=sqrt(tSpin[0]*tSpin[0] + tSpin[1]*tSpin[1] + tSpin[2]*tSpin[2]);
						int cN = cnfin - cnini;
						snprintf(shortBufer,80,"%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,\n",tPositin[0],tPositin[1],tPositin[2],tSpin[0]/cN,tSpin[1]/cN,tSpin[2]/cN,modulus/cN);
						fputs (shortBufer,pFile);
					}
				}
			}else if(WhichAverageMode==ALONG_B){
				for (an = anini; an<anfin; an++){
					for (cn = cnini; cn<cnfin; cn++){
						tSpin[0]=0;
						tSpin[1]=0;
						tSpin[2]=0;
						tPositin[0]=abc[0][0]*0.5+abc[1][0]*0.5+abc[2][0]*0.5;
						tPositin[1]=abc[0][1]*0.5+abc[1][1]*0.5+abc[2][1]*0.5;
						tPositin[2]=abc[0][2]*0.5+abc[1][2]*0.5+abc[2][2]*0.5;
						for (bn = bnini; bn<bnfin; bn++){
							n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
							n = n*AtomsPerBlock;//index of the first spin in the block
							for (atom=0; atom<AtomsPerBlock; atom++){
							    N = n + atom;
							    tSpin[0]+=bSx[N];
							    tSpin[1]+=bSy[N];
							    tSpin[2]+=bSz[N];
							}
						}
						tPositin[0]+=abc[0][0]*an+abc[2][0]*cn;
						tPositin[1]+=abc[0][1]*an+abc[2][1]*cn;
						tPositin[2]+=abc[0][2]*an+abc[2][2]*cn;
						double modulus=sqrt(tSpin[0]*tSpin[0] + tSpin[1]*tSpin[1] + tSpin[2]*tSpin[2]);
						int bN = bnfin - bnini;
						snprintf(shortBufer,80,"%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,\n",tPositin[0],tPositin[1],tPositin[2],tSpin[0]/bN,tSpin[1]/bN,tSpin[2]/bN,modulus/bN);
						fputs (shortBufer,pFile);
					}
				}
			}else if(WhichAverageMode==ALONG_A){
				for (cn = cnini; cn<cnfin; cn++){
					for (bn = bnini; bn<bnfin; bn++){
						tSpin[0]=0;
						tSpin[1]=0;
						tSpin[2]=0;
						tPositin[0]=abc[0][0]*0.5+abc[1][0]*0.5+abc[2][0]*0.5;
						tPositin[1]=abc[0][1]*0.5+abc[1][1]*0.5+abc[2][1]*0.5;
						tPositin[2]=abc[0][2]*0.5+abc[1][2]*0.5+abc[2][2]*0.5;
						for (an = anini; an<anfin; an++){
							n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
							n = n*AtomsPerBlock;//index of the first spin in the block
							for (atom=0; atom<AtomsPerBlock; atom++){
							    N = n + atom;
							    tSpin[0]+=bSx[N];
							    tSpin[1]+=bSy[N];
							    tSpin[2]+=bSz[N];
							}
						}
						tPositin[0]+=abc[2][0]*cn+abc[1][0]*bn;
						tPositin[1]+=abc[2][1]*cn+abc[1][1]*bn;
						tPositin[2]+=abc[2][2]*cn+abc[1][2]*bn;
						double modulus=sqrt(tSpin[0]*tSpin[0] + tSpin[1]*tSpin[1] + tSpin[2]*tSpin[2]);
						int aN = anfin - anini;
						snprintf(shortBufer,80,"%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,\n",tPositin[0],tPositin[1],tPositin[2],tSpin[0]/aN,tSpin[1]/aN,tSpin[2]/aN,modulus/aN);
						fputs (shortBufer,pFile);
					}
				}
			printf("averaged over c-->%s done!\n",outputfilename);
			}
		    fclose (pFile);
		}
  }
}

void TW_CALL CB_ReadCSV( void *clientData )
{
    FILE * pFile = fopen(inputfilename, "r");
    if(pFile) {
		char c;
		char line[120];
		float px,py,pz,sx,sy,sz;
		int pos=0;
    	do{ // read titles of the columns 
		    	c = (char)fgetc(pFile);//get char and move pointer to the next position
		    	if (c != EOF) {
		    		 	line[pos++] = c;	
		    	}
		    }while(c != EOF && c != '\n');

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
			//printf("%f,%f,%f,%f,%f,%f\n", px,py,pz,sx,sy,sz);
			if (i<NOS) 
			{
				// Px[i]=px;
				// Py[i]=py;
				// Pz[i]=pz;
				bSx[i]=Sx[i]=sx;
				bSy[i]=Sy[i]=sy;
				bSz[i]=Sz[i]=sz;
		 	}
		 	i++;
		}while(c != EOF); 
    fclose(pFile);           
    }
	ChangeVectorMode(1);
}

int ReadHeaderLine(FILE * fp, char * line){
	char c;//single character 
	int pos=0;
	do{ c = (char)fgetc(fp);//get current char and move pointer to the next position
    	if (c != EOF && c != '\n') { line[pos++] = c;}//if it's not the end of the file
	}while(c != EOF && c != '\n');//if it's not the end of the file or end of the line
	line[pos] = 0;//complite the readed line
	if ((pos==0 || line[0]!='#') && c != EOF){
		return ReadHeaderLine(fp, line);// recursive call for ReadHeaderLine if the current line is empty
	} 
	return pos-1;// the last simbol is the line end simbol
}

void ReadDataLine(FILE * fp, char * line){
	char c;//single character 
	int pos=0;
	do{ c = (char)fgetc(fp);
    	if (c != EOF && c != '\n') { line[pos++] = c;}
	}while(c != EOF && c != '\n');
	line[pos] = 0;
}

void TW_CALL CB_ReadOVF( void *clientData )
{
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
    FILE * FilePointer = fopen(inputfilename, "rb");
	if(FilePointer!=NULL) {	
		lineLength=ReadHeaderLine(FilePointer, line);//read and check the first nonempty line which starts with '#'
		if (lineLength==-1) {// if there are no one line which starts with '#'
			printf("%s has a wrong file format! \n", inputfilename);
		}else{
		    sscanf(line, "# %s %s %s", keyW1, keyW2, keyW3 );
		    if(strncmp(keyW1, "OOMMF",5)!=0 || strncmp(keyW2, "OVF",  3)!=0 || strncmp(keyW3, "2.0",  3)!=0){
		        //if the first line isn't "OOMMF OFV 2.0"
		    	printf("%s has wrong header of wrong file format! \n", inputfilename);
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
			//if (xnodes>uABC[0]) {imax = uABC[0];}else{imax = xnodes;}
			//if (ynodes>uABC[1]) {jmax = uABC[1];}else{jmax = ynodes;}
			//if (znodes>uABC[2]) {kmax = uABC[2];}else{kmax = znodes;}
			int n;
			if (strncmp(keyW2, "Text",4)==0){
				//Text data format
				printf("...reading data in text format: %s \n", inputfilename);
				for (int k=0; k<znodes; k++){
					for (int j=0; j<ynodes; j++){
						for (int i=0; i<xnodes; i++){
							ReadDataLine(FilePointer, line);
							if (k<uABC[2] && j<uABC[1] && i<uABC[0]){
								n = i + j*uABC[0] + k*uABC[0]*uABC[1];
								sscanf(line, "%lf %lf %lf", &bSx[n],&bSy[n],&bSz[n]);
								Sx[n]=bSx[n]; Sy[n]=bSy[n];Sz[n]=bSz[n];
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
				printf("...reading data of binary (%d) format: %s \n", binType, inputfilename);
				// fread (&bSx[0],binType,1,FilePointer);
				// //printf("%f\n",nx[0]);
				// for (int k=0; k<znodes; k++){
				// 	for (int j=0; j<ynodes; j++){
				// 		for (int i=0; i<xnodes; i++){
				// 			int n = i + j*xnodes + k*xnodes*ynodes;
				// 			fread (&bSx[n],binType,1,FilePointer);
				// 			fread (&bSy[n],binType,1,FilePointer);
				// 			fread (&bSz[n],binType,1,FilePointer);
				// 			Sx[n]=bSx[n]; Sy[n]=bSy[n];Sz[n]=bSz[n];
				// 		}
				// 	}
				// }
				if(fread(&bSx[0],binType,1,FilePointer)) {
				//printf("%f \n",bSx[0]);	

					for (int k=0; k<znodes; k++){
						for (int j=0; j<ynodes; j++){
							for (int i=0; i<xnodes; i++){
								if (k<uABC[2] && j<uABC[1] && i<uABC[0]){
									n = i + j*xnodes + k*xnodes*ynodes; //index of the block!
									//printf("n=%d\n", n);
									if (binType==4){
										if(!fread(&temp4_x,binType,1,FilePointer)) break;
										if(!fread(&temp4_y,binType,1,FilePointer)) break;
										if(!fread(&temp4_z,binType,1,FilePointer)) break;
										for (int t=0; t<AtomsPerBlock; t++){
											int I=n*AtomsPerBlock+t;
											Sx[I]=bSx[I]=(double)temp4_x; 
											Sy[I]=bSy[I]=(double)temp4_y;
											Sz[I]=bSz[I]=(double)temp4_z;		
										}
									}else{
										if(!fread(&temp8_x,binType,1,FilePointer)) break;
										if(!fread(&temp8_y,binType,1,FilePointer)) break;
										if(!fread(&temp8_z,binType,1,FilePointer)) break;
										for (int t=0; t<AtomsPerBlock; t++){
											int I=n*AtomsPerBlock+t;
											Sx[I]=bSx[I]=temp8_x; 
											Sy[I]=bSy[I]=temp8_y;
											Sz[I]=bSz[I]=temp8_z;		
										}
									}
								}	
							}
						}
					}
				}else{printf("problem\n");}
			}else{
				printf("Do not know what to do with \"%s\" data format in %s\n", keyW2, inputfilename);
			}
		}else{
			printf("%s has wrong data format or dimentionality!\n", inputfilename);
		}       
		// when everything is done
		printf("Done!\n");
		fclose(FilePointer);
	}else{printf("Cannot open file: %s \n", inputfilename);}
	ChangeVectorMode(1);
}

void readConfigFile()
{
	char  configfilename[64] = "magnoom.cfg";
	char  line[256];//whole line of header should be not longer then 256 characters
	int   lineLength=0;
    char  keyW1 [256];//key word 1
    char  keyW2 [256];//key word 2
    char  keyW3 [256];//key word 3

    FILE * FilePointer = fopen(configfilename, "rb");
	if(FilePointer!=NULL) {	
		lineLength=ReadHeaderLine(FilePointer, line);//read and check the first nonempty line which starts with '#'
		if (lineLength==-1) {// if there are no one line which starts with '#'
			printf("%s has a wrong format! \n", configfilename);
		    printf("new magnoom.cfg will be created.\n");
		}else{
		    sscanf(line, "# %s %s %s", keyW1, keyW2, keyW3 );
		    if(strncmp(keyW1, "begin",5)!=0 || strncmp(keyW2, "magnoom",7)!=0 || strncmp(keyW3, "config",  6)!=0){
		        //if the first line isn't "# magnoom config file"
		    	printf("%s has wrong header of wrong file format! \n", configfilename);
		    	lineLength=-1;
		    }
		}
		//READING HEADER
		if (lineLength!=-1){
			do{
				lineLength = ReadHeaderLine(FilePointer, line);
				sscanf(line, "# %s %s %s", keyW1, keyW2, keyW3 );
				//printf("%s %s %s\n", keyW1, keyW2, keyW3);
				if (strncmp(keyW1, "ax:",3)==0) {
					sscanf(keyW2, "%f", &abc[0][0] );
					printf("ax=%f\n", abc[0][0]);					
				}else if (strncmp(keyW1, "ay:",3)==0) {
					sscanf(keyW2, "%f", &abc[0][1] );
					printf("ay=%f\n", abc[0][1]);
				}else if (strncmp(keyW1, "az:",3)==0) {
					sscanf(keyW2, "%f", &abc[0][2] );
					printf("az=%f\n", abc[0][2]);					
				}else if (strncmp(keyW1, "bx:",3)==0) {
					sscanf(keyW2, "%f", &abc[1][0] );
					printf("bx=%f\n", abc[1][0]);					
				}else if (strncmp(keyW1, "by:",3)==0) {
					sscanf(keyW2, "%f", &abc[1][1] );
					printf("by=%f\n", abc[1][1]);
				}else if (strncmp(keyW1, "bz:",3)==0) {
					sscanf(keyW2, "%f", &abc[1][2] );
					printf("bz=%f\n", abc[1][2]);					
				}else if (strncmp(keyW1, "cx:",3)==0) {
					sscanf(keyW2, "%f", &abc[2][0] );
					printf("cx=%f\n", abc[2][0]);					
				}else if (strncmp(keyW1, "cy:",3)==0) {
					sscanf(keyW2, "%f", &abc[2][1] );
					printf("cy=%f\n", abc[2][1]);
				}else if (strncmp(keyW1, "cz:",3)==0) {
					sscanf(keyW2, "%f", &abc[2][2] );
					printf("cz=%f\n", abc[2][2]);					
				}else if (strncmp(keyW1, "Na:",3)==0) {
					sscanf(keyW2, "%d", &uABC[0] );
					printf("Na=%d\n", uABC[0]);					
				}else if (strncmp(keyW1, "Nb:",3)==0) {
					sscanf(keyW2, "%d", &uABC[1] );
					printf("Nb=%d\n", uABC[1]);
				}else if (strncmp(keyW1, "Nc:",3)==0) {
					sscanf(keyW2, "%d", &uABC[2] );
					printf("Nc=%d\n", uABC[2]);					
				}else if (strncmp(keyW1, "Shells:",7)==0) {
					sscanf(keyW2, "%d", &ShellNumber );
					printf("Shells=%d\n", ShellNumber);					
				}else if (strncmp(keyW1, "BCa:",3)==0) {
					sscanf(keyW2, "%d", &Boundary[0] );
					printf("BCa=%d\n", Boundary[0]);					
				}else if (strncmp(keyW1, "BCb:",3)==0) {
					sscanf(keyW2, "%d", &Boundary[1] );
					printf("BCc=%d\n", Boundary[1]);					
				}else if (strncmp(keyW1, "BCc:",3)==0) {
					sscanf(keyW2, "%d", &Boundary[2] );
					printf("BCc=%d\n", Boundary[2]);					
				}
			}while(strncmp(keyW1, "end",3)!=0 || strncmp(keyW2, "magnoom",7)!=0 || strncmp(keyW3, "config",  6)!=0);
		}     
			AtomsPerBlock = sizeof(Block)/sizeof(float)/3;
			free(RadiusOfShell);
			RadiusOfShell = (float *)calloc(ShellNumber , sizeof(float));  
			free(NeighborsPerAtom);
			NeighborsPerAtom = (int *)calloc(AtomsPerBlock, sizeof(int));
			// total number of neighbour pairs per whole map of neighbours
 			NOS=AtomsPerBlock*uABC[0]*uABC[1]*uABC[2]; // number of spins
			NOS_AL=AtomsPerBlock*uABC[1]*uABC[2]; // number of spins per A layer
			NOS_BL=AtomsPerBlock*uABC[0]*uABC[2]; // number of spins per B layer
			NOS_CL=AtomsPerBlock*uABC[0]*uABC[1]; // number of spins per C layer

			iNOS = 1.0/NOS;

			NOB=uABC[0]*uABC[1]*uABC[2]; // number of Blocks
			NOB_AL=uABC[1]*uABC[2]; // number of spins per A layer
			NOB_BL=uABC[0]*uABC[2]; // number of spins per B layer
			NOB_CL=uABC[0]*uABC[1]; // number of spins per C layer
		// when everything is done
		printf("Done!\n");
		fclose(FilePointer);
	}else{
		printf("Cannot open file: %s \n", configfilename);
		printf("new magnoom.cfg will be created.\n");
	}
}


void setupTweakBar()
{
/*  Global settings for the bar-menu  */
	TwDefine(" GLOBAL iconpos=topleft "); // icons go to top-left corner of the window
	TwDefine(" GLOBAL iconalign=horizontal "); // icons will be aligned horizontally
	TwDefine(" GLOBAL contained=true "); // bars cannot move outside of the window
	TwDefine(" GLOBAL iconmargin='1 8'  "); // icons will be displayed at 1 and 16 pixels from the horizontal and vertical window borders respectively
/*  Help Bar F1 */
	help_bar = TwGetBarByIndex(0);
	TwDefine(" TW_HELP size='440 530' color='70 100 100'");
	TwDefine(" TW_HELP help='F1: show/hide (this) Help bar' "); // change default tweak bar size and color

/*  View Bar F2 */
    view_bar = TwNewBar("View");
    TwDefine(" View iconified=true "); 
    TwDefine(" View size='220 530' color='100 100 70' alpha=200 "); // change default tweak bar size and color
    TwDefine(" View help='F2: show/hide View bar' "); // change default tweak bar size and color
	{
	TwEnumVal		enColorsTw[] = { {WHITE,"White"}, {BLACK, "Black"}, {RED, "Red"}, {GREEN, "Green"}, {BLUE, "Blue"} };
	TwType			TV_TYPE_COLOR = TwDefineEnum("BG_Color", enColorsTw, 5);
	TwAddVarRW(view_bar, "Background color", TV_TYPE_COLOR, &WhichColor, "help='Background color for 3D scene' ");
	}

	{
	TwEnumVal		enProjectionsTw[] = { {ORTHO, "Orthogonal"}, {PERSP, "Perspective"} };
	TwType			TV_TYPE_PROJ = TwDefineEnum("ProjectionType", enProjectionsTw, 2);
	TwAddVarRW(view_bar, "Projection", TV_TYPE_PROJ, &WhichProjection, "keyIncr='p' help='Type of 3D projection' ");
	}

	{
	TwEnumVal		enVectorModeTw[] = {{ARROW1, "Arrows"} 
										,{CONE1,   "Cones"} 
										,{CANE,    "Canes"} 
										,{uPOINT,  "Points"}
										,{BOX1,    "Boxes"}
									};
	TwType			TV_TYPE_VEC_MOD = TwDefineEnum("Type_of_vectors", enVectorModeTw, 5);
	TwAddVarCB(view_bar, "Type of vectors", TV_TYPE_VEC_MOD, CB_SetVectorMode, CB_GetVectorMode, &WhichVectorMode, "keyIncr='v' keyDecr='V' help='Type of 3D vectors' ");
	}

	TwAddVarCB(view_bar, "Pivot", TW_TYPE_FLOAT, CB_SetPivot, CB_GetPivot,  &Pivot, " min=0 max=1 step=0.01 help='Pivot of 3D arrow.' ");
	TwAddVarCB(view_bar, "Faces", TW_TYPE_INT32, CB_SetFaces, CB_GetFaces,  &arrowFaces, " min=3 max=20 step=1 help='Number of faces for 3D arrow.' ");
	TwAddVarCB(view_bar, "Scale", TW_TYPE_FLOAT, CB_SetScale, CB_GetScale,  &Scale, " min=0.1 max=2.5 step=0.01 keyIncr='+' keyDecr='-' help='Scale the vectors.' ");

	TwAddVarRW(view_bar, "Show basis", TW_TYPE_BOOL32, &AxesOn, " key=CTRL+o ");
	TwAddVarRW(view_bar, "Show box", TW_TYPE_BOOL32, &BoxOn, " key=CTRL+b ");

	TwAddSeparator(view_bar, "view_sep1", NULL);
	{
	TwEnumVal		enSliceModeTw[] = { 
										{A_AXIS, "a-axis"}, 
										{B_AXIS, "b-axis"}, 
										{C_AXIS, "c-axis"},
										{FILTER, "filter"}
									  };
	TwType			TV_TYPE_VEC_MOD = TwDefineEnum("Slicing", enSliceModeTw, 4);
	TwAddVarCB(view_bar, "Slicing mode", TV_TYPE_VEC_MOD, CB_SetSliceMode, CB_GetSliceMode, &WhichSliceMode, "keyIncr='/' keyDecr='?' help='Slising plane perpenticulat to the choosen axis' ");
	}
	TwAddVarCB(view_bar, "T_max", TW_TYPE_FLOAT, CB_SetThetaMax, CB_GetThetaMax, &theta_max1, " label='Theta max' min=0 max=3.141592 step=0.01 help='max value for polar angle theta'");
	TwAddVarCB(view_bar, "T_min", TW_TYPE_FLOAT, CB_SetThetaMin, CB_GetThetaMin, &theta_min1, " label='Theta min' min=0 max=3.141592 step=0.01 help='min value for polar angle theta'");

	TwAddVarCB(view_bar, "F_max", TW_TYPE_INT32, CB_SetPhiMax, CB_GetPhiMax, &phi_max1, " label='Phi max' min=0 max=360 step=1 help='max value for azimuthal angle phi'");
	TwAddVarCB(view_bar, "F_min", TW_TYPE_INT32, CB_SetPhiMin, CB_GetPhiMin, &phi_min1, " label='Phi min' min=0 max=360 step=1 help='min value for azimuthal angle phi'");

	TwAddSeparator(view_bar, "view_sep2", NULL);

	TwAddVarRW(view_bar, "Intensity", TW_TYPE_FLOAT, &g_LightMultiplier, 
	" label='Light intensity' min=0.1 max=4 step=0.02 help='Increase/decrease the light power.' group='Light' ");
	TwAddVarRW(view_bar, "LightDir", TW_TYPE_DIR3F, &g_LightDirection, 
	" label='Light direction' opened=false help='Change the light direction.' group='Light'");
	temp_color[0] = 230;
	temp_color[1] = 230;
	temp_color[2] = 255;
	TwSetParam(view_bar, "LightDir", "arrowcolor", TW_PARAM_INT32, 3, temp_color);

	TwAddVarRW(view_bar, "Light_On_Off", TW_TYPE_BOOL32, &Light_On, 
	" label='Light On/Off' key=l help='Reflectins' group='Light'");

	TwAddVarRW(view_bar, "CamAng", TW_TYPE_FLOAT, &PerspSet[0], 
	" label='camera angle' min=1 max=120 help='camera angle' group='Camera control'");
	TwAddVarRW(view_bar, "PosX", TW_TYPE_FLOAT, &TransXYZ[0], 
	" label='position in x' min=-1000 max=1000 help='camera position along X-axis' group='Camera control'");
	TwAddVarRW(view_bar, "PosY", TW_TYPE_FLOAT, &TransXYZ[1], 
	" label='position in y' min=-1000 max=1000 help='camera position along Y-axis' group='Camera control'");
	TwAddVarRW(view_bar, "PosZ", TW_TYPE_FLOAT, &TransXYZ[2], 
	" label='position in z' min=-1000 max=1000 help='camera position along Z-axis' group='Camera control'");

	TwAddVarRW(view_bar, "RotX", TW_TYPE_FLOAT, &Rot[0], 
	" label='turn around x' min=-360 max=360 help='rotate camera around X-axis' group='Camera control'");
	TwAddVarRW(view_bar, "RotY", TW_TYPE_FLOAT, &Rot[1], 
	" label='turn around y' min=-360 max=360 help='rotate camera around Y-axis' group='Camera control'");
	TwAddVarRW(view_bar, "RotZ", TW_TYPE_FLOAT, &Rot[2], 
	" label='turn around z' min=-360 max=360 help='rotate camera around Z-axis' group='Camera control'");

	TwAddVarRW(view_bar, "CamBank", TW_TYPE_INT32, &CurrentCameraPositionBank, 
	" label='Current camera' min=0 max=4 group='Camera read/write'");
	TwAddButton(view_bar, "Read Camera", CB_GetCameraPosition, NULL, "label='read camera pos.' ");
	TwAddButton(view_bar, "Write Camera", CB_SaveCameraPosition, NULL, "label='save camera pos.' ");


	TwAddVarCB(view_bar, "ColSh", TW_TYPE_INT32, CB_SetColorShift, CB_GetColorShift, &ColorShift, 
	" label='Rotate hue' min=0 max=360 help='rotate color hue in xy-plane' group='HSV color map'");
	TwAddVarCB(view_bar, "InvHue", TW_TYPE_BOOL32, CB_SetInvHue, CB_GetInvHue, &InvertHue, 
	" label='Invert hue' help='invert RGB to RBG color hue' group='HSV color map'");
	TwAddVarCB(view_bar, "InvVal", TW_TYPE_BOOL32, CB_SetInvVal, CB_GetInvVal, &InvertValue, 
	" label='Invert value' help='invert black to white' group='HSV color map'");
	{
	TwEnumVal		enWhichColorSchemeTw[] = { {RGB, "RGB"}, {RYGB, "RYGB"}};
	TwType			TV_TYPE_COL_SCH = TwDefineEnum("ColorSh", enWhichColorSchemeTw, 2);
	TwAddVarCB(view_bar, "Color scheme", TV_TYPE_COL_SCH, CB_SetColorScheme, CB_GetColorScheme, &WhichColorScheme, "help='Type of 3D vectors' group='HSV color map' ");
	}

/*  Hamiltonian parameters&controls F3 */
	control_bar = TwNewBar("Parameters&Controls");
    TwDefine(" Parameters&Controls iconified=true "); 
    TwDefine(" Parameters&Controls size='220 530' color='100 70 100' alpha=200 "); // change default tweak bar size and color
    TwDefine(" Parameters&Controls help='F3: show/hide Control bar' "); // change default tweak bar size and color

	TwAddButton(control_bar, "Run", CB_Run, NULL, "key='r' label='RUN simulation' ");
	TwAddVarRW(control_bar, "Record", TW_TYPE_BOOL32, &Record, 
	"label='Recording' true='Rec.' false='Pause' help='Recording <sx>, <sy>, <sz> in each iteration'");
	TwAddVarRW(control_bar, "Rec_Iteration", TW_TYPE_INT32, &rec_iteration, 
	"label='Every iteration' min=1 max=1000 step=1 ");
	TwAddButton(control_bar, "Clean the record", CB_CleanSxSySzFile, NULL, "label= 'Clean the record' help='clean the output file with <sx>, <sy>, <sz>' ");
	TwAddButton(control_bar, "Reset iterations", CB_ResetIterations, NULL, "label='Reset iterations' ");

	TwAddSeparator(control_bar, "sep-3", NULL);

	TwAddVarRW(control_bar, "BCinA", TW_TYPE_BOOL32, &Boundary[0], 
	"label='along a' group='Boundary conditions' true='periodic' false='open' help='set boundary conditions along translation vector 'a' '");
	TwAddVarRW(control_bar, "BCinB", TW_TYPE_BOOL32, &Boundary[1], 
	"label='along b' group='Boundary conditions' true='periodic' false='open' help='set boundary conditions along translation vector 'b' '");
	TwAddVarRW(control_bar, "BCinC", TW_TYPE_BOOL32, &Boundary[2], 
	"label='along c' group='Boundary conditions' true='periodic' false='open' help='set boundary conditions along translation vector 'c' '");
    
    TwAddVarRW(control_bar, "Preces", TW_TYPE_BOOL32, &Precession, 
	"label='precession' group='LLG' true='On' false='Off' help='On/Off precession'");
	TwAddVarRW(control_bar, "Damping", TW_TYPE_FLOAT, &damping, 
	"label='Damping' min=0 max=100 step=0.01 group='LLG' ");
	TwAddVarRW(control_bar, "Time_step", TW_TYPE_FLOAT, &t_step, 
	"label='Time step' min=0 max=0.5 step=0.001   group='LLG' ");
	TwAddVarRW(control_bar, "temperature", TW_TYPE_FLOAT, &Temperature, 
	"label='k_b*T' min=0 max=100 step=0.01 group='LLG' ");

    TwAddSeparator(control_bar, "sep-2", NULL);
	TwAddVarCB(control_bar, "FieldTheta", TW_TYPE_FLOAT, CB_SetHfieldTheta, CB_GetHfieldTheta, &VHtheta, "label='H theta'  step=1 help='Change the direction of applied field' ");
	TwAddVarCB(control_bar, "FieldPhi", TW_TYPE_FLOAT, CB_SetHfieldPhi, CB_GetHfieldPhi, &VHphi, "label='H phi' step=1  help='Change the direction of applied field' ");
	// TwAddVarCB(control_bar, "FieldDir", TW_TYPE_DIR3F, CB_SetHfieldDir, CB_GetHfieldDir, VHf, 
	// "label='Field direction' opened=true help='Change the direction of applied field' ");
	// temp_color[0] = 55;
	// temp_color[1] = 55;
	// temp_color[2] = 155;
	// TwSetParam(control_bar, "FieldDir", "arrowcolor", TW_PARAM_INT32, 3, temp_color);
	// // TwAddVarRW(control_bar, "Field", TW_TYPE_FLOAT, &Hf, 
	// // "label='Field' help='The value of uniaxial anisotropy' ");
	TwAddVarCB(control_bar, "Field", TW_TYPE_FLOAT, CB_SetHfield, CB_GetHfield, &Hf, 
	"label='Field'  min=0 step=0.00001 help='The value of uniaxial anisotropy' ");

	TwAddSeparator(control_bar, "control_sep1", NULL);
	TwAddVarCB(control_bar, "FieldDir", TW_TYPE_DIR3F, CB_SetHfieldXYZ, CB_GetHfieldXYZ, &VHf, 
	" label='Field direction' opened=false help='Change the applied field direction.' ");

	TwAddSeparator(control_bar, "control_sep2", NULL);

    TwAddSeparator(control_bar, "sep-1", NULL);


	TwAddVarRW(control_bar, "KudDir", TW_TYPE_DIR3F, &VKu, 
	"label='Ku' opened=true help='The direction of uniaxial anisotropy' ");
	temp_color[0] = 55;
	temp_color[1] = 155;
	temp_color[2] = 55;
	TwSetParam(control_bar, "KudDir", "arrowcolor", TW_PARAM_INT32, 3, temp_color);
	TwAddVarRW(control_bar, "Ku", TW_TYPE_FLOAT, &Ku, 
	"label='Ku' help='The value of uniaxial anisotropy' ");

	//////////////////////////////////////////
	TwAddSeparator(control_bar, "sep0", NULL);
    //////////////////////////////////////////
	TwAddVarRW(control_bar, "Kcub", TW_TYPE_FLOAT, &Kc, 
	"label='Kc' help='The value of cubic anisotropy' ");
    //////////////////////////////////////////
    TwAddSeparator(control_bar, "sep1", NULL);
    //////////////////////////////////////////
	for(int s=0; s<ShellNumber; s++)
	{
	snprintf(shortBufer,50,"J%1i",s);
	TwAddVarRW(control_bar, shortBufer, TW_TYPE_FLOAT, &Jij[s], 
	"help='Heisenberg exchange' ");
	TwSetParam(control_bar, shortBufer, "label",  TW_PARAM_CSTRING , 1, shortBufer);
	}
	//////////////////////////////////////////
	TwAddSeparator(control_bar, "sep2", NULL);
	//////////////////////////////////////////	
	for(int s=0; s<ShellNumber; s++)
	{
	snprintf(shortBufer,50,"B%1i",s);
	TwAddVarRW(control_bar, shortBufer, TW_TYPE_FLOAT, &Bij[s], 
	"help='Bi-quadratic exchange' ");
	TwSetParam(control_bar, shortBufer, "label",  TW_PARAM_CSTRING , 1, shortBufer);
	}
	//////////////////////////////////////////
	TwAddSeparator(control_bar, "sep3", NULL);
	//////////////////////////////////////////
	for(int s=0; s<ShellNumber; s++)
	{
	snprintf(shortBufer,50,"D%1i",s);
	TwAddVarRW(control_bar, shortBufer, TW_TYPE_FLOAT, &Dij[s], 
	"help='Dzyaloshinskii-Moriya' ");
	TwSetParam(control_bar, shortBufer, "label",  TW_PARAM_CSTRING , 1, shortBufer);
	}
	//////////////////////////////////////////
	TwAddSeparator(control_bar, "sep4", NULL);
	//////////////////////////////////////////
	TwAddVarRW(control_bar, "ST param.", TW_TYPE_FLOAT, &Cu, 
	"help='Dzyaloshinskii-Moriya' ");
	TwAddVarRW(control_bar, "CurrentDir", TW_TYPE_DIR3F, &VCu, 
	"label='cur. dir.' opened=true help='The polarization direction of electric current' ");

/*  Initial state F4 */
	initial_bar = TwNewBar("Initial_State");
	TwDefine(" Initial_State iconified=true "); 
	TwDefine(" Initial_State size='220 530' color='70 70 100'  alpha=200"); // change default tweak bar size and color
	TwDefine(" Initial_State help='F4: show/hide Initial state bar' "); // change default tweak bar size and color
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
										{SKYRMION_L,"Sk. lattice"	        }
									};
	TwType			TV_TYPE_INI_STATE = TwDefineEnum("IniState", enIniStateTw, 13);
	TwAddVarRW(initial_bar, "Choose ini. state", TV_TYPE_INI_STATE, &WhichInitialState, "help='Choose initial spin configuration'");
	}

	TwAddVarRW(initial_bar, "chSize", TW_TYPE_FLOAT,  &chSize, 
	" min=0 max=100000 step=0.5 help='characteristic size of modulated state: Skyrmion diameter or spiral period' ");

	TwAddVarRW(initial_bar, "chDir", TW_TYPE_DIR3F,  &chDir, 
	" help='direction of modulations e.g. k-vector of the spin spiral.' ");

	TwAddButton(initial_bar, "Set initial", CB_SetInitial, NULL, "key=I label='insert initial state' ");

	TwAddSeparator(initial_bar, "sep01", NULL);
	TwAddVarRW(initial_bar, "Degrees", TW_TYPE_FLOAT,  &RotateAllSpins, 
	" min=-360 max=360 step=1 help='Rotate all spins about characteristic direction' ");
	TwAddButton(initial_bar, "Rotate spins", CB_RotateAllSpins, NULL, "label='rotate all spins' ");

	TwAddSeparator(initial_bar, "sep1", NULL);
	TwAddButton(initial_bar, "Invert X", CB_InvertX, NULL, "label='invert n_x component' ");
	TwAddButton(initial_bar, "Invert Y", CB_InvertY, NULL, "label='invert n_y component' ");
	TwAddButton(initial_bar, "Invert Z", CB_InvertZ, NULL, "label='invert n_z component' ");

	TwAddSeparator(initial_bar, "sep02", NULL);
	TwAddVarRW(initial_bar, "Input file name:", TW_TYPE_CSSTRING(sizeof(inputfilename)), inputfilename, "");
	TwAddButton(initial_bar, "Read from CSV", CB_ReadCSV, NULL, "label='read from *.csv file' ");
	TwAddButton(initial_bar, "Read from OVF", CB_ReadOVF, NULL, "label='read from *.ovf file' ");

	TwAddSeparator(initial_bar, "sep2", NULL);
	TwAddVarRW(initial_bar, "Save slice", TW_TYPE_BOOL32, &save_slice, " label='Save slice' help='Save current slice only' ");
	{
	TwEnumVal		enAverage_mode[] = {{ALONG_A, 	"Along a-axis"	}, 
										{ALONG_B, 	"Along b-axis"	}, 
										{ALONG_C, 	"Along c-axis"	}, 
										{ALONG_0, 	"No avearge  "	}
									};

	TwType			TV_TYPE_AVERAGE_MODE = TwDefineEnum("Average mode", enAverage_mode, 4);
	TwAddVarRW(initial_bar, "Choose average mode", TV_TYPE_AVERAGE_MODE, &WhichAverageMode, "help='Choose type of average mode'");}

	TwAddVarRW(initial_bar, "Output file name:", TW_TYPE_CSSTRING(sizeof(outputfilename)), outputfilename, ""); 
	TwAddButton(initial_bar, "Write to CSV", CB_SaveCSV, NULL, "label='write to *.csv file' ");		

/*  AC field F5 */
	ac_field_bar = TwNewBar("AC_Field");
	TwDefine(" AC_Field iconified=true "); 
	TwDefine(" AC_Field size='220 530' color='100 70 70'  alpha=200"); // change default tweak bar size and color
	TwDefine(" AC_Field help='F5: show/hide AC Field bar' "); // change default tweak bar size and color

	{
	TwEnumVal		enenACFieldTw[] = { {SIN_FIELD,  	"AC Sin(omega*t) "		},
										{GAUSSIAN_FIELD,"Gaussian field pulse" }};
	TwType			TV_TYPE_INI_STATE = TwDefineEnum("AC-field", enenACFieldTw, 2);
	TwAddVarRW(ac_field_bar, "Type of AC field", TV_TYPE_INI_STATE, &WhichACField, "help='Choose type of signal for time dependent magnetic field'");
	}

	TwAddVarRW(ac_field_bar, "ACfieldDir", TW_TYPE_DIR3F, &VHac, 
	"label='Field direction' opened=true help='Change the direction of applied field' ");
	temp_color[0] = 55;
	temp_color[1] = 55;
	temp_color[2] = 155;
	TwSetParam(ac_field_bar, "ACfieldDir", "arrowcolor", TW_PARAM_INT32, 3, temp_color);
	TwAddVarRW(ac_field_bar, "AC field", TW_TYPE_FLOAT, &Hac, 
	"label='AC field' help='The value of AC field amplitude' ");

	TwAddVarCB(ac_field_bar, "Period/Width", TW_TYPE_FLOAT, CB_SetACPeriod, CB_GetACPeriod,  &Period_dc, "min=0 help='period of sin-field or width of gaussian pulse field' ");
	TwAddVarCB(ac_field_bar, "Omega=2*pi*P", TW_TYPE_FLOAT, CB_SetOmega, CB_GetOmega,  &Omega_dc, "min=0 help='period of sin-field or width of gaussian pulse field' ");
	TwAddVarRW(ac_field_bar, "AC field on/off", TW_TYPE_BOOL32, &AC_FIELD_ON, 
	"keyIncr='f' label='AC/DC on/off' true='on' false='off' help='On/off ac field'");

/*  Info bar F11 */
	info_bar = TwNewBar("Info");
	TwDefine(" Info refresh=0.5 ");
	TwDefine(" Info iconified = false movable = false alwaysbottom=true resizable=false fontstyle=fixed fontsize=2"); 
	TwDefine(" Info help='F11: show/hide info-bar' "); // change default tweak bar size and color
	TwDefine(" Info color='10 10 10' alpha=0 "); // change default tweak bar size and color
	TwDefine(" Info help='F11: show/hide info-bar' "); // change default tweak bar size and color
	TwDefine(" Info position = '1 30' size ='200 500' valueswidth=120"); // change default tweak bar size and color
	TwAddVarRO(info_bar, "R/S", TW_TYPE_BOOL32,  &Play, "true='RUNING' false='STOPED' ");
	TwAddVarRO(info_bar, "Rec.", TW_TYPE_BOOL32,  &Record, "true='On' false='Off' ");
	TwAddVarRO(info_bar, "ACF", TW_TYPE_BOOL32,  &AC_FIELD_ON, "true='On' false='Off' help='AC filed on/off'");
	TwAddSeparator(info_bar, "sep-0", NULL);
	TwAddVarRO(info_bar, "NPB", TW_TYPE_INT32,  &AtomsPerBlock, "help='number of atoms per block' ");
	TwAddVarRO(info_bar, "N_a", TW_TYPE_INT32,  &uABC[0], "help='translations along a' ");
	TwAddVarRO(info_bar, "N_b", TW_TYPE_INT32,  &uABC[1], "help='translations along b' ");
	TwAddVarRO(info_bar, "N_c", TW_TYPE_INT32,  &uABC[2], "help='translations along c' ");	
	TwAddVarRO(info_bar, "NOS", TW_TYPE_INT32,  &NOS, "help='Number of spins' ");
	TwAddSeparator(info_bar, "sep", NULL);
	TwAddVarRO(info_bar, "FPS", TW_TYPE_FLOAT,  &FPS, "help='Frame per second' ");
	TwAddVarRO(info_bar, "ITR", TW_TYPE_INT32,  &currentIteration, "help='Total number of iterations' ");
	TwAddVarRO(info_bar, "IPS", TW_TYPE_FLOAT,  &IPS, "help='Iterations per secon' ");

	TwAddSeparator(info_bar, "sep0", NULL);
	TwAddVarRO(info_bar, "Etot", TW_TYPE_DOUBLE, &totalEnergy, " help='Total energy' ");
	TwAddVarRO(info_bar, "e", TW_TYPE_DOUBLE, &perSpEnergy, " help='Energy density per spin'");

	TwAddSeparator(info_bar, "sep01", NULL);
	TwAddVarRO(info_bar, "e0", TW_TYPE_DOUBLE, &totalEnergyFerro, " help='Energy density for ferromagnetic state' ");
	TwAddVarRO(info_bar, "e-e0", TW_TYPE_DOUBLE, &perSpEnergyMinusFerro, " help='Energy density per spin wrt ferromagnetic state'");

	TwAddSeparator(info_bar, "sep1", NULL);	
	TwAddVarRO(info_bar, "M_x", TW_TYPE_DOUBLE, &Mtot[0], " help='x-component of total moment' ");
	TwAddVarRO(info_bar, "M_y", TW_TYPE_DOUBLE, &Mtot[1], " help='y-component of total moment' ");
	TwAddVarRO(info_bar, "M_z", TW_TYPE_DOUBLE, &Mtot[2], " help='z-component of total moment' ");

	TwAddSeparator(info_bar, "sep2", NULL);
	TwAddVarRO(info_bar, "m_x", TW_TYPE_DOUBLE, &mtot[0], " help='x-component of average moment per spin' ");
	TwAddVarRO(info_bar, "m_y", TW_TYPE_DOUBLE, &mtot[1], " help='y-component of average moment per spin' ");
	TwAddVarRO(info_bar, "m_z", TW_TYPE_DOUBLE, &mtot[2], " help='z-component of average moment per spin' ");
	TwAddSeparator(info_bar, "sep3", NULL);
	TwAddVarRO(info_bar, "max_torque", TW_TYPE_DOUBLE, &MAX_TORQUE, " help='maximum torque acting on the spin' ");
}


void MouseButton( int button, int state, int x, int y ) // called when the mouse button transitions down or up
{
	int b = 0;			// LEFT, MIDDLE, or RIGHT
if( !TwEventMouseButtonGLUT(button, state, x, y) )  // send event to AntTweakBar
    { // event has not been handled by AntTweakBar
      // your code here to handle the event
		switch( button )
		{
			case GLUT_LEFT_BUTTON:		// 0
				b = LEFT;		break;

			case GLUT_MIDDLE_BUTTON:	// 1
				b = MIDDLE;		break;

			case GLUT_RIGHT_BUTTON:		// 2
				b = RIGHT;		break;

			case 3:		                // scrol up
			TransXYZ[2]+=0.5;
			//glui_window->sync_live( );	
			break;

			case 4:		                // scrol down
			TransXYZ[2]-=0.5;
			//glui_window->sync_live( );	
			break;

			default:
				b = 0;
				fprintf( stderr, "Unknown mouse button: %d\n", button );
		}

		// button down sets the bit, up clears the bit:
		if( state == GLUT_DOWN )
		{
			Xmouse = x;
			Ymouse = y;
			ActiveButton |= b;		// set the proper bit
		}
		else
		{
			ActiveButton &= ~b;		// clear the proper bit
		}
	}
}


void MouseMotion( int x, int y ) // called when the mouse moves while a button is down
{
	int dx = x - Xmouse;		// change in mouse coords
	int dy = y - Ymouse;
if( !TwEventMouseMotionGLUT(x, y) )  // send event to AntTweakBar
    { // event has not been handled by AntTweakBar
      // your code here to handle the event
      if( ( ActiveButton & LEFT ) != 0 )
		{
			Rot[2] += ( dx )*0.1f;
			Rot[0] += ( dy )*0.1f;
		}

		if( ( ActiveButton & MIDDLE ) != 0 )
		{
			TransXYZ[0]+=dx*0.05;
			TransXYZ[1]-=dy*0.05;	
		}

		if( ( ( ActiveButton & RIGHT ) != 0 ) & (true))//WhichProjection == PERSP
		{
			TransXYZ[2]+=(dx-dy)*0.05;
		}

		Xmouse = x;			// new current position
		Ymouse = y;	
		// synchronize the GLUI display with the variables:
		//glui_window->sync_live( );		
		// glutSetWindow( GLUT_window );
		// glutPostRedisplay( );
    }
}

// the keyboard callback:
void Keyboard( unsigned char c, int x, int y )
{
if( !TwEventKeyboardGLUT(c, x, y) )  // send event to AntTweakBar 
  { // event has not been handled by AntTweakBar
    // your code here to handle the event	
	// if( false ) fprintf( stderr, "Keyboard: '%c' (0x%0x)\n", c, c );

		switch( c )
		{
			case '<':
				switch(WhichSliceMode)
                {case A_AXIS:
					if (A_layer_min>1) {A_layer_min-=1; ChangeVectorMode(0);}
				 break;
				 case B_AXIS:
				    if (B_layer_min>1) {B_layer_min-=1; ChangeVectorMode(0);}
				 break;
				 case C_AXIS:
				 	if (C_layer_min>1) {C_layer_min-=1; ChangeVectorMode(0);}
				 break;	
				 default:
				 break;		 
				}
			break;

			case '>':
				switch(WhichSliceMode)
                {case A_AXIS:
					if (A_layer_min<A_layer_max) {A_layer_min+=1; ChangeVectorMode(0);}
				 break;
				 case B_AXIS:
				    if (B_layer_min<B_layer_max) {B_layer_min+=1; ChangeVectorMode(0);}
				 break;
				 case C_AXIS:
				 	if (C_layer_min<C_layer_max) {C_layer_min+=1; ChangeVectorMode(0);}
				 break;	
				 default:
				 break;		 
				}
			break;

			case ',':
				switch(WhichSliceMode)
                {case A_AXIS:
					if (A_layer_max>A_layer_min) {A_layer_max-=1; ChangeVectorMode (0);}
				 break;
				 case B_AXIS:
				    if (B_layer_max>B_layer_min) {B_layer_max-=1; ChangeVectorMode (0);}
				 break;
				 case C_AXIS:
				 	if (C_layer_max>C_layer_min) {C_layer_max-=1; ChangeVectorMode (0);}
				 break;
				 default:
				 break;			 
				}
			break;

			case '.':
				switch(WhichSliceMode)
                {case A_AXIS:
					if (A_layer_max<uABC[0]) {A_layer_max+=1; ChangeVectorMode (0);}
				 break;
				 case B_AXIS:
				    if (B_layer_max<uABC[1]) {B_layer_max+=1; ChangeVectorMode (0);}
				 break;
				 case C_AXIS:
				 	if (C_layer_max<uABC[2]) {C_layer_max+=1; ChangeVectorMode (0);}
				 break;	
				 default:
				 break;		 
				}
			break;

			case 'q':
			case 'Q':
				Rot[2] -= 0.25;
			break;

			case 'e':
			case 'E':
				Rot[2] += 0.25;
			break;

			case 'w':
			case 'W':
				//Rot[0] += 0.25;
				TransXYZ[2]-=0.5;
				break;
			case 's':
			case 'S':
				//Rot[0] -= 0.25;
				TransXYZ[2]+=0.5;
				break;
			case 'a':
			case 'A':
				Rot[0] -= 0.25;
				//TransXYZ[0]-=1;	
				break;
			case 'd':
			case 'D':
				Rot[0] += 0.25;
				//TransXYZ[0]+=1;	
				break;

			case 'x':
			case 'X':
				Buttons( XUP );
			break;

			case 'y':
			case 'Y':
				Buttons( YUP );
			break;

			case 'z':
			case 'Z':
				Buttons( ZUP );
			break;

			case ESCAPE:
				Buttons( QUIT );
				break;

			default:
				fprintf( stderr, "Don't know what to do with keyboard hit: '%c' (0x%0x)\n", c, c );
		}

	}
}



void KeyboardAdd( int key, int x, int y ){
int isiconified;
if( !TwEventSpecialGLUT(key, x, y) )  // send event to AntTweakBar TwEventSpecialGLUT
    { // event has not been handled by AntTweakBar
      // your code here to handle the event	
		switch( key )
		{
			case GLUT_KEY_UP:
				switch (WhichSliceMode)
				{case A_AXIS:
					if (A_layer_max < uABC[0]) {A_layer_max+=1; A_layer_min+=1; ChangeVectorMode(1);}
				 break;
				 case B_AXIS:
				    if (B_layer_max < uABC[1]) {B_layer_max+=1; B_layer_min+=1; ChangeVectorMode(1);}
				 break;
				 case C_AXIS:
				 	if (C_layer_max < uABC[2]) {C_layer_max+=1; C_layer_min+=1; ChangeVectorMode(1);}
				 break;	
				 default:
				 break;		 
				}
				break;
			case GLUT_KEY_DOWN:
				switch (WhichSliceMode)
				{case A_AXIS:
					if (A_layer_min>1) {A_layer_min-=1; A_layer_max-=1; ChangeVectorMode(1);}
				 break;
				 case B_AXIS:
				    if (B_layer_min>1) {B_layer_min-=1; B_layer_max-=1; ChangeVectorMode(1);}
				 break;
				 case C_AXIS:
				 	if (C_layer_min>1) {C_layer_min-=1; C_layer_max-=1; ChangeVectorMode(1);} 
				 break;
				 default:
				 break;
				}
				break;
			case  GLUT_KEY_F1:
				TwGetParam(help_bar, NULL, "iconified", TW_PARAM_INT32, 1, &isiconified);
				if (isiconified){
			    	TwDefine(" TW_HELP iconified=false ");				
				}else{
					TwDefine(" TW_HELP iconified=true ");
				}
				break;
			case  GLUT_KEY_F2:
				TwGetParam(view_bar, NULL, "iconified", TW_PARAM_INT32, 1, &isiconified);
				if (isiconified){
					TwDefine(" View iconified=false ");				
				}else{
					TwDefine(" View iconified=true ");
				}
				break;
			case  GLUT_KEY_F3:
				TwGetParam(control_bar, NULL, "iconified", TW_PARAM_INT32, 1, &isiconified);
				if (isiconified){
					TwDefine(" Parameters&Controls iconified=false ");				
				}else{
					TwDefine(" Parameters&Controls iconified=true ");
				}
				break;
			case  GLUT_KEY_F4:
				TwGetParam(initial_bar, NULL, "iconified", TW_PARAM_INT32, 1, &isiconified);
				if (isiconified){
					TwDefine(" Initial_State iconified=false ");				
				}else{
					TwDefine(" Initial_State iconified=true ");
				}
				break;
			case  GLUT_KEY_F5:
				TwGetParam(ac_field_bar, NULL, "iconified", TW_PARAM_INT32, 1, &isiconified);
				if (isiconified){
					TwDefine(" AC_Field iconified=false ");				
				}else{
					TwDefine(" AC_Field iconified=true ");
				}
				break;
			case  GLUT_KEY_F11:
				TwGetParam(info_bar, NULL, "iconified", TW_PARAM_INT32, 1, &isiconified);
				if (isiconified){
					TwDefine(" Info iconified=false ");				
				}else{
					TwDefine(" Info iconified=true ");
				}
				break;
			default:
				fprintf( stderr, "Don't know what to do with special key: '%c' (0x%0x)\n", key, key );
		}
	}
}

void Buttons( int id )
{
	switch( id )
	{
		case XUP:
			Xup( );
			glutSetWindow( GLUT_window );
			glutPostRedisplay( );
			break;

		case YUP:
			Yup( );
			glutSetWindow( GLUT_window );
			glutPostRedisplay( );
			break;

		case ZUP:
			Zup( );
			glutSetWindow( GLUT_window );
			glutPostRedisplay( );
			break;

		case QUIT:
			glutSetWindow( GLUT_window );	//
			glFinish( );					// gracefully close out the graphics
			glutDestroyWindow( GLUT_window );// gracefully close the graphics window
			exit( 0 );						// gracefully exit the program
			break;

		default:
			fprintf( stderr, "Don't know what to do with Button ID %d\n", id );
	}
}

void ChangeInitialState ( int id )
{
	InitSpinComponents( Px, Py, Pz, Sx, Sy, Sz, id );
	for (int i=0;i<NOS;i++) { bSx[i]=Sx[i]; bSy[i]=Sy[i]; bSz[i]=Sz[i];}
	ChangeVectorMode ( 1 );
}

void ChangeVectorMode( int id )
{
	switch ( id )
	{
		case 0: // change of mode e.g. from point to arrow
		ReallocateArrayDrawing();
		// Fill array for prototype (arrow or cane) array 
		UpdatePrototypeVerNorInd(vertexProto, normalProto, indicesProto, arrowFaces, WhichVectorMode,0); 
		// Fill big array for indecies for all arrows, cans, cones or boxes 
		UpdateIndices(indicesProto , IdNumProto, indices, IdNum, VCNumProto);
		case 1:	// change only the layer(s) for drawing
		UpdateVerticesNormalsColors(vertexProto, normalProto, VCNumProto, vertices, normals, colors, VCNum, Px, Py, Pz, bSx, bSy, bSz, WhichVectorMode);
		UpdateVBO(&vboIdV, &vboIdN, &vboIdC, &iboIdI,vertices, normals, colors, indices);
		break;
	}
}

void ChangeColorMap ( int id )
{
	switch (id)
	{
		case 0:
		if(WhichColorScheme==RGB)
		{
		HueMap[0]=HueMapRGB[0]+ColorShift;
		HueMap[1]=HueMapRGB[1]+ColorShift;
		HueMap[2]=HueMapRGB[2]+ColorShift;
		HueMap[3]=HueMapRGB[3]+ColorShift;
		HueMap[4]=HueMapRGB[4]+ColorShift;
		HueMap[5]=HueMapRGB[5]+ColorShift;
		} 	
		else
		{
		HueMap[0]=HueMapRYGB[0]+ColorShift;
		HueMap[1]=HueMapRYGB[1]+ColorShift;
		HueMap[2]=HueMapRYGB[2]+ColorShift;
		HueMap[3]=HueMapRYGB[3]+ColorShift;
		HueMap[4]=HueMapRYGB[4]+ColorShift;
		HueMap[5]=HueMapRYGB[5]+ColorShift;
		}
		InitRGB(RHue, GHue, BHue, HueMap);
		case 1:
		ChangeVectorMode(1);
	}
}

void HSVtoRGB(float Vec[3], float RGB[3], int inV, int inH )
{
	int H = inH*359+(1-2*inH)*((int) atan2int (Vec[1]+0.01,Vec[0])); //it's fast atan function [int deg], see MATH.cpp

    float maxV, minV, dV;
    if ((1-2*inV)*Vec[2]<0)
    {
    	minV = 0;
    	maxV = 1 - fabs(Vec[2]);
    }else{
    	minV = fabs(Vec[2]);
    	maxV = 1;    	
    }
    dV = maxV-minV;
	RGB[0] = RHue[H]*dV+minV;
	RGB[1] = GHue[H]*dV+minV;
	RGB[2] = BHue[H]*dV+minV;
}

void InitRGB(float* R, float* G, float* B, int *hueMap)
{//It is a slow function  because of many "if" statements but we are not going to call it very often
	int Red		= hueMap[0];
	int Yellow	= hueMap[1];
	int Green	= hueMap[2];
	int Cyan	= hueMap[3];
	int Blue	= hueMap[4]; 
	int Magenta	= hueMap[5];
	int I; 
	for (int i=0;i<720;i++)
	{
		I=i;
	if( (I>=Red) && (I< Yellow) )
		{
			R[i%360] =1.f;
			G[i%360] =1.f/(Yellow - Red) * (I - Red);
			B[i%360] =0.f;
		}
	    else if ( (I>=Yellow) && (I < Green) )
		{
			R[i%360] =1.f/(Green - Yellow) * (Yellow-I)+1.f;
			G[i%360] =1.f;
			B[i%360] =0.f;
		}
	    else if ( (I>=Green) && (I< Cyan) )
		{
			R[i%360] =0.f;
			G[i%360] =1.f;
			B[i%360] =1.f/(Cyan - Green) * (I - Green);
		}
	    else if  ( (I>=Cyan) && (I < Blue) )
		{
			R[i%360] =0.f;
			G[i%360] =1.f/(Blue - Cyan) * (Cyan-I)+1.f;
			B[i%360] =1.f;
		}
	    else if ( (I>=Blue) && (I< Magenta) )
		{
			R[i%360] =1.f/(Magenta - Blue) * (I - Blue);
			G[i%360] =0.f;
			B[i%360] =1.f;
		}
	    else if ( (I>=Magenta) && (I< (Red+360)) )
		{
			R[i%360] =1.f;
			G[i%360] =0.f;
			B[i%360] =1.f/(Red+360 - Magenta) * (Magenta-I)+1.f;
		}
	} 
}


void ReallocateArrayDrawing()
{
	free(vertexProto); free(normalProto); free(indicesProto);
	free(vertices); free(normals); free(colors); free(indices);

	int NOS_L=0;
	int NOB_L=0;
	switch( WhichSliceMode)
	{
		case A_AXIS:
		NOS_L=NOS_AL * (1+A_layer_max - A_layer_min);
		if (A_layer_max - A_layer_min<2){
			NOB_L=NOB_AL * (1+A_layer_max - A_layer_min);
		}else{
			NOB_L=2*NOB_AL+2*(A_layer_max - A_layer_min-1)*(uABC[1]-1+uABC[2]-1);
		}
		break;
		case B_AXIS:
		NOS_L=NOS_BL * (1+B_layer_max - B_layer_min);
		if (B_layer_max - B_layer_min<2){
			NOB_L=NOB_BL * (1+B_layer_max - B_layer_min);
		}else{
			NOB_L=2*NOB_BL + 2*(B_layer_max - B_layer_min-1)*(uABC[0]-1+uABC[2]-1);
		}
		break;
		case C_AXIS:
		NOS_L=NOS_CL * (1+C_layer_max - C_layer_min);
		if (C_layer_max - C_layer_min<2){
			NOB_L=NOB_CL * (1+C_layer_max - C_layer_min);
		}else{
			NOB_L=2*NOB_CL + 2*(C_layer_max - C_layer_min-1)*(uABC[0]-1+uABC[1]-1);
		}
		break;
		case FILTER:
		NOS_L=NOS;
		NOB_L=NOS/AtomsPerBlock;
		break;
	}

	switch( WhichVectorMode )
	{
		case ARROW1:		
			// arrowFaces - number of arrow faces
			ElNumProto = 5*arrowFaces-4; // number of triangles per arrow
			IdNumProto = 3*ElNumProto; // number of indixes per arrow
			VCNumProto = 3*(2*(1+arrowFaces)-2+4*arrowFaces+3*arrowFaces);
			// Allocate memory for arrow prototype  
			vertexProto  = (float  *)malloc(VCNumProto * sizeof( float  ));
			normalProto  = (float  *)malloc(VCNumProto * sizeof( float  ));
			indicesProto = (GLuint *)malloc(IdNumProto * sizeof( GLuint ));	
			ElNum = NOS_L * ElNumProto;
			IdNum = NOS_L * IdNumProto;
			VCNum = NOS_L * VCNumProto;
			// Allocate memory for all ARROW1 (spins) 
			vertices	= (float  *)malloc(VCNum * sizeof( float  ));
			normals 	= (float  *)malloc(VCNum * sizeof( float  ));
			colors 		= (float  *)malloc(VCNum * sizeof( float  ));
			indices		= (GLuint *)malloc(IdNum * sizeof( GLuint ));	
			break;

		case CONE1:		
			// arrowFaces - number of cone faces
			ElNumProto = 2*arrowFaces-2; // number of triangles per cone
			IdNumProto = 3*ElNumProto; // number of indixes per cone
			VCNumProto = 3*((arrowFaces)+3*arrowFaces);// number of vertecies per cone			
			// Allocate memory for arrow prototype  
			vertexProto  = (float  *)malloc(VCNumProto * sizeof( float  ));
			normalProto  = (float  *)malloc(VCNumProto * sizeof( float  ));
			indicesProto = (GLuint *)malloc(IdNumProto * sizeof( GLuint ));	
			ElNum = NOS_L * ElNumProto;
			IdNum = NOS_L * IdNumProto;
			VCNum = NOS_L * VCNumProto;	
			// Allocate memory for all CONES (spins) 
			vertices	= (float  *)malloc(VCNum * sizeof( float  ));
			normals 	= (float  *)malloc(VCNum * sizeof( float  ));
			colors 		= (float  *)malloc(VCNum * sizeof( float  ));
			indices		= (GLuint *)malloc(IdNum * sizeof( GLuint ));	
			break;

		case BOX1:		
			ElNumProto = 6*2; // number of triangles per box
			IdNumProto = 3*ElNumProto; // number of indixes per box
			VCNumProto = 6*4*3; // number of vertecies per box			
			// Allocate memory for box prototype  
			vertexProto  = (float  *)malloc(VCNumProto* sizeof( float  ));
			normalProto  = (float  *)malloc(VCNumProto* sizeof( float  ));
			indicesProto = (GLuint *)malloc(IdNumProto* sizeof( GLuint ));	
			ElNum = NOB_L * ElNumProto;
			IdNum = NOB_L * IdNumProto;
			VCNum = NOB_L * VCNumProto;	
			// Allocate memory for all boxes (spins averaged over Block) 
			vertices	= (float  *)malloc(VCNum * sizeof( float  ));
			normals 	= (float  *)malloc(VCNum * sizeof( float  ));
			colors 		= (float  *)malloc(VCNum * sizeof( float  ));
			indices		= (GLuint *)malloc(IdNum * sizeof( GLuint ));	
			break;

		case uPOINT:		
			// arrowFaces - number of arrow faces
			IdNumProto = 1; // number of indixes per point
			VCNumProto = 3; // number component of vertices per point
			// Allocate memory for point prototype
			vertexProto  = (float  *)malloc(VCNumProto * sizeof( float  ));
			normalProto  = NULL;
			indicesProto = (GLuint *)malloc(IdNumProto * sizeof( GLuint ));	
			IdNum = NOS_L * IdNumProto;
			VCNum = NOS_L * VCNumProto; // total number of all components of all vertices 
			// Allocate memory for all points (spins) 
			vertices	= (float  *)malloc(VCNum * sizeof( float  ));
			normals 	= NULL;
			colors 		= (float  *)malloc(VCNum * sizeof( float  ));
			indices		= (GLuint *)malloc(IdNum * sizeof( GLuint ));	
			break;

		case CANE:	// is default vector mode

		default:
			IdNumProto = 2; // number of indixes per cane
			VCNumProto = 3*2; // number component of vertices per cane
			// Allocate memory for cane prototype
			vertexProto  = (float  *)malloc(VCNumProto * sizeof( float  ));
			normalProto  = NULL;
			indicesProto = (GLuint *)malloc(IdNumProto * sizeof( GLuint ));	
			IdNum = NOS_L * IdNumProto;
			VCNum = NOS_L * VCNumProto; // total number of all components of all vertices 
			// Allocate memory for all CANES (spins) 
			vertices	= (float  *)malloc(VCNum * sizeof( float  ));
			normals 	= NULL;
			colors 		= (float  *)malloc(VCNum * sizeof( float  ));
			indices		= (GLuint *)malloc(IdNum * sizeof( GLuint ));
	}
}

void ReallocateArrayDrawing_H()
{
	free(vertexProto_H); free(normalProto_H); free(indicesProto_H);	
	free(vertices_H); free(normals_H); free(colors_H); free(indices_H);			
	// arrowFaces - number of arrow faces
	int ElNum_H = 5*arrowFaces_H-4; // number of triangles per arrow
	IdNum_H = 3*ElNum_H; // number of indixes per arrow
	VCNum_H = 3*(2*(1+arrowFaces_H)-2+4*arrowFaces_H+3*arrowFaces_H);
	// Allocate memory for H arrow prototype  
	vertexProto_H  	= (float  *)malloc(VCNum_H * sizeof( float  ));
	normalProto_H  	= (float  *)malloc(VCNum_H * sizeof( float  ));
	indicesProto_H 	= (GLuint *)malloc(IdNum_H * sizeof( GLuint ));	
	// Allocate memory for H arrow 
	vertices_H		= (float  *)malloc(VCNum_H * sizeof( float  ));
	normals_H 		= (float  *)malloc(VCNum_H * sizeof( float  ));
	colors_H 		= (float  *)malloc(VCNum_H * sizeof( float  ));
	indices_H		= (GLuint *)malloc(IdNum_H * sizeof( GLuint ));				
}

void ReallocateArrayDrawing_BOX()
{
	free(vertices_BOX); free(normals_BOX); free(colors_BOX); free(indices_BOX);			
	int ElNum_BOX = 6*2*12; // number of triangles 
	IdNum_BOX = 3*ElNum_BOX; // number of indixes
	VCNum_BOX = 6*4*3*12;
	vertices_BOX	= (float  *)malloc(VCNum_BOX * sizeof( float  ));
	normals_BOX 	= (float  *)malloc(VCNum_BOX * sizeof( float  ));
	colors_BOX 		= (float  *)malloc(VCNum_BOX * sizeof( float  ));
	indices_BOX		= (GLuint *)malloc(IdNum_BOX * sizeof( GLuint ));				
}

void ReallocateArrayDrawing_BASIS()
{
	free(vertices_BASIS); free(normals_BASIS); free(colors_BASIS); free(indices_BASIS);			
	int ElNum_BASIS = 6*2*4; // number of triangles 
	IdNum_BASIS = 3*ElNum_BASIS; // number of indixes
	VCNum_BASIS = 6*4*3*4;
	vertices_BASIS	= (float  *)malloc(VCNum_BASIS * sizeof( float  ));
	normals_BASIS	= (float  *)malloc(VCNum_BASIS * sizeof( float  ));
	colors_BASIS	= (float  *)malloc(VCNum_BASIS * sizeof( float  ));
	indices_BASIS	= (GLuint *)malloc(IdNum_BASIS * sizeof( GLuint ));				
}

void ReallocateArrayDrawing_PBC()
{
	free(vertices_PBC_A); free(normals_PBC_A); free(colors_PBC_A); free(indices_PBC_A);		
	free(vertices_PBC_B); free(normals_PBC_B); free(colors_PBC_B); free(indices_PBC_B);		
	free(vertices_PBC_C); free(normals_PBC_C); free(colors_PBC_C); free(indices_PBC_C);			
	int ElNum_PBC = 6*2*16; // number of triangles 
	IdNum_PBC = 3*ElNum_PBC; // number of indixes 
	VCNum_PBC = 6*4*3*16;
	vertices_PBC_A	= (float  *)malloc(VCNum_PBC * sizeof( float  ));
	vertices_PBC_B	= (float  *)malloc(VCNum_PBC * sizeof( float  ));
	vertices_PBC_C	= (float  *)malloc(VCNum_PBC * sizeof( float  ));

	normals_PBC_A 	= (float  *)malloc(VCNum_PBC * sizeof( float  ));
	normals_PBC_B 	= (float  *)malloc(VCNum_PBC * sizeof( float  ));
	normals_PBC_C 	= (float  *)malloc(VCNum_PBC * sizeof( float  ));

	colors_PBC_A 	= (float  *)malloc(VCNum_PBC * sizeof( float  ));
	colors_PBC_B 	= (float  *)malloc(VCNum_PBC * sizeof( float  ));
	colors_PBC_C 	= (float  *)malloc(VCNum_PBC * sizeof( float  ));

	indices_PBC_A	= (GLuint *)malloc(VCNum_PBC * sizeof( GLuint ));
	indices_PBC_B	= (GLuint *)malloc(VCNum_PBC * sizeof( GLuint ));
	indices_PBC_C	= (GLuint *)malloc(VCNum_PBC * sizeof( GLuint ));					
}

void UpdatePrototypeVerNorInd(float * V, float * N, GLuint * I, int faces, int mode, int style)//faces = arrowFaces
{
	int   i, j;
	float H = 1.00f;	// total length
	float dF=D2R*360.f/float(faces);
	float R = H*0.20f;	// big radius
	float r = H*0.06f;	// small radius 
	float h = H*0.40f;	// H - head
	float P = H*Pivot;	// pivot - central point on which the arrow turns (shift up in Z from 0 to 1)
	float cosF0, cosF1, cosF2, sinF0, sinF1, sinF2;
	float tmp0[3], tmp1[3], tmp2[3], tmp3[3];
	float v[8][3]={	{0,0,0},
					{abc[0][0],abc[0][1],abc[0][2]},
					{abc[1][0],abc[1][1],abc[1][2]},
					{abc[0][0]+abc[1][0],abc[0][1]+abc[1][1],abc[0][2]+abc[1][2]},
					{abc[2][0],abc[2][1],abc[2][2]},
					{abc[2][0]+abc[0][0],abc[2][1]+abc[0][1],abc[2][2]+abc[0][2]},
					{abc[2][0]+abc[1][0],abc[2][1]+abc[1][1],abc[2][2]+abc[1][2]},
					{abc[2][0]+abc[0][0]+abc[1][0],
						abc[2][1]+abc[0][1]+abc[1][1],
							abc[2][2]+abc[0][2]+abc[1][2]},
				};
	if (style==1)
	{
		R = H*0.1f;	// big radius
		r = H*0.03f;	// small radius 
		h = H*0.65f;	// H - head
	}

	//float h=H-H/GoldenRatio;	//H - head
	//float R=h/GoldenRatio;	//big radius
	//float r=R-R/GoldenRatio;	//small radius 
	// Arrow //////////////////////////////////// Cane//////////////////////////////////////////////
	//          Top____________                //                                                 //
	//     H    / \           ^                //             O v1(x1,y1,z1)                      // 
	//     e   /   \          |                //            /                                    //
	//     a  /     \         |                //           /                                     //
	//     d /__R/2__\____    H    Z           //          / GL_LINE                              //
	//       T|     |    ^    |    ^           //         /                                       //
	//	     a|     |    h    |    |           //        /                                        //
	//       i|     |    |    |    |           //       /                                         //
	//       l|_____|____v____v_   o --- > X   //      o v0(x0,y0,z0)                             //
	//          r/2                            //                                                 //
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

/*			i++; V[i] = 0.f;   N[i] = 0;	
			i++; V[i] = 0.f;   N[i] = 0;	
			i++; V[i] = H-P;   N[i] = 0;	

			i++; V[i] = cosF2; N[i] = cosF2;	
			i++; V[i] = sinF2; N[i] = sinF2;	
			i++; V[i] = h-P;   N[i] = 0;

			i++; V[i] = cosF1; N[i] = cosF1;
			i++; V[i] = sinF1; N[i] = sinF1;	 
			i++; V[i] = h-P;   N[i] = 0;*/	

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
			cosF0 = 1*cos(dF*n+dF/2.f); sinF0 = 1*sin(dF*n+dF/2.f);	// for normals
			//cosF0 = 1*cos(dF*n); sinF0 = 1*sin(dF*n);	// for normals
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

			i++; V[i] = tmp0[0] = 0.f;		
			i++; V[i] = tmp0[1] = 0.f;		
			i++; V[i] = tmp0[2] = H-P;		

			i++; V[i] = tmp2[0] = cosF2;	
			i++; V[i] = tmp2[1] = sinF2;	
			i++; V[i] = tmp2[2] = - P;	

			i++; V[i] = tmp1[0] = cosF1;	
			i++; V[i] = tmp1[1] = sinF1;	 
			i++; V[i] = tmp1[2] = - P;		

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

			V[++i] = v[0][0]; N[i] = tmp1[0]; V[i+12] = v[4][0]; N[i+12] = tmp1[0];
			V[++i] = v[0][0]; N[i] = tmp1[1]; V[i+12] = v[4][1]; N[i+12] = tmp1[1];
			V[++i] = v[0][0]; N[i] = tmp1[2]; V[i+12] = v[4][2]; N[i+12] = tmp1[2];

			V[++i] = v[1][0]; N[i] = tmp1[0]; V[i+12] = v[6][0]; N[i+12] = tmp1[0];
			V[++i] = v[1][1]; N[i] = tmp1[1]; V[i+12] = v[6][1]; N[i+12] = tmp1[1];
			V[++i] = v[1][2]; N[i] = tmp1[2]; V[i+12] = v[6][2]; N[i+12] = tmp1[2];

			V[++i] = v[2][0]; N[i] = tmp1[0]; V[i+12] = v[5][0]; N[i+12] = tmp1[0];
			V[++i] = v[2][1]; N[i] = tmp1[1]; V[i+12] = v[5][1]; N[i+12] = tmp1[1];
			V[++i] = v[2][2]; N[i] = tmp1[2]; V[i+12] = v[5][2]; N[i+12] = tmp1[2];
			
			V[++i] = v[3][0]; N[i] = tmp1[0]; V[i+12] = v[7][0]; N[i+12] = tmp1[0];
			V[++i] = v[3][1]; N[i] = tmp1[1]; V[i+12] = v[7][1]; N[i+12] = tmp1[1];
			V[++i] = v[3][2]; N[i] = tmp1[2]; V[i+12] = v[7][2]; N[i+12] = tmp1[2];	
            i+=12;
			I[++j] = 0; I[++j] = 1; I[++j] = 2; // first  triangle bottom
			I[++j] = 2; I[++j] = 1; I[++j] = 3; // second triangle bottom
			I[++j] = 4; I[++j] = 5; I[++j] = 6; // third  triangle top
			I[++j] = 6; I[++j] = 5; I[++j] = 7; // fourth triangle top

            // right/left
            V[++i] = v[0][0]; N[i] = tmp2[0]; V[i+12] = v[3][0]; N[i+12] = tmp2[0];
            V[++i] = v[0][0]; N[i] = tmp2[1]; V[i+12] = v[3][1]; N[i+12] = tmp2[1];
            V[++i] = v[0][0]; N[i] = tmp2[2]; V[i+12] = v[3][2]; N[i+12] = tmp2[2];

            V[++i] = v[4][0]; N[i] = tmp2[0]; V[i+12] = v[7][0]; N[i+12] = tmp2[0];
            V[++i] = v[4][1]; N[i] = tmp2[1]; V[i+12] = v[7][1]; N[i+12] = tmp2[1];
            V[++i] = v[4][2]; N[i] = tmp2[2]; V[i+12] = v[7][2]; N[i+12] = tmp2[2];

            V[++i] = v[1][0]; N[i] = tmp2[0]; V[i+12] = v[2][0]; N[i+12] = tmp2[0];
            V[++i] = v[1][1]; N[i] = tmp2[1]; V[i+12] = v[2][1]; N[i+12] = tmp2[1];
            V[++i] = v[1][2]; N[i] = tmp2[2]; V[i+12] = v[2][2]; N[i+12] = tmp2[2];
            
            V[++i] = v[5][0]; N[i] = tmp2[0]; V[i+12] = v[6][0]; N[i+12] = tmp2[0];
            V[++i] = v[5][1]; N[i] = tmp2[1]; V[i+12] = v[6][1]; N[i+12] = tmp2[1];
            V[++i] = v[5][2]; N[i] = tmp2[2]; V[i+12] = v[6][2]; N[i+12] = tmp2[2];  
            i+=12;
			I[++j] = 8+0; I[++j] = 8+1; I[++j] = 8+2; // first  triangle 
			I[++j] = 8+2; I[++j] = 8+1; I[++j] = 8+3; // second triangle 
			I[++j] = 8+4; I[++j] = 8+5; I[++j] = 8+6; // third  triangle 
			I[++j] = 8+6; I[++j] = 8+5; I[++j] = 8+7; // fourth triangle 

            // front/back
            V[++i] = v[0][0]; N[i] = tmp3[0]; V[i+12] = v[1][0]; N[i+12] = tmp3[0];
            V[++i] = v[0][0]; N[i] = tmp3[1]; V[i+12] = v[1][1]; N[i+12] = tmp3[1];
            V[++i] = v[0][0]; N[i] = tmp3[2]; V[i+12] = v[1][2]; N[i+12] = tmp3[2];

            V[++i] = v[2][0]; N[i] = tmp3[0]; V[i+12] = v[5][0]; N[i+12] = tmp3[0];
            V[++i] = v[2][1]; N[i] = tmp3[1]; V[i+12] = v[5][1]; N[i+12] = tmp3[1];
            V[++i] = v[2][2]; N[i] = tmp3[2]; V[i+12] = v[5][2]; N[i+12] = tmp3[2];

            V[++i] = v[4][0]; N[i] = tmp3[0]; V[i+12] = v[3][0]; N[i+12] = tmp3[0];
            V[++i] = v[4][1]; N[i] = tmp3[1]; V[i+12] = v[3][1]; N[i+12] = tmp3[1];
            V[++i] = v[4][2]; N[i] = tmp3[2]; V[i+12] = v[3][2]; N[i+12] = tmp3[2];
            
            V[++i] = v[6][0]; N[i] = tmp3[0]; V[i+12] = v[7][0]; N[i+12] = tmp3[0];
            V[++i] = v[6][1]; N[i] = tmp3[1]; V[i+12] = v[7][1]; N[i+12] = tmp3[1];
            V[++i] = v[6][2]; N[i] = tmp3[2]; V[i+12] = v[7][2]; N[i+12] = tmp3[2];  

			I[++j] = 2*8+0; I[++j] = 2*8+1; I[++j] = 2*8+2; // first  triangle 
			I[++j] = 2*8+2; I[++j] = 2*8+1; I[++j] = 2*8+3; // second triangle 
			I[++j] = 2*8+4; I[++j] = 2*8+5; I[++j] = 2*8+6; // third  triangle 
			I[++j] = 2*8+6; I[++j] = 2*8+5; I[++j] = 2*8+7; // fourth triangle 
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


void UpdateIndices(GLuint * Iinp , int Kinp, GLuint * Iout, int Kout, int VerN)
{	// VerN number of vertesiec components per one prototype arrow
	int i,j;
	int NOS_L=0;
	int NOB_L=0;
	switch( WhichSliceMode)
	{
		case A_AXIS:
		NOS_L=NOS_AL * (1+A_layer_max - A_layer_min);
		if (A_layer_max - A_layer_min<2){
			NOB_L=NOB_AL * (1+A_layer_max - A_layer_min);
		}else{
			NOB_L=2*NOB_AL+2*(A_layer_max - A_layer_min-1)*(uABC[1]-1+uABC[2]-1);
		}
		break;
		case B_AXIS:
		NOS_L=NOS_BL * (1+B_layer_max - B_layer_min);
		if (B_layer_max - B_layer_min<2){
			NOB_L=NOB_BL * (1+B_layer_max - B_layer_min);
		}else{
			NOB_L=2*NOB_BL + 2*(B_layer_max - B_layer_min-1)*(uABC[0]-1+uABC[2]-1);
		}
		break;
		case C_AXIS:
		NOS_L=NOS_CL * (1+C_layer_max - C_layer_min);
		if (C_layer_max - C_layer_min<2){
			NOB_L=NOB_CL * (1+C_layer_max - C_layer_min);
		}else{
			NOB_L=2*NOB_CL + 2*(C_layer_max - C_layer_min-1)*(uABC[0]-1+uABC[1]-1);
		}
		break;
		case FILTER:
		NOS_L=NOS;
		NOB_L=NOS/AtomsPerBlock;
		break;
	}

	if (WhichVectorMode==BOX1){
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



void UpdateVerticesNormalsColors (float * Vinp, float * Ninp, int Kinp, 
							float * Vout, float * Nout, float * Cout, int Kout, 
							float * Px, float * Py, float * Pz,
							double * Sx, double * Sy, double * Sz, int mode)
{
	//float tmpV1[3], tmpV2[3], tmpV3[3], RGB[3], U, A;
	float S[3], RGB[3], U, A;
	float vlength;
	int i,j,n,N;
	int anini=0;
	int anfin=uABC[0];
	int bnini=0;
	int bnfin=uABC[1];
	int cnini=0;
	int cnfin=uABC[2];

	switch( WhichSliceMode)
	{
		case A_AXIS:
			anini=A_layer_min-1;
	        anfin=A_layer_max;
		break;
		case B_AXIS:
			bnini=B_layer_min-1;
	        bnfin=B_layer_max;
		break;
		case C_AXIS:
			cnini=C_layer_min-1;
	        cnfin=C_layer_max;
		break;
		default:
		break;
	}
	j=-1;
	switch (WhichVectorMode)
	{
	case ARROW1:
	case CONE1:
		if (WhichSliceMode==FILTER)
			{
			N_filter=0;
			for (int an = 0; an<uABC[0]; an++) 
			{
			for (int bn = 0; bn<uABC[1]; bn++) 
			{
			for (int cn = 0; cn<uABC[2]; cn++) 
			{
				n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
				n = n*AtomsPerBlock;//index of the first spin in the block
				for (int atom=0; atom<AtomsPerBlock; atom++)//atom is the index of the atom in block
				{
				    N = n + atom;
					//slow version is commented but easy to read: 
					S[0] = Sx[N];
					S[1] = Sy[N];
					S[2] = Sz[N];
					int phi = atan2int( S[1], S[0] );// return integer angle phi 0 - 360
					//if (S[2]>=Sz_min1 && S[2]<=Sz_max1)
					if ( (S[2]>=Sz_min1 && S[2]<=Sz_max1) && (phi>=phi_min1 && phi<=phi_max1) )
					{
						vlength = Unitf(S,S);
						N_filter++;
				        HSVtoRGB(S, RGB, InvertValue, InvertHue);
						j++;
						for (int k=0; k<Kinp/3; k++) // k runs over vertices of the arrow/cone 
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

							Vout[i+0] =(-S[1]*A + Vinp[3*k+0]*S[2] + S[0]*Vinp[3*k+2]			)*Scale*vlength + Px[N];
							Vout[i+1] =( S[0]*A + Vinp[3*k+1]*S[2] + S[1]*Vinp[3*k+2]			)*Scale*vlength + Py[N];
							Vout[i+2] =( Vinp[3*k+2]*S[2] - (S[0]*Vinp[3*k+0]+S[1]*Vinp[3*k+1])	)*Scale*vlength + Pz[N];	

							//slow version is commented but easy to read:
							// tmpV1[0] = Ninp[3*k+0];
							// tmpV1[1] = Ninp[3*k+1];
							// tmpV1[2] = Ninp[3*k+2];
							//NewBasisCartesian(tmpV1, S, tmpV3);	// to find vertices normals w.r.t basis based on Sx,Sy,Sz
							// Nout[i+0] = tmpV3[0];		// x-component of vertex normal
							// Nout[i+1] = tmpV3[1];		// y-component of vertex normal
							// Nout[i+2] = tmpV3[2];		// z-component of vertex normal

							A = (-S[1]*Vinp[3*k+0] + S[0]*Vinp[3*k+1])*(1. - S[2])*U; 

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
		}else{//IF SLICING MODE = A-AXIS, B-AXIS, C-AXIS 
			for (int an = anini; an<anfin; an++) 
			{
			for (int bn = bnini; bn<bnfin; bn++) 
			{
			for (int cn = cnini; cn<cnfin; cn++) 
			{
				n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
				n = n*AtomsPerBlock;//index of the first spin in the block
				for (int atom=0; atom<AtomsPerBlock; atom++)//atom is the index of the atom in block
				{
				    N = n + atom;
					//slow version is commented but easy to read: 
					S[0] = Sx[N];
					S[1] = Sy[N];
					S[2] = Sz[N];
					vlength = Unitf(S,S);
			        HSVtoRGB( S, RGB, InvertValue, InvertHue);
					j++;
					if (S[2]==-1){
					for (int k=0; k<Kinp/3; k++){// k runs over vertices 
							i = j*Kinp + 3*k;
							Vout[i+0] = (-Vinp[3*k+0])*Scale*vlength + Px[N];
							Vout[i+1] = ( Vinp[3*k+1])*Scale*vlength + Py[N];
							Vout[i+2] = (-Vinp[3*k+2])*Scale*vlength + Pz[N];	

							Nout[i+0] = -Ninp[3*k+0];
							Nout[i+1] =  Ninp[3*k+1];
							Nout[i+2] = -Ninp[3*k+2];

							Cout[i+0] = RGB[0];			// x-component of vertex normal
							Cout[i+1] = RGB[1];			// y-component of vertex normal
							Cout[i+2] = RGB[2];			// z-component of vertex normal
						}
					}else{
					for (int k=0; k<Kinp/3; k++) // k runs over vertices of the arrow/cone 
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

							Vout[i+0] =(-S[1]*A + Vinp[3*k+0]*S[2] + S[0]*Vinp[3*k+2]			)*Scale*vlength + Px[N];
							Vout[i+1] =( S[0]*A + Vinp[3*k+1]*S[2] + S[1]*Vinp[3*k+2]			)*Scale*vlength + Py[N];
							Vout[i+2] =( Vinp[3*k+2]*S[2] - (S[0]*Vinp[3*k+0]+S[1]*Vinp[3*k+1])	)*Scale*vlength + Pz[N];

							//slow version is commented but easy to read:
							// tmpV1[0] = Ninp[3*k+0];
							// tmpV1[1] = Ninp[3*k+1];
							// tmpV1[2] = Ninp[3*k+2];
							//NewBasisCartesian(tmpV1, S, tmpV3);	// to find vertices normals w.r.t basis based on Sx,Sy,Sz
							// Nout[i+0] = tmpV3[0];		// x-component of vertex normal
							// Nout[i+1] = tmpV3[1];		// y-component of vertex normal
							// Nout[i+2] = tmpV3[2];		// z-component of vertex normal

							A = (-S[1]*Vinp[3*k+0] + S[0]*Vinp[3*k+1])*(1. - S[2])*U;

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
		if (WhichSliceMode==FILTER)
		{
			N_filter=0;
			for (int an = 0; an<uABC[0]; an++) {
			for (int bn = 0; bn<uABC[1]; bn++) {
			for (int cn = 0; cn<uABC[2]; cn++) {	
				S[0] = 0;
				S[1] = 0;
				S[2] = 0;

				n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
				for (int atom=0; atom<AtomsPerBlock; atom++)//atom is the index of the atom in block
				{
				    N = atom + n*AtomsPerBlock;//n*AtomsPerBlock=index of the first spin in the block;
	 
					S[0]+= Sx[N];
					S[1]+= Sy[N];
					S[2]+= Sz[N];
			    }
				int phi = atan2int( S[1], S[0] );// return integer angle phi 0 - 360
				//if (S[2]>=Sz_min1 && S[2]<=Sz_max1)
				if ( (S[2]>=Sz_min1 && S[2]<=Sz_max1) && (phi>=phi_min1 && phi<=phi_max1) )
				{
					N_filter++;
				    (void)Unitf(S,S);
				    HSVtoRGB( S, RGB, InvertValue, InvertHue);

					j++;
					for (int k=0; k<Kinp/3; k++) // k runs over vertices of the box 
					{
						i = j*Kinp + 3*k;	// vertex index

						Vout[i+0] = Vinp[3*k+0] + BPx[n];
						Vout[i+1] = Vinp[3*k+1] + BPy[n];
						Vout[i+2] = Vinp[3*k+2] + BPz[n];	

						Nout[i+0] = Ninp[3*k+0];
						Nout[i+1] = Ninp[3*k+1];
						Nout[i+2] = Ninp[3*k+2];

						Cout[i+0] = RGB[0];			
						Cout[i+1] = RGB[1];			
						Cout[i+2] = RGB[2];	
					}
				}//if	
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

					n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
					for (int atom=0; atom<AtomsPerBlock; atom++)//atom is the index of the atom in block
					{
					    N = atom + n*AtomsPerBlock;//n*AtomsPerBlock=index of the first spin in the block;
		 
						S[0]+= Sx[N];
						S[1]+= Sy[N];
						S[2]+= Sz[N];
				    }
				    (void)Unitf(S,S);
				    HSVtoRGB( S, RGB, InvertValue, InvertHue);

					j++;
					for (int k=0; k<Kinp/3; k++) // k runs over vertices of the box 
					{
						i = j*Kinp + 3*k;	// vertex index

						Vout[i+0] = Vinp[3*k+0] + BPx[n];
						Vout[i+1] = Vinp[3*k+1] + BPy[n];
						Vout[i+2] = Vinp[3*k+2] + BPz[n];	

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
		if (WhichSliceMode==FILTER)
		{
			N_filter=0;
			for (int an = 0; an<uABC[0]; an++) // n runs over spins 
			{
			for (int bn = 0; bn<uABC[1]; bn++) // n runs over spins 
			{
			for (int cn = 0; cn<uABC[2]; cn++) // n runs over spins 
			{
				n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
				n = n*AtomsPerBlock;//index of the first spin in the block
				for (int atom=0; atom<AtomsPerBlock; atom++)//atom is the index of the atom in block
				{
					S[0] = Sx[n+atom];
					S[1] = Sy[n+atom];
					S[2] = Sz[n+atom];
					vlength = Unitf(S,S);
					int phi = atan2int( S[1], S[0] );// return integer angle phi 0 - 360
					//if (S[2]>=Sz_min1 && S[2]<=Sz_max1)
					if ( (S[2]>=Sz_min1 && S[2]<=Sz_max1) && (phi>=phi_min1 && phi<=phi_max1) )
					{
						N_filter++;
				        HSVtoRGB( S, RGB, InvertValue, InvertHue);
				        j++;
						i = j*Kinp;			// index of first cane vertex 
						Vout[i+0] = Px[n+atom];	// new x-component of vertex + translation
						Vout[i+1] = Py[n+atom];	// new y-component of vertex + translation
						Vout[i+2] = Pz[n+atom];	// new z-component of vertex + translation
						Cout[i+0] = RGB[0];	// x-component of vertex normal
						Cout[i+1] = RGB[1];	// y-component of vertex normal
						Cout[i+2] = RGB[2];	// z-component of vertex normal
				    }
				}	
			}
			}
			}
		}else{
			for (int an = anini; an<anfin; an++) // n runs over spins 
			{
			for (int bn = bnini; bn<bnfin; bn++) // n runs over spins 
			{
			for (int cn = cnini; cn<cnfin; cn++) // n runs over spins 
			{
				n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
				n = n*AtomsPerBlock;//index of the first spin in the block
				for (int atom=0; atom<AtomsPerBlock; atom++)//atom is the index of the atom in block
				{
					S[0] = Sx[n+atom];
					S[1] = Sy[n+atom];
					S[2] = Sz[n+atom];
					vlength = Unitf(S,S);
			        HSVtoRGB( S, RGB, InvertValue, InvertHue);
			        j++;
					i = j*Kinp;			// index of first cane vertex 
					Vout[i+0] = Px[n+atom];	// new x-component of vertex + translation
					Vout[i+1] = Py[n+atom];	// new y-component of vertex + translation
					Vout[i+2] = Pz[n+atom];	// new z-component of vertex + translation
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
		if (WhichSliceMode==FILTER)
		{
			N_filter=0;
			for (int an = 0; an<uABC[0]; an++) // an runs over slice along A_AXIS
			{
			for (int bn = 0; bn<uABC[1]; bn++) // bn runs over slice along B_AXIS
			{
			for (int cn = cnini; cn<uABC[2]; cn++) // cn runs over slice along C_AXIS
			{
				n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
				n = n*AtomsPerBlock;//index of the first spin in the block
				for (int atom=0; atom<AtomsPerBlock; atom++)//atom is the index of the atom in block
				{
					S[0] = Sx[n+atom];
					S[1] = Sy[n+atom];
					S[2] = Sz[n+atom];
					vlength = Unitf(S,S);
					int phi = atan2int( S[1], S[0] );// return integer angle phi 0 - 360
					//if (S[2]>=Sz_min1 && S[2]<=Sz_max1)
					if ( (S[2]>=Sz_min1 && S[2]<=Sz_max1) && (phi>=phi_min1 && phi<=phi_max1) )
					{
						N_filter++;
				        HSVtoRGB( S, RGB, InvertValue, InvertHue);
				        j++;
						//i = (n-nini)*Kinp;							// index of ferst cane vertex 
						i = j*Kinp;
						Vout[i+0] = S[0]*(1-Pivot)*Scale*vlength + Px[n+atom];	// new x-component of vertex + translation
						Vout[i+1] = S[1]*(1-Pivot)*Scale*vlength + Py[n+atom];	// new y-component of vertex + translation
						Vout[i+2] = S[2]*(1-Pivot)*Scale*vlength + Pz[n+atom];	// new z-component of vertex + translation
						//i = n*Kinp/3*4;		// colors contains 4 floats
						Cout[i+0] = RGB[0];					// x-component of vertex normal
						Cout[i+1] = RGB[1];					// y-component of vertex normal
						Cout[i+2] = RGB[2];					// z-component of vertex normal
						//Cout[i+3] = 1.f;
						//printf( "|V1=%f,%f,%f \n",Vout[i+0],Vout[i+1],Vout[i+2]);
						//i = (n-nini)*Kinp + 3*1;					// index of ferst cane vertex 
						i = j*Kinp+ 3*1;
						Vout[i+0] = -S[0]*(Pivot)*Scale*vlength + Px[n+atom];		// new x-component of vertex + translation
						Vout[i+1] = -S[1]*(Pivot)*Scale*vlength + Py[n+atom];		// new y-component of vertex + translation
						Vout[i+2] = -S[2]*(Pivot)*Scale*vlength + Pz[n+atom];		// new z-component of vertex + translation
						//i = n*Kinp/3*4+4;			// colors contains 4 floats
						Cout[i+0] = RGB[0];					// x-component of vertex normal
						Cout[i+1] = RGB[1];					// y-component of vertex normal
						Cout[i+2] = RGB[2];					// z-component of vertex normal
						//Cout[i+3] = 1.f;
						//printf( "|V1=%f,%f,%f \n",Vout[i+0],Vout[i+1],Vout[i+2]);
					}
				}
			}
			}
			}
		}else{
			for (int an = anini; an<anfin; an++) // an runs over slice along A_AXIS
			{
			for (int bn = bnini; bn<bnfin; bn++) // bn runs over slice along B_AXIS
			{
			for (int cn = cnini; cn<cnfin; cn++) // cn runs over slice along C_AXIS
			{
				n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
				n = n*AtomsPerBlock;//index of the first spin in the block
				for (int atom=0; atom<AtomsPerBlock; atom++)//atom is the index of the atom in block
				{
					S[0] = Sx[n+atom];
					S[1] = Sy[n+atom];
					S[2] = Sz[n+atom];
					vlength = Unitf(S,S);
			        HSVtoRGB( S, RGB, InvertValue, InvertHue);
			        j++;
					//i = (n-nini)*Kinp;							// index of ferst cane vertex 
					i = j*Kinp;
					Vout[i+0] = S[0]*(1-Pivot)*Scale*vlength + Px[n+atom];	// new x-component of vertex + translation
					Vout[i+1] = S[1]*(1-Pivot)*Scale*vlength + Py[n+atom];	// new y-component of vertex + translation
					Vout[i+2] = S[2]*(1-Pivot)*Scale*vlength + Pz[n+atom];	// new z-component of vertex + translation
					//i = n*Kinp/3*4;		// colors contains 4 floats
					Cout[i+0] = RGB[0];					// x-component of vertex normal
					Cout[i+1] = RGB[1];					// y-component of vertex normal
					Cout[i+2] = RGB[2];					// z-component of vertex normal
					//Cout[i+3] = 1.f;
					//printf( "|V1=%f,%f,%f \n",Vout[i+0],Vout[i+1],Vout[i+2]);
					//i = (n-nini)*Kinp + 3*1;					// index of ferst cane vertex 
					i = j*Kinp+ 3*1;
					Vout[i+0] = -S[0]*(Pivot)*Scale*vlength + Px[n+atom];		// new x-component of vertex + translation
					Vout[i+1] = -S[1]*(Pivot)*Scale*vlength + Py[n+atom];		// new y-component of vertex + translation
					Vout[i+2] = -S[2]*(Pivot)*Scale*vlength + Pz[n+atom];		// new z-component of vertex + translation
					//i = n*Kinp/3*4+4;			// colors contains 4 floats
					Cout[i+0] = RGB[0];					// x-component of vertex normal
					Cout[i+1] = RGB[1];					// y-component of vertex normal
					Cout[i+2] = RGB[2];					// z-component of vertex normal
					//Cout[i+3] = 1.f;
					//printf( "|V1=%f,%f,%f \n",Vout[i+0],Vout[i+1],Vout[i+2]);
				}
			}
			}
			}
		}
	}
}

void UpdateVerticesNormalsColors_H(float * Vinp, float * Ninp, int Kinp, 
							float * Vout, float * Nout, float * Cout,
							float Px, float Py, float Pz,
							float Sx, float Sy, float Sz)
{
	int i;
	float U,A;
	if (Sz==-1){
		for (int k=0; k<Kinp/3; k++){// k runs over vertices 
			i = 3*k;	// vertex index
			U = 1.0f/(Sx*Sx + Sy*Sy+(1e-37f)); 		
			A = (-Sy*Vinp[3*k+0] + Sx*Vinp[3*k+1])*(1. - Sz)*U; 
			Vout[i+0] = (-Vinp[i+0])*Hf*Scale_H + Px;
			Vout[i+1] = ( Vinp[i+1])*Hf*Scale_H + Py;
			Vout[i+2] = (-Vinp[i+2])*Hf*Scale_H + Pz;	

			A = (-Sy*Ninp[3*k+0] + Sx*Ninp[3*k+1])*(1. - Sz)*U; 		
			Nout[i+0] = -Ninp[i+0];
			Nout[i+1] =  Ninp[i+1];
			Nout[i+2] = -Ninp[i+2];
			if (WhichColor == BLACK){
				Cout[i+0] = 0.9;			// x-component of vertex normal
				Cout[i+1] = 0.9;			// y-component of vertex normal
				Cout[i+2] = 0.9;			// z-component of vertex normal
			}else{
				Cout[i+0] = 0.3;			// x-component of vertex normal
				Cout[i+1] = 0.3;			// y-component of vertex normal
				Cout[i+2] = 0.3;			// z-component of vertex normal			
			}
		}
	}else{	
		for (int k=0; k<Kinp/3; k++){// k runs over vertices 
			i = 3*k;	// vertex index
			U = 1.0f/(Sx*Sx + Sy*Sy+(1e-37f)); 		
			A = (-Sy*Vinp[3*k+0] + Sx*Vinp[3*k+1])*(1. - Sz)*U; 
			Vout[i+0] = (-Sy*A + Vinp[3*k+0]*Sz + Sx*Vinp[3*k+2]			)*Hf*Scale_H + Px;
			Vout[i+1] = ( Sx*A + Vinp[3*k+1]*Sz + Sy*Vinp[3*k+2]			)*Hf*Scale_H + Py;
			Vout[i+2] = ( Vinp[3*k+2]*Sz - (Sx*Vinp[3*k+0]+Sy*Vinp[3*k+1])	)*Hf*Scale_H + Pz;	

			A = (-Sy*Ninp[3*k+0] + Sx*Ninp[3*k+1])*(1. - Sz)*U; 		
			Nout[i+0] =-Sy*A + Ninp[3*k+0]*Sz + Sx*Ninp[3*k+2];
			Nout[i+1] = Sx*A + Ninp[3*k+1]*Sz + Sy*Ninp[3*k+2];
			Nout[i+2] = Ninp[3*k+2]*Sz - (Sx*Ninp[3*k+0]+Sy*Ninp[3*k+1]);
			if (WhichColor == BLACK){
				Cout[i+0] = 0.9;			// x-component of vertex normal
				Cout[i+1] = 0.9;			// y-component of vertex normal
				Cout[i+2] = 0.9;			// z-component of vertex normal
			}else{
				Cout[i+0] = 0.3;			// x-component of vertex normal
				Cout[i+1] = 0.3;			// y-component of vertex normal
				Cout[i+2] = 0.3;			// z-component of vertex normal			
			}
		}
	}
}

void
parallelepiped( float abc[][3], float tr[3], float scale1, float scale2, float scale3, 
	float color[3], int offset_index, float* V, float* N, float* C, GLuint* I)
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
	p0[0]  +=tr[0]; p0[1]  +=tr[1]; p0[2]  +=tr[2];
	p1[0]  +=tr[0]; p1[1]  +=tr[1]; p1[2]  +=tr[2];
	p2[0]  +=tr[0]; p2[1]  +=tr[1]; p2[2]  +=tr[2];
	p3[0]  +=tr[0]; p3[1]  +=tr[1]; p3[2]  +=tr[2];
	p12[0] +=tr[0]; p12[1] +=tr[1]; p12[2] +=tr[2];
	p13[0] +=tr[0]; p13[1] +=tr[1]; p13[2] +=tr[2];
	p23[0] +=tr[0]; p23[1] +=tr[1]; p23[2] +=tr[2];
	p123[0]+=tr[0]; p123[1]+=tr[1]; p123[2]+=tr[2];

	int i = offset_index*6*4*3-1;//vertex component counter
	int j = offset_index*6*2*3-1;//vertex index counter 6 sides, 2 triangles, 3 vertices per triangle
	int k = 0;//used with j

	//top vertices + normals, two triangles: p3-p123-p13, p3-p23-p123
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
	//top indices p3-p123-p13:
	k = 0;
	I[++j] = offset_index*6*4 + k * 4 + 0; //p3
	I[++j] = offset_index*6*4 + k * 4 + 1; //p123
	I[++j] = offset_index*6*4 + k * 4 + 2; //p13
    //p3-p23-p123
	I[++j] = offset_index*6*4 + k * 4 + 0; //p3
	I[++j] = offset_index*6*4 + k * 4 + 3; //p23
	I[++j] = offset_index*6*4 + k * 4 + 1; //p123

	//bottom vertices + normals, two triangles: p0-p1-p12, p0-p12-p2
	V[++i] = p0[0]; N[i] = normal3[0]; C[i]=color[0];
	V[++i] = p0[1]; N[i] = normal3[1]; C[i]=color[1];
	V[++i] = p0[2]; N[i] = normal3[2]; C[i]=color[2];

	V[++i] = p12[0]; N[i] = normal3[0]; C[i]=color[0];
	V[++i] = p12[1]; N[i] = normal3[1]; C[i]=color[1];
	V[++i] = p12[2]; N[i] = normal3[2]; C[i]=color[2];

	V[++i] = p1[0]; N[i] = normal3[0]; C[i]=color[0];
	V[++i] = p1[1]; N[i] = normal3[1]; C[i]=color[1];
	V[++i] = p1[2]; N[i] = normal3[2]; C[i]=color[2];

	V[++i] = p2[0]; N[i] = normal3[0]; C[i]=color[0];
	V[++i] = p2[1]; N[i] = normal3[1]; C[i]=color[1];
	V[++i] = p2[2]; N[i] = normal3[2]; C[i]=color[2];
	//bottom indices p0-p1-p12:
	k = 1;
	I[++j] = offset_index*6*4 + k * 4 + 2; //p0
	I[++j] = offset_index*6*4 + k * 4 + 1; //p1
	I[++j] = offset_index*6*4 + k * 4 + 0; //p12
    //p0-p12-p2
	I[++j] = offset_index*6*4 + k * 4 + 1; //p0
	I[++j] = offset_index*6*4 + k * 4 + 3; //p12
	I[++j] = offset_index*6*4 + k * 4 + 0; //p2

	//front vertices + normals, two triangles: p0-p3-p13, p0-p13-p1
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
	//front indices p0-p3-p13:
	k = 2;
	I[++j] = offset_index*6*4 + k * 4 + 0; //p0
	I[++j] = offset_index*6*4 + k * 4 + 1; //p3
	I[++j] = offset_index*6*4 + k * 4 + 2; //p13
    //p0-p13-p1
	I[++j] = offset_index*6*4 + k * 4 + 0; //p0
	I[++j] = offset_index*6*4 + k * 4 + 2; //p13
	I[++j] = offset_index*6*4 + k * 4 + 3; //p1

	//back vertices + normals, two triangles: p2-p12-p123, p2-p123-p23
	V[++i] = p2[0]; N[i] = normal2[0]; C[i]=color[0];
	V[++i] = p2[1]; N[i] = normal2[1]; C[i]=color[1];
	V[++i] = p2[2]; N[i] = normal2[2]; C[i]=color[2];

	V[++i] = p12[0]; N[i] = normal2[0]; C[i]=color[0];
	V[++i] = p12[1]; N[i] = normal2[1]; C[i]=color[1];
	V[++i] = p12[2]; N[i] = normal2[2]; C[i]=color[2];

	V[++i] = p123[0]; N[i] = normal2[0]; C[i]=color[0];
	V[++i] = p123[1]; N[i] = normal2[1]; C[i]=color[1];
	V[++i] = p123[2]; N[i] = normal2[2]; C[i]=color[2];

	V[++i] = p23[0]; N[i] = normal2[0]; C[i]=color[0];
	V[++i] = p23[1]; N[i] = normal2[1]; C[i]=color[1];
	V[++i] = p23[2]; N[i] = normal2[2]; C[i]=color[2];
	//back indices p2-p12-p123:
	k = 3;
	I[++j] = offset_index*6*4 + k * 4 + 0; //p2
	I[++j] = offset_index*6*4 + k * 4 + 1; //p12
	I[++j] = offset_index*6*4 + k * 4 + 2; //p123
    //p2-p123-p23
	I[++j] = offset_index*6*4 + k * 4 + 0; //p2
	I[++j] = offset_index*6*4 + k * 4 + 2; //p123
	I[++j] = offset_index*6*4 + k * 4 + 3; //p23
 
	//righ vertices + normals, two triangles: p1-p13-p123, p1-p123-p12
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
	//righ indices p1-p13-p123:
	k = 4;
	I[++j] = offset_index*6*4 + k * 4 + 0; //p1
	I[++j] = offset_index*6*4 + k * 4 + 1; //p13
	I[++j] = offset_index*6*4 + k * 4 + 2; //p123
    //p1-p123-p12
	I[++j] = offset_index*6*4 + k * 4 + 0; //p1
	I[++j] = offset_index*6*4 + k * 4 + 2; //p123
	I[++j] = offset_index*6*4 + k * 4 + 3; //p12	
 
	//left vertices + normals, two triangles: p0-p2-p23, p0-p23-p3
	V[++i] = p0[0]; N[i] = normal1[0]; C[i]=color[0];
	V[++i] = p0[1]; N[i] = normal1[1]; C[i]=color[1];
	V[++i] = p0[2]; N[i] = normal1[2]; C[i]=color[2];

	V[++i] = p2[0]; N[i] = normal1[0]; C[i]=color[0];
	V[++i] = p2[1]; N[i] = normal1[1]; C[i]=color[1];
	V[++i] = p2[2]; N[i] = normal1[2]; C[i]=color[2];

	V[++i] = p23[0]; N[i] = normal1[0]; C[i]=color[0];
	V[++i] = p23[1]; N[i] = normal1[1]; C[i]=color[1];
	V[++i] = p23[2]; N[i] = normal1[2]; C[i]=color[2];

	V[++i] = p3[0]; N[i] = normal1[0]; C[i]=color[0];
	V[++i] = p3[1]; N[i] = normal1[1]; C[i]=color[1];
	V[++i] = p3[2]; N[i] = normal1[2]; C[i]=color[2];
	//left indices p0-p2-p23:
	k = 5;
	I[++j] = offset_index*6*4 + k * 4 + 0; //p0
	I[++j] = offset_index*6*4 + k * 4 + 1; //p2
	I[++j] = offset_index*6*4 + k * 4 + 2; //p23
    //p0-p23-p3
	I[++j] = offset_index*6*4 + k * 4 + 0; //p0
	I[++j] = offset_index*6*4 + k * 4 + 2; //p23
	I[++j] = offset_index*6*4 + k * 4 + 3; //p3	
}


void UpdateVerticesNormalsColors_BOX(float * vertices, float * normals, float * colors, GLuint * indices, float box[3][3])
{
	float 	d = WireWidth;
	float 	Tr[3] = {-(box[0][0]+box[1][0]+box[2][0])/2.f,
					 -(box[0][1]+box[1][1]+box[2][1])/2.f,
					 -(box[0][2]+box[1][2]+box[2][2])/2.f 
					};
	Tr[0] -= d/2;
	Tr[1] -= d/2;
	Tr[2] -= d/2;
	float tr[3] = {Tr[0], Tr[1], Tr[2]};
	float length_a = uABC[0] * sqrt(abc[0][0]*abc[0][0] + abc[0][1]*abc[0][1] + abc[0][2]*abc[0][2])+d;
	float length_b = uABC[1] * sqrt(abc[1][0]*abc[1][0] + abc[1][1]*abc[1][1] + abc[1][2]*abc[1][2])+d;
	float length_c = uABC[2] * sqrt(abc[2][0]*abc[2][0] + abc[2][1]*abc[2][1] + abc[2][2]*abc[2][2])+d;
	float color[3]={0.7,0.7,0.7};

	parallelepiped( abc, tr, length_a, d, d, color, 0, vertices, normals, colors, indices );//(0,0,0)-->(1,0,0)
	parallelepiped( abc, tr, d, length_b, d, color, 1, vertices, normals, colors, indices );//(0,0,0)-->(0,1,0)
	parallelepiped( abc, tr, d, d, length_c, color, 2, vertices, normals, colors, indices );//(0,0,0)-->(0,0,1)

	tr[0] = Tr[0]+abc[1][0]*uABC[1]; tr[1] = Tr[1]+abc[1][1]*uABC[1]; tr[2] = Tr[2]+abc[1][2]*uABC[1];
	parallelepiped( abc, tr, length_a, d, d, color, 3, vertices, normals, colors, indices );//(0,1,0)-->(1,1,0)

	tr[0]=Tr[0]+abc[2][0]*uABC[2]; tr[1]=Tr[1]+abc[2][1]*uABC[2]; tr[2]=Tr[2]+abc[2][2]*uABC[2];
	parallelepiped( abc, tr, length_a, d, d, color, 4, vertices, normals, colors, indices );//(0,0,1)-->(0,1,1)

	tr[0]+=abc[1][0]*uABC[1]; tr[1]+=abc[1][1]*uABC[1]; tr[2]+=abc[1][2]*uABC[1];
	parallelepiped( abc, tr, length_a, d, d, color, 5, vertices, normals, colors, indices );//(1,0,1)-->(1,1,1)

	tr[0]=Tr[0]+abc[0][0]*uABC[0]; tr[1]=Tr[1]+abc[0][1]*uABC[0]; tr[2]=Tr[2]+abc[0][2]*uABC[0];
	parallelepiped( abc, tr, d, length_b, d, color, 6, vertices, normals, colors, indices );//(1,0,0)-->(1,1,0)

	tr[0]=Tr[0]+abc[2][0]*uABC[2]; tr[1]=Tr[1]+abc[2][1]*uABC[2]; tr[2]=Tr[2]+abc[2][2]*uABC[2];
	parallelepiped( abc, tr, d, length_b, d, color, 7, vertices, normals, colors, indices );//(0,0,1)-->(0,1,1)

	tr[0]+=abc[0][0]*uABC[0]; tr[1]+=abc[0][1]*uABC[0]; tr[2]+=abc[0][2]*uABC[0];
	parallelepiped( abc, tr, d, length_b, d, color, 8, vertices, normals, colors, indices );//(0,1,1)-->(1,1,1)

	tr[0]=Tr[0]+abc[0][0]*uABC[0]; tr[1]=Tr[1]+abc[0][1]*uABC[0]; tr[2]=Tr[2]+abc[0][2]*uABC[0];
	parallelepiped( abc, tr, d, d, length_c, color, 9, vertices, normals, colors, indices );//(1,0,0)-->(1,0,1)

	tr[0]=Tr[0]+abc[1][0]*uABC[1]; tr[1]=Tr[1]+abc[1][1]*uABC[1]; tr[2]=Tr[2]+abc[1][2]*uABC[1];
	parallelepiped( abc, tr, d, d, length_c, color, 10, vertices, normals, colors, indices );//(0,1,0)-->(0,1,1)

	tr[0]+=abc[0][0]*uABC[0]; tr[1]+=abc[0][1]*uABC[0]; tr[2]+=abc[0][2]*uABC[0];
	parallelepiped( abc, tr, d, d, length_c, color, 11, vertices, normals, colors, indices );//(1,1,0)-->(1,1,1)
}

void UpdateVerticesNormalsColors_BASIS(float * vertices, float * normals, float * colors, GLuint * indices, float box[3][3])
{
	float 	d = WireWidth;
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
	
	parallelepiped( cube, tr, length, d, d, colorR, 0, vertices, normals, colors, indices );//X
	parallelepiped( cube, tr, d, length, d, colorG, 1, vertices, normals, colors, indices );//Y
	parallelepiped( cube, tr, d, d, length, colorB, 2, vertices, normals, colors, indices );//Z
	d *= 1.5; tr[0]=-d/2; tr[1]=-d/2; tr[2]=-d/2;

	parallelepiped( cube, tr, d, d, d, color0, 3, vertices, normals, colors, indices );//O

}

void UpdateVerticesNormalsColors_PBC(int K0, float * vertices, float * normals, float * colors, GLuint * indices, float box[3][3])
{
	float 	d = WireWidth;
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

	tr[0] = Tr[0]+abc[K0][0]*(uABC[K0]+6*d); 
	tr[1] = Tr[1]+abc[K0][1]*(uABC[K0]+6*d);
	tr[2] = Tr[2]+abc[K0][2]*(uABC[K0]+6*d);

	parallelepiped( abc, tr, length_a, length_b, length_c, color, 0, vertices, normals, colors, indices );

	tr[0]=Tr[0]+abc[K0][0]*(uABC[K0]+16*d); 
	tr[1]=Tr[1]+abc[K0][1]*(uABC[K0]+16*d);
	tr[2]=Tr[2]+abc[K0][2]*(uABC[K0]+16*d);
    parallelepiped( abc, tr, length_a, length_b, length_c, color, 1, vertices, normals, colors, indices );

	tr[0]=Tr[0]+abc[K0][0]*(-10*d); 
	tr[1]=Tr[1]+abc[K0][1]*(-10*d);
	tr[2]=Tr[2]+abc[K0][2]*(-10*d);
	parallelepiped( abc, tr, length_a, length_b, length_c, color, 2, vertices, normals, colors, indices );

	tr[0]=Tr[0]+abc[K0][0]*(-20*d); 
	tr[1]=Tr[1]+abc[K0][1]*(-20*d);
	tr[2]=Tr[2]+abc[K0][2]*(-20*d);
	parallelepiped( abc, tr, length_a, length_b, length_c, color, 3, vertices, normals, colors, indices );

	tr[0]=Tr[0]+abc[K1][0]*uABC[K1]+abc[K0][0]*(uABC[K0]+6*d);
	tr[1]=Tr[1]+abc[K1][1]*uABC[K1]+abc[K0][1]*(uABC[K0]+6*d);
	tr[2]=Tr[2]+abc[K1][2]*uABC[K1]+abc[K0][2]*(uABC[K0]+6*d);
	parallelepiped( abc, tr, length_a, length_b, length_c, color, 4, vertices, normals, colors, indices );

	tr[0]=Tr[0]+abc[K1][0]*uABC[K1]+abc[K0][0]*(uABC[K0]+16*d);
	tr[1]=Tr[1]+abc[K1][1]*uABC[K1]+abc[K0][1]*(uABC[K0]+16*d);
	tr[2]=Tr[2]+abc[K1][2]*uABC[K1]+abc[K0][2]*(uABC[K0]+16*d);
	parallelepiped( abc, tr, length_a, length_b, length_c, color, 5, vertices, normals, colors, indices );

	tr[0]=Tr[0]+abc[K1][0]*uABC[K1]+abc[K0][0]*(-10*d);
	tr[1]=Tr[1]+abc[K1][1]*uABC[K1]+abc[K0][1]*(-10*d);
	tr[2]=Tr[2]+abc[K1][2]*uABC[K1]+abc[K0][2]*(-10*d);	
	parallelepiped( abc, tr, length_a, length_b, length_c, color, 6, vertices, normals, colors, indices );

	tr[0]=Tr[0]+abc[K1][0]*uABC[K1]+abc[K0][0]*(-20*d);
	tr[1]=Tr[1]+abc[K1][1]*uABC[K1]+abc[K0][1]*(-20*d);
	tr[2]=Tr[2]+abc[K1][2]*uABC[K1]+abc[K0][2]*(-20*d);	
	parallelepiped( abc, tr, length_a, length_b, length_c, color, 7, vertices, normals, colors, indices );

	tr[0]=Tr[0]+abc[K2][0]*uABC[K2]+abc[K0][0]*(uABC[K0]+6*d);
	tr[1]=Tr[1]+abc[K2][1]*uABC[K2]+abc[K0][1]*(uABC[K0]+6*d);
	tr[2]=Tr[2]+abc[K2][2]*uABC[K2]+abc[K0][2]*(uABC[K0]+6*d);
	parallelepiped( abc, tr, length_a, length_b, length_c, color, 8, vertices, normals, colors, indices );

	tr[0]=Tr[0]+abc[K2][0]*uABC[K2]+abc[K0][0]*(uABC[K0]+16*d);
	tr[1]=Tr[1]+abc[K2][1]*uABC[K2]+abc[K0][1]*(uABC[K0]+16*d);
	tr[2]=Tr[2]+abc[K2][2]*uABC[K2]+abc[K0][2]*(uABC[K0]+16*d);

	parallelepiped( abc, tr, length_a, length_b, length_c, color, 9, vertices, normals, colors, indices );
	tr[0]=Tr[0]+abc[K2][0]*uABC[K2]+abc[K0][0]*(-10*d);
	tr[1]=Tr[1]+abc[K2][1]*uABC[K2]+abc[K0][1]*(-10*d);
	tr[2]=Tr[2]+abc[K2][2]*uABC[K2]+abc[K0][2]*(-10*d);
	parallelepiped( abc, tr, length_a, length_b, length_c, color, 10, vertices, normals, colors, indices );
	tr[0]=Tr[0]+abc[K2][0]*uABC[K2]+abc[K0][0]*(-20*d);
	tr[1]=Tr[1]+abc[K2][1]*uABC[K2]+abc[K0][1]*(-20*d);
	tr[2]=Tr[2]+abc[K2][2]*uABC[K2]+abc[K0][2]*(-20*d);
	parallelepiped( abc, tr, length_a, length_b, length_c, color, 11, vertices, normals, colors, indices );	

	tr[0]=Tr[0]+abc[K1][0]*uABC[K1]+abc[K2][0]*uABC[K2]+abc[K0][0]*(uABC[K0]+6*d);
	tr[1]=Tr[1]+abc[K1][1]*uABC[K1]+abc[K2][1]*uABC[K2]+abc[K0][1]*(uABC[K0]+6*d);
	tr[2]=Tr[2]+abc[K1][2]*uABC[K1]+abc[K2][2]*uABC[K2]+abc[K0][2]*(uABC[K0]+6*d);
	parallelepiped( abc, tr, length_a, length_b, length_c, color, 12, vertices, normals, colors, indices );

	tr[0]=Tr[0]+abc[K1][0]*uABC[K1]+abc[K2][0]*uABC[K2]+abc[K0][0]*(uABC[K0]+16*d);
	tr[1]=Tr[1]+abc[K1][1]*uABC[K1]+abc[K2][1]*uABC[K2]+abc[K0][1]*(uABC[K0]+16*d);
	tr[2]=Tr[2]+abc[K1][2]*uABC[K1]+abc[K2][2]*uABC[K2]+abc[K0][2]*(uABC[K0]+16*d);
	parallelepiped( abc, tr, length_a, length_b, length_c, color, 13, vertices, normals, colors, indices );

	tr[0]=Tr[0]+abc[K1][0]*uABC[K1]+abc[K2][0]*uABC[K2]+abc[K0][0]*(-10*d);
	tr[1]=Tr[1]+abc[K1][1]*uABC[K1]+abc[K2][1]*uABC[K2]+abc[K0][1]*(-10*d);
	tr[2]=Tr[2]+abc[K1][2]*uABC[K1]+abc[K2][2]*uABC[K2]+abc[K0][2]*(-10*d);
	parallelepiped( abc, tr, length_a, length_b, length_c, color, 14, vertices, normals, colors, indices );

	tr[0]=Tr[0]+abc[K1][0]*uABC[K1]+abc[K2][0]*uABC[K2]+abc[K0][0]*(-20*d);
	tr[1]=Tr[1]+abc[K1][1]*uABC[K1]+abc[K2][1]*uABC[K2]+abc[K0][1]*(-20*d);
	tr[2]=Tr[2]+abc[K1][2]*uABC[K1]+abc[K2][2]*uABC[K2]+abc[K0][2]*(-20*d);
	parallelepiped( abc, tr, length_a, length_b, length_c, color, 15, vertices, normals, colors, indices );
}

void UpdateSpinPositions(float abc[][3], int uABC[3], float BD[][3], int NBD, float box[][3], float * Px, float * Py, float * Pz)
{
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
				BPx[bi] = abc[0][0]*J + abc[1][0]*K + abc[2][0]*L-Tr[0]; 
				BPy[bi] = abc[0][1]*J + abc[1][1]*K + abc[2][1]*L-Tr[1];
				BPz[bi] = abc[0][2]*J + abc[1][2]*K + abc[2][2]*L-Tr[2];				
			}
		}
	}	
}
///////////////////////////////////////////////////////////////////////////////////////////////////

void CreateNewVBO( ){
	glGenBuffers(1, &vboIdV);
	glGenBuffers(1, &vboIdN);
	glGenBuffers(1, &vboIdC);
	glGenBuffers(1, &iboIdI);
}

void CreateNewVBO_H( ){
	glGenBuffers(1, &vboIdV_H);
	glGenBuffers(1, &vboIdN_H);
	glGenBuffers(1, &vboIdC_H);
	glGenBuffers(1, &iboIdI_H);
}

void CreateNewVBO_BOX( ){
	glGenBuffers(1, &vboIdV_BOX);
	glGenBuffers(1, &vboIdN_BOX);
	glGenBuffers(1, &vboIdC_BOX);
	glGenBuffers(1, &iboIdI_BOX);
}

void CreateNewVBO_BASIS( ){
	glGenBuffers(1, &vboIdV_BASIS);
	glGenBuffers(1, &vboIdN_BASIS);
	glGenBuffers(1, &vboIdC_BASIS);
	glGenBuffers(1, &iboIdI_BASIS);
}

void CreateNewVBO_PBC( ){
	glGenBuffers(1, &vboIdV_PBC_A);
	glGenBuffers(1, &vboIdN_PBC_A);
	glGenBuffers(1, &vboIdC_PBC_A);
	glGenBuffers(1, &iboIdI_PBC_A);

	glGenBuffers(1, &vboIdV_PBC_B);
	glGenBuffers(1, &vboIdN_PBC_B);
	glGenBuffers(1, &vboIdC_PBC_B);
	glGenBuffers(1, &iboIdI_PBC_B);

	glGenBuffers(1, &vboIdV_PBC_C);
	glGenBuffers(1, &vboIdN_PBC_C);
	glGenBuffers(1, &vboIdC_PBC_C);
	glGenBuffers(1, &iboIdI_PBC_C);
}

void UpdateVBO(GLuint * V, GLuint * N, GLuint * C, GLuint * I, float * ver, float * nor, float * col, GLuint * ind)
{	//ver, nor, col and ind pointer to arrays of vertxcies components, norlamls, colors and indecies 
	if (WhichSliceMode==FILTER){
			IdNum = N_filter * IdNumProto;
			VCNum = N_filter * VCNumProto;
	}
	switch (WhichVectorMode)
	{
		case ARROW1:
		case CONE1:
		case BOX1:			
				glBindBuffer(GL_ARRAY_BUFFER, *V);
				glBufferData(GL_ARRAY_BUFFER, VCNum * sizeof(float), NULL, GL_DYNAMIC_DRAW);
				glBufferSubData(GL_ARRAY_BUFFER, 0, VCNum * sizeof(float), ver); 

				glBindBuffer(GL_ARRAY_BUFFER, *N);
				glBufferData(GL_ARRAY_BUFFER, VCNum * sizeof(float), NULL, GL_DYNAMIC_DRAW);
				glBufferSubData(GL_ARRAY_BUFFER, 0, VCNum * sizeof(float), nor);

				glBindBuffer(GL_ARRAY_BUFFER, *C);
				glBufferData(GL_ARRAY_BUFFER, VCNum * sizeof(float), NULL, GL_DYNAMIC_DRAW);
				glBufferSubData(GL_ARRAY_BUFFER, 0, VCNum * sizeof(float), col);

				glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, *I);
				glBufferData(GL_ELEMENT_ARRAY_BUFFER, 1*IdNum * sizeof(GLuint), NULL, GL_DYNAMIC_DRAW);//***?1<->2?
				glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, IdNum * sizeof(GLuint), ind);
		break;
		case uPOINT:
				glBindBuffer(GL_ARRAY_BUFFER, *V);
				glBufferData(GL_ARRAY_BUFFER, VCNum * sizeof(float), 0, GL_DYNAMIC_DRAW);
				glBufferSubData(GL_ARRAY_BUFFER, 0, VCNum * sizeof(float), ver); 
				// for cane and point modes we do not need normals	
				glBindBuffer(GL_ARRAY_BUFFER, *C);
				glBufferData(GL_ARRAY_BUFFER, VCNum * sizeof(float), 0, GL_DYNAMIC_DRAW);
				glBufferSubData(GL_ARRAY_BUFFER, 0, VCNum * sizeof(float), col);

				glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, *I);
				glBufferData(GL_ELEMENT_ARRAY_BUFFER, IdNum * sizeof(GLuint), 0, GL_DYNAMIC_DRAW);
				glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, IdNum * sizeof(GLuint), ind);
		case CANE: //Cane is default vector mode
		default:
				glBindBuffer(GL_ARRAY_BUFFER, *V);
				glBufferData(GL_ARRAY_BUFFER, VCNum * sizeof(float), 0, GL_DYNAMIC_DRAW);
				glBufferSubData(GL_ARRAY_BUFFER, 0, VCNum * sizeof(float), ver); 
				// for cane and point modes we do not need normals	
				glBindBuffer(GL_ARRAY_BUFFER, *C);
				glBufferData(GL_ARRAY_BUFFER, VCNum * sizeof(float), 0, GL_DYNAMIC_DRAW);
				glBufferSubData(GL_ARRAY_BUFFER, 0, VCNum * sizeof(float), col);

				glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, *I);
				glBufferData(GL_ELEMENT_ARRAY_BUFFER, IdNum * sizeof(GLuint), 0, GL_DYNAMIC_DRAW);
				glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, IdNum * sizeof(GLuint), ind);
	}
}

void UpdateVBO_H(GLuint * V, GLuint * N, GLuint * C, GLuint * I, float * ver, float * nor, float * col, GLuint * ind)
{	//ver, nor, col and ind pointer to arrays of vertxcies components, norlamls, colors and indecies 
	glBindBuffer(GL_ARRAY_BUFFER, *V);
	glBufferData(GL_ARRAY_BUFFER, VCNum_H * sizeof(float), NULL, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, VCNum_H* sizeof(float), ver); 

	glBindBuffer(GL_ARRAY_BUFFER, *N);
	glBufferData(GL_ARRAY_BUFFER, VCNum_H * sizeof(float), NULL, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, VCNum_H * sizeof(float), nor);

	glBindBuffer(GL_ARRAY_BUFFER, *C);
	glBufferData(GL_ARRAY_BUFFER, VCNum_H * sizeof(float), NULL, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, VCNum_H* sizeof(float), col);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, *I);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, IdNum_H * sizeof(GLuint), NULL, GL_STATIC_DRAW);
	glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, IdNum_H * sizeof(GLuint), ind);
}

void UpdateVBO_BOX(GLuint * V, GLuint * N, GLuint * C, GLuint * I, float * ver, float * nor, float * col, GLuint * ind)
{	
	glBindBuffer(GL_ARRAY_BUFFER, *V);
	glBufferData(GL_ARRAY_BUFFER, VCNum_BOX * sizeof(float), NULL, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, VCNum_BOX * sizeof(float), ver); 

	glBindBuffer(GL_ARRAY_BUFFER, *N);
	glBufferData(GL_ARRAY_BUFFER, VCNum_BOX * sizeof(float), NULL, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, VCNum_BOX * sizeof(float), nor);

	glBindBuffer(GL_ARRAY_BUFFER, *C);
	glBufferData(GL_ARRAY_BUFFER, VCNum_BOX * sizeof(float), NULL, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, VCNum_BOX * sizeof(float), col);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, *I);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, IdNum_BOX * sizeof(GLuint), NULL, GL_STATIC_DRAW);
	glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, IdNum_BOX * sizeof(GLuint), ind);
}

void UpdateVBO_BASIS(GLuint * V, GLuint * N, GLuint * C, GLuint * I, float * ver, float * nor, float * col, GLuint * ind)
{	
	glBindBuffer(GL_ARRAY_BUFFER, *V);
	glBufferData(GL_ARRAY_BUFFER, VCNum_BASIS * sizeof(float), NULL, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, VCNum_BASIS * sizeof(float), ver); 

	glBindBuffer(GL_ARRAY_BUFFER, *N);
	glBufferData(GL_ARRAY_BUFFER, VCNum_BASIS * sizeof(float), NULL, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, VCNum_BASIS * sizeof(float), nor);

	glBindBuffer(GL_ARRAY_BUFFER, *C);
	glBufferData(GL_ARRAY_BUFFER, VCNum_BASIS * sizeof(float), NULL, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, VCNum_BASIS * sizeof(float), col);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, *I);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, IdNum_BASIS * sizeof(GLuint), NULL, GL_STATIC_DRAW);
	glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, IdNum_BASIS * sizeof(GLuint), ind);
}

void UpdateVBO_PBC(GLuint * V, GLuint * N, GLuint * C, GLuint * I, float * ver, float * nor, float * col, GLuint * ind)
{	
	glBindBuffer(GL_ARRAY_BUFFER, *V);
	glBufferData(GL_ARRAY_BUFFER, VCNum_PBC * sizeof(float), NULL, GL_DYNAMIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, VCNum_PBC * sizeof(float), ver); 

	glBindBuffer(GL_ARRAY_BUFFER, *N);
	glBufferData(GL_ARRAY_BUFFER, VCNum_PBC* sizeof(float), NULL, GL_DYNAMIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, VCNum_PBC * sizeof(float), nor);

	glBindBuffer(GL_ARRAY_BUFFER, *C);
	glBufferData(GL_ARRAY_BUFFER, VCNum_PBC * sizeof(float), NULL, GL_DYNAMIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, VCNum_PBC * sizeof(float), col);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, *I);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, IdNum_PBC * sizeof(GLuint), NULL, GL_DYNAMIC_DRAW);
	glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, IdNum_PBC * sizeof(GLuint), ind);	
}


void drawVBO()
{
	switch (WhichVectorMode)
	{
		case BOX1:
		case ARROW1:
		case CONE1:
			if (Light_On) {
				glEnable(GL_LIGHTING);
			}else{
				glDisable(GL_LIGHTING);
			}
			glBindBuffer(GL_ARRAY_BUFFER, vboIdC);		glColorPointer(3, GL_FLOAT, 0, (void*)0);
			glBindBuffer(GL_ARRAY_BUFFER, vboIdN);		glNormalPointer(GL_FLOAT, 0, (void*)0);
			glBindBuffer(GL_ARRAY_BUFFER, vboIdV);		glVertexPointer(3, GL_FLOAT, 0, (void*)0);	

			glEnableClientState(GL_COLOR_ARRAY);		// enable color arrays
			glEnableClientState(GL_NORMAL_ARRAY);		// enable normal arrays
			glEnableClientState(GL_VERTEX_ARRAY);		// enable vertex arrays	

			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, iboIdI);
			glDrawElements(GL_TRIANGLES, IdNum, GL_UNSIGNED_INT, (void*)(0));

			glDisableClientState(GL_VERTEX_ARRAY);		// disable vertex arrays
			glDisableClientState(GL_NORMAL_ARRAY);		// disable normal arrays
			glDisableClientState(GL_COLOR_ARRAY);		// disable color arrays

			glBindBuffer(GL_ARRAY_BUFFER,			0);	// disable vertex arrays
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,	0);	// disable normal arrays
			glEnable(GL_LIGHTING);
		break;

		case uPOINT:
		    glDisable(GL_LIGHTING);
			glBindBuffer(GL_ARRAY_BUFFER, vboIdC);		glColorPointer(3, GL_FLOAT, 0, (void*)0);
			glBindBuffer(GL_ARRAY_BUFFER, vboIdV);		glVertexPointer(3, GL_FLOAT, 0, (void*)0);	

			glEnableClientState(GL_COLOR_ARRAY);		// enable normal arrays
			glEnableClientState(GL_VERTEX_ARRAY);		// enable vertex arrays	

			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, iboIdI);
			
			glPointSize(10.f*Scale);
			//Note, the range specifies the range of index values in the index region being rendered from.
			//glDrawElements(GL_POINTS, iNum, GL_UNSIGNED_INT, (void*)(iStart*sizeof(GLuint)));
			glDrawElements(GL_POINTS, IdNum, GL_UNSIGNED_INT, (void*)0);//draw whole VBO

			glDisableClientState(GL_VERTEX_ARRAY);		// disable vertex arrays
			glDisableClientState(GL_COLOR_ARRAY);		// disable normal arrays

			glBindBuffer(GL_ARRAY_BUFFER,			0);	// disable vertex arrays
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,	0);	// disable normal arrays
			glEnable(GL_LIGHTING);
		break;

		case CANE: // CANE is default vector mode
		default:
			glDisable(GL_LIGHTING);
			glBindBuffer(GL_ARRAY_BUFFER, vboIdC);		glColorPointer(3, GL_FLOAT, 0, (void*)0);
			glBindBuffer(GL_ARRAY_BUFFER, vboIdV);		glVertexPointer(3, GL_FLOAT, 0, (void*)0);	

			glEnableClientState(GL_COLOR_ARRAY);		// enable normal arrays
			glEnableClientState(GL_VERTEX_ARRAY);		// enable vertex arrays	

			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, iboIdI);
			
			glLineWidth(5.0f*Scale);
			//glDrawElements(GL_LINES, iNum, GL_UNSIGNED_INT, (void*)(iStart*sizeof(GLuint)));// draw a stick VCNum
			glDrawElements(GL_LINES, IdNum, GL_UNSIGNED_INT, (void*)0);// draw a stick VCNum
			glPointSize(4.0f*Scale);
			//glDrawElements(GL_POINTS, iNum, GL_UNSIGNED_INT, (void*)(iStart*sizeof(GLuint)));// draw a ball (point at the end of the stick)
			glDrawElements(GL_POINTS, IdNum, GL_UNSIGNED_INT, (void*)0);// draw a ball (point at the end of the stick)

			glDisableClientState(GL_VERTEX_ARRAY);		// disable vertex arrays
			glDisableClientState(GL_COLOR_ARRAY);		// disable normal arrays

			glBindBuffer(GL_ARRAY_BUFFER,			0);	// disable vertex arrays
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,	0);	// disable normal arrays

			glBindBuffer(GL_ARRAY_BUFFER, vboIdC);		glColorPointer(3, GL_FLOAT, 6*sizeof(float), (void*)0);
			glBindBuffer(GL_ARRAY_BUFFER, vboIdV);		glVertexPointer(3, GL_FLOAT, 6*sizeof(float), (void*)0);

			glEnableClientState(GL_COLOR_ARRAY);		// enable normal arrays
			glEnableClientState(GL_VERTEX_ARRAY);		// enable vertex arrays		

			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, iboIdI);
			
			glPointSize(10.5f*Scale);
			//glDrawElements(GL_POINTS, iNum/2, GL_UNSIGNED_INT, (void*)(iStart/2*sizeof(GLuint)));
			glDrawElements(GL_POINTS, IdNum/2, GL_UNSIGNED_INT, (void*)0);// draw a ball (point at the end of the stick)

			glDisableClientState(GL_VERTEX_ARRAY);		// disable vertex arrays
			glDisableClientState(GL_COLOR_ARRAY);		// disable normal arrays

			glBindBuffer(GL_ARRAY_BUFFER,			0);	// disable vertex arrays
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,	0);	// disable normal arrays
			glEnable(GL_LIGHTING);
	}
}

void drawVBO_H()
{
	glBindBuffer(GL_ARRAY_BUFFER, vboIdC_H);		glColorPointer(3, GL_FLOAT, 0, (void*)0);
	glBindBuffer(GL_ARRAY_BUFFER, vboIdN_H);		glNormalPointer(GL_FLOAT, 0, (void*)0);
	glBindBuffer(GL_ARRAY_BUFFER, vboIdV_H);		glVertexPointer(3, GL_FLOAT, 0, (void*)0);	

	glEnableClientState(GL_COLOR_ARRAY);		// enable color arrays
	glEnableClientState(GL_NORMAL_ARRAY);		// enable normal arrays
	glEnableClientState(GL_VERTEX_ARRAY);		// enable vertex arrays	

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, iboIdI_H);
	glDrawElements(GL_TRIANGLES, IdNum_H, GL_UNSIGNED_INT, (void*)(0));

	glDisableClientState(GL_VERTEX_ARRAY);		// disable vertex arrays
	glDisableClientState(GL_NORMAL_ARRAY);		// disable normal arrays
	glDisableClientState(GL_COLOR_ARRAY);		// disable color arrays

	glBindBuffer(GL_ARRAY_BUFFER,			0);	// disable vertex arrays
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,	0);	// disable normal arrays
}

void drawVBO_BOX()
{
	glBindBuffer(GL_ARRAY_BUFFER, vboIdC_BOX);		glColorPointer(3, GL_FLOAT, 0, (void*)0);
	glBindBuffer(GL_ARRAY_BUFFER, vboIdN_BOX);		glNormalPointer(GL_FLOAT, 0, (void*)0);
	glBindBuffer(GL_ARRAY_BUFFER, vboIdV_BOX);		glVertexPointer(3, GL_FLOAT, 0, (void*)0);	

	glEnableClientState(GL_COLOR_ARRAY );		// enable color arrays
	glEnableClientState(GL_NORMAL_ARRAY);		// enable normal arrays
	glEnableClientState(GL_VERTEX_ARRAY);		// enable vertex arrays	

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, iboIdI_BOX);
	glDrawElements(GL_TRIANGLES, IdNum_BOX, GL_UNSIGNED_INT, (void*)(0));

	glDisableClientState(GL_VERTEX_ARRAY);		// disable vertex arrays
	glDisableClientState(GL_NORMAL_ARRAY);		// disable normal arrays
	glDisableClientState(GL_COLOR_ARRAY );		// disable color arrays

	glBindBuffer(GL_ARRAY_BUFFER,			0);	// disable vertex arrays
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,	0);	// disable normal arrays
}

void drawVBO_BASIS()
{
	glBindBuffer(GL_ARRAY_BUFFER, vboIdC_BASIS);		glColorPointer(3, GL_FLOAT, 0, (void*)0);
	glBindBuffer(GL_ARRAY_BUFFER, vboIdN_BASIS);		glNormalPointer(GL_FLOAT, 0, (void*)0);
	glBindBuffer(GL_ARRAY_BUFFER, vboIdV_BASIS);		glVertexPointer(3, GL_FLOAT, 0, (void*)0);	

	glEnableClientState(GL_COLOR_ARRAY );		// enable color arrays
	glEnableClientState(GL_NORMAL_ARRAY);		// enable normal arrays
	glEnableClientState(GL_VERTEX_ARRAY);		// enable vertex arrays	

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, iboIdI_BASIS);
	glDrawElements(GL_TRIANGLES, IdNum_BASIS, GL_UNSIGNED_INT, (void*)(0));

	glDisableClientState(GL_VERTEX_ARRAY);		// disable vertex arrays
	glDisableClientState(GL_NORMAL_ARRAY);		// disable normal arrays
	glDisableClientState(GL_COLOR_ARRAY );		// disable color arrays

	glBindBuffer(GL_ARRAY_BUFFER,			0);	// disable vertex arrays
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,	0);	// disable normal arrays
}

void drawVBO_PBC_A()
{
	glBindBuffer(GL_ARRAY_BUFFER, vboIdC_PBC_A);		glColorPointer(3, GL_FLOAT, 0, (void*)0);
	glBindBuffer(GL_ARRAY_BUFFER, vboIdN_PBC_A);		glNormalPointer(GL_FLOAT, 0, (void*)0);
	glBindBuffer(GL_ARRAY_BUFFER, vboIdV_PBC_A);		glVertexPointer(3, GL_FLOAT, 0, (void*)0);	

	glEnableClientState(GL_COLOR_ARRAY);		// enable color arrays
	glEnableClientState(GL_NORMAL_ARRAY);		// enable normal arrays
	glEnableClientState(GL_VERTEX_ARRAY);		// enable vertex arrays	

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, iboIdI_PBC_A);
	glDrawElements(GL_TRIANGLES, IdNum_PBC, GL_UNSIGNED_INT, (void*)(0));

	glDisableClientState(GL_VERTEX_ARRAY);		// disable vertex arrays
	glDisableClientState(GL_NORMAL_ARRAY);		// disable normal arrays
	glDisableClientState(GL_COLOR_ARRAY);		// disable color arrays

	glBindBuffer(GL_ARRAY_BUFFER,			0);	// disable vertex arrays
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,	0);	// disable normal arrays
}

void drawVBO_PBC_B()
{
	glBindBuffer(GL_ARRAY_BUFFER, vboIdC_PBC_B);		glColorPointer(3, GL_FLOAT, 0, (void*)0);
	glBindBuffer(GL_ARRAY_BUFFER, vboIdN_PBC_B);		glNormalPointer(GL_FLOAT, 0, (void*)0);
	glBindBuffer(GL_ARRAY_BUFFER, vboIdV_PBC_B);		glVertexPointer(3, GL_FLOAT, 0, (void*)0);	

	glEnableClientState(GL_COLOR_ARRAY);		// enable color arrays
	glEnableClientState(GL_NORMAL_ARRAY);		// enable normal arrays
	glEnableClientState(GL_VERTEX_ARRAY);		// enable vertex arrays	

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, iboIdI_PBC_B);
	glDrawElements(GL_TRIANGLES, IdNum_PBC, GL_UNSIGNED_INT, (void*)(0));

	glDisableClientState(GL_VERTEX_ARRAY);		// disable vertex arrays
	glDisableClientState(GL_NORMAL_ARRAY);		// disable normal arrays
	glDisableClientState(GL_COLOR_ARRAY);		// disable color arrays

	glBindBuffer(GL_ARRAY_BUFFER,			0);	// disable vertex arrays
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,	0);	// disable normal arrays
}

void drawVBO_PBC_C()
{
	glBindBuffer(GL_ARRAY_BUFFER, vboIdC_PBC_C);		glColorPointer(3, GL_FLOAT, 0, (void*)0);
	glBindBuffer(GL_ARRAY_BUFFER, vboIdN_PBC_C);		glNormalPointer(GL_FLOAT, 0, (void*)0);
	glBindBuffer(GL_ARRAY_BUFFER, vboIdV_PBC_C);		glVertexPointer(3, GL_FLOAT, 0, (void*)0);	

	glEnableClientState(GL_COLOR_ARRAY);		// enable color arrays
	glEnableClientState(GL_NORMAL_ARRAY);		// enable normal arrays
	glEnableClientState(GL_VERTEX_ARRAY);		// enable vertex arrays	

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, iboIdI_PBC_C);
	glDrawElements(GL_TRIANGLES, IdNum_PBC, GL_UNSIGNED_INT, (void*)(0));

	glDisableClientState(GL_VERTEX_ARRAY);		// disable vertex arrays
	glDisableClientState(GL_NORMAL_ARRAY);		// disable normal arrays
	glDisableClientState(GL_COLOR_ARRAY);		// disable color arrays

	glBindBuffer(GL_ARRAY_BUFFER,			0);	// disable vertex arrays
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,	0);	// disable normal arrays
}