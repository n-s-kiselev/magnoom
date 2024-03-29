GLuint 			GLUT_window;
const char *	WINDOWTITLE = { "Magnoom v1.0" };
int				window_width	= 1400;//1024;//400;//
int				window_height	= 800;//768;//300;//
float 			asp_rat			= (float)( ((double)window_width)/((double)window_height) );
float 			asp_rat_inv		= (float)( ((double)window_height)/((double)window_width) );

float			CameraEye[3]	= { 0.0, 0.0, 150.0}; // "camera position"
float			CameraC[3]		= { 0.0, 0.0,  0.0}; // "look at point"
float			CameraUp[3]		= { 0.0, 1.0,  0.0}; // "where is up direction"

float           axisX[] = { 1, 0, 0 };
float           axisY[] = { 0, 1, 0 };
float           axisZ[] = { 0, 0, 1 };
// light parameters
// GLfloat			light_ambient[]  = {0.1, 0.1, 0.1, 1.0};
// GLfloat			light_diffuse[]  = {0.9, 0.9, 0.9, 1.0};
// GLfloat			light_specular[] = {0.0, 0.0, 0.0, 1.0};



// material parameters
GLfloat ambient[] = {0.33, 0.22, 0.03, 0.0};
GLfloat diffuse[] = {0.78, 0.57, 0.11, 1.0};
GLfloat specular[] = {0.1, 0.1, 0.08, 1.0};
GLfloat shininess = 200.0;

float			PerspSet[4]		= {60.0, asp_rat, 1, 5000}; // {Setings: Field of view vertical, apect ratio, zNear, zFar}  

typedef enum	{ORTHO,	PERSP} enProjections;	// declare new enum type for projections
enProjections	WhichProjection = PERSP; // PERSP by default 	

typedef enum	{RND, HOMO, SKYRM1, SKYRM2, SKYRM3, BOBBER_T, BOBBER_B, BOBBER_L, BOBBER_L_T, BOBBER_L_B, 
HOPFION1, SPIRAL, SKYRMION_L, GLOBULA, MultyQ, NORM} enIniState; // which mode
enIniState		WhichInitialState = RND;	// RND by default 

typedef enum 	{DEFAULT_G, CILINDER_G, SPHERE_G} enGeom; // which mode
enGeom         WhichGeometry = DEFAULT_G;

float			chSizeG = 50; // characteristic size of simulated domain in units of "a"




// which button:
typedef enum 	{XUP, YUP, ZUP, ADD_SET, RESET, QUIT, PLAY, RECORD} enButton;

// the color numbers, this order must match the radio button order
typedef enum	{WHITE, BLACK, RED, GREEN, BLUE, MANUAL} enColors;
int my_background_color[3] = { 55, 55, 155 };
enColors		WhichBackgroundColor = MANUAL;	// index into Colors[ ]

int temp_color[3] = { 55, 55, 155 };

// the color definitions, this order must match the radio button order
GLfloat 	BackgroundColors[6][3] = 
{
				{ 1, 1, 1 },		// white
				{ 0.1, 0.1, 0.1 },	// black
				{ 1., 0.8, 0.8 },	// red
				{ 0.8, 1., 0.8 },	// green
				{ 0.8, 0.8, 1. },	// blue
				{ 0.3, 0.3, 0.3 },	// manual
};

typedef enum	{ARROW1 
				,CONE1 
				,CANE
				,uPOINT 
				,BOX1
				} enVectorMode; // which mode
enVectorMode	WhichVectorMode	= BOX1;	// CANE by default 






// non-constant global variables:
int				ActiveButton;	// current mous button that is down
int				Xmouse, Ymouse;	// mouse values

float 			RotAxis[3]={0,0,0};//NSK
// Shape orientation (stored as a quaternion)
float 			q_Rotation[] = { 0.0f, 0.0f, 0.0f, 1.0f };


float			Rot[3]={0,0,0};	// rotation angles in degrees
float			dRot[3]={0,0,0};// rotation angles +
float 			RotSpeed=1;// rotation speed
float			TransXYZ[3]={0,0,0};// set by glui translation widgets
float			dTransXYZ[3]={0,0,0};
float			TransSpeed=3;

int             CurrentCameraPositionBank=0;
float           CameraPosition[NumCamPosSave][7];// array which contains camera positions 
float			Scale = 1.5f;	// scaling factors for arrows [0.1-2] 
float			Pivot = 0.55f;
float 			WireWidth = 0.2;

float			Scale_H = 10;//(float)(uABC[0]+uABC[1]+uABC[2]);	// scaling factors for arrows [0.1-2] 

// // Slicing parameters
// int				A_layer_min = 1;	// which layer (along a tr. vect. ) to show max=uABC[0]
// int				B_layer_min = 1;	// which layer (along b tr. vect. ) to show max=uABC[1]
// int				C_layer_min = 1;	// which layer (along c tr. vect. ) to show max=uABC[2]
// int				A_layer_max = 1;	// which layer (along a tr. vect. ) to show max=uABC[0]
// int				B_layer_max = 1;	// which layer (along b tr. vect. ) to show max=uABC[1]
// int				C_layer_max = 1;	// which layer (along c tr. vect. ) to show max=uABC[2]

// typedef enum	{A_AXIS, B_AXIS, C_AXIS, FILTER} enSliceMode; // which mode
// enSliceMode	    WhichSliceMode	= C_AXIS;	

int   N_filter=0;

int SpinFilter1=1;
int SpinFilter2=0;
int SpinFilter3=0;
// int PhiInvert1=0;
// int PhiInvert2=0;

int theta_max1 = 90+5; //0.01;//
float  Sz_min1 = cos(theta_max1*PI/180.0);
int theta_min1 = 90-5; //0;//   
float  Sz_max1 = cos(theta_min1*PI/180.0);

int theta_max2 = 180; //0.01;//
float  Sz_min2 = cos(theta_max2*PI/180.0);
int theta_min2 = 160; //0;//   
float  Sz_max2 = cos(theta_min2*PI/180.0);

int theta_max3 = 10; //0.01;//
float  Sz_min3 = cos(theta_max3*PI/180.0);
int theta_min3 = 0; //0;//   
float  Sz_max3 = cos(theta_min3*PI/180.0);

int phi_max1=360;
int phi_min1=0;

int phi_max2=360;
int phi_min2=0;

int phi_max3=360;
int phi_min3=0;

int GreedFilter=0;
int GreedFilterInvert=0;

int GreedFilterMaxA=uABC[0]-1;// redefined in readConfigFile()
int GreedFilterMinA=0;

int GreedFilterMaxB=uABC[1]-1;// redefined in readConfigFile()
int GreedFilterMinB=0;

int GreedFilterMaxC=uABC[2]-1;// redefined in readConfigFile()
int GreedFilterMinC=0;

// Parameters for initial state 
float			chSize = 12; // characteristic size of initial state in units of "a"
float			chDir[3] = {0,1,0}; // characteristic size of initial state in units of "a"

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
TwBar *my_window; // !!!!!!
float axis[] = { 0.7f, 0.7f, 0.0f }; // initial model rotation
float angle = 0.8f;

// Light parameter
float g_LightMultiplier = 1.0f;
float g_LightDirection[] = { 0.5f, 0.5f, -0.5f };
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

GLfloat*    vertices_AC_phase  = NULL; // 
GLfloat*    colors_AC_phase    = NULL; // 
GLuint*     indices_AC_phase   = NULL; // 


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
int 		N_Multisample	= 16; // NUMBER OF MULTISAMPLE only in freeglut=Linux

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

GLuint      vboIdV_AC_phase;   // ID of VBO for vertex arrays
GLuint      vboIdC_AC_phase;   // ID of VBO for color arrays
GLuint      iboIdI_AC_phase;   // ID of IBO for index arrays

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

int         IdNum_AC_phase;
int         VCNum_AC_phase;

int			IdNum_PBC;
int			VCNum_PBC;



void			ChangeVectorMode( int );
void			ChangeColorMap( int );
void			ChangeInitialState( int );
void			Buttons( int );
void			keyboardDown( unsigned char, int, int );
void			keyboardUp( unsigned char, int, int );
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
void            drawVBO_AC_phase();
void			drawVBO_PBC_A();
void			drawVBO_PBC_B();
void			drawVBO_PBC_C();
void			idle();
void			setupTweakBar();



// return the number of seconds since the start of the program:
float ElapsedSeconds( )	{
	int ms = glutGet( GLUT_ELAPSED_TIME );	// get # of milliseconds since the start of the program
	return (float)ms / 1000.0f;				// convert it to seconds:
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

void Zup( )
{
    q_Rotation[0]=q_Rotation[1]=q_Rotation[2]=0.;
    q_Rotation[3]=1.;
    Rot[0] = Rot[1] = Rot[2] = TransXYZ[0] = TransXYZ[1] = 0.;
}

void Xup( )
{   float quat1[4];
    Zup( );
    SetQuaternionFromAxisAngle(axisZ, -D2R*90, quat1);
    MultiplyQuaternions( quat1, q_Rotation, q_Rotation);
    SetQuaternionFromAxisAngle(axisX, -D2R*90, quat1);
    MultiplyQuaternions( quat1, q_Rotation, q_Rotation);
    GetEulerFromQuaternion(q_Rotation, Rot);
    TransXYZ[0] = TransXYZ[1] = 0.;
}

void Yup( )//metka
{   float quat1[4];
    Zup( );
    SetQuaternionFromAxisAngle(axisZ, D2R*180, quat1);
    MultiplyQuaternions( quat1, q_Rotation, q_Rotation);
    SetQuaternionFromAxisAngle(axisX, -D2R*90, quat1);
    MultiplyQuaternions( quat1, q_Rotation, q_Rotation);
    GetEulerFromQuaternion(q_Rotation, Rot);
	TransXYZ[0] = TransXYZ[1] = 0.;
}



// use glut to display a string of characters using a raster font:
void DoRasterString( float x, float y, float z, char const *s )
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
void DoStrokeString( float x, float y, float z, float ht, char const *s )
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
	float mat[4*4]; // rotation matrix
    float   axisX[] = { 1, 0, 0 };
    float   axisY[] = { 0, 1, 0 };
    float   axisZ[] = { 0, 0, 1 };
    float   quat1[4];//NSK
    float   quat2[4];//NSK 
    float   quat3[4];//NSK 

	GLdouble Hight;
	float Vtemp[3];

	glutSetWindow( GLUT_window );// set which window we want to do the graphics into
	glClearColor( BackgroundColors[WhichBackgroundColor][0], BackgroundColors[WhichBackgroundColor][1], BackgroundColors[WhichBackgroundColor][2], 0. );// setup the clear values
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

	//NSK glMatrixMode( GL_MODELVIEW ); glLoadIdentity( );


	// set the eye position, look-at position, and up-vector:
	gluLookAt(	CameraEye[0],	CameraEye[1],	CameraEye[2],   // position  
				CameraC[0],		CameraC[1],		CameraC[2],     // look at
				CameraUp[0],	CameraUp[1],	CameraUp[2]);   // up

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_NORMALIZE);
	float v[4]; 
    v[0] = v[1] = v[2] = g_LightMultiplier*0.4f; v[3] = 1.0f;
    glLightfv(GL_LIGHT0, GL_AMBIENT, v);
    v[0] = v[1] = v[2] = g_LightMultiplier*0.8f; v[3] = 1.0f;
    glLightfv(GL_LIGHT0, GL_DIFFUSE, v);

	//Add ambient light
	// GLfloat lightPos0[] = {CameraEye[0], CameraEye[1], CameraEye[2], 1.0f}; //Positioned at (4, 0, 8)
	// glLightfv(GL_LIGHT0, GL_POSITION, lightPos0);
	v[0] = -g_LightDirection[0]; v[1] = -g_LightDirection[1]; v[2] = -g_LightDirection[2]; v[3] = 0.0f;
	glLightfv(GL_LIGHT0, GL_POSITION, v);
	
	//Add directed light
	// GLfloat lightColor1[] = {0.2f, 0.2f, 0.2f, 1.0f}; //Color (0.5, 0.2, 0.2)
	// glLightfv(GL_LIGHT1, GL_DIFFUSE, lightColor1);
	// v[0] = -g_LightDirection[0]; v[1] = -g_LightDirection[1]; v[2] = -g_LightDirection[2]; v[3] = 0.0f;
	// glLightfv(GL_LIGHT1, GL_POSITION, v);

    glMaterialfv(GL_FRONT, GL_AMBIENT, ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
    glMaterialf(GL_FRONT, GL_SHININESS, shininess);

	// translate the scene:
	TransXYZ[0]+=dTransXYZ[0];
	TransXYZ[1]+=dTransXYZ[1];
	TransXYZ[2]+=dTransXYZ[2];

	glTranslatef( (GLfloat)TransXYZ[0], (GLfloat)TransXYZ[1], -(GLfloat)TransXYZ[2] );
    //new value for the Euler angles due to mose rotation:
    GetEulerFromQuaternion(q_Rotation, Rot); 
    //new directions of Cartesian axes:
    RotateVectorByQuaternion(axisX, q_Rotation);
    RotateVectorByQuaternion(axisY, q_Rotation);
    RotateVectorByQuaternion(axisZ, q_Rotation);
    //adding the rotation about each axes (from keyboard), to the rotation by mouse:
    SetQuaternionFromAxisAngle(axisX, D2R*dRot[0], quat1);
    SetQuaternionFromAxisAngle(axisY, D2R*dRot[1], quat2);
    SetQuaternionFromAxisAngle(axisZ, D2R*dRot[2], quat3);
    //combining all rotations by means of quaternion multiplication:
    MultiplyQuaternions( quat2, q_Rotation, q_Rotation);
    MultiplyQuaternions( quat3, q_Rotation, q_Rotation);
    MultiplyQuaternions( quat1, q_Rotation, q_Rotation);
    //
	ConvertQuaternionToMatrix(q_Rotation, mat);
    //
    glMultMatrixf(mat);


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

	glPopMatrix();//NSK

    glDisable(GL_LIGHTING);

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
    glMatrixMode( GL_MODELVIEW);
    glLoadIdentity();
    
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
	#if !defined(__APPLE__)
	glutSetOption(GLUT_MULTISAMPLE, N_Multisample);
	#endif
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

	// InitRGB(RHue, GHue, BHue, HueMapRGB);
	// Set the GLUT callback functions
	glutDisplayFunc( Display );// DisplayFunc -- redraws the OpenGl main window
	glutReshapeFunc( Resize );// ReshapeFunc -- handles the user resizing the window
	glutKeyboardFunc(keyboardDown);// KeyboardFunc -- handles a keyboard input
	glutKeyboardUpFunc(keyboardUp);
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
	//glutTimerFunc( 40, idle, 0 );// TimerFunc -- trigger something to happen a certain time from now
	glutIdleFunc( idle );// IdleFunc -- what to do when nothing else is going on

	// Initialize AntTweakBar
    TwInit(TW_OPENGL, NULL);
    // - Send 'glutGetModifers' function pointer to AntTweakBar;
    //   required because the GLUT key event functions do not report key modifiers states.
    TwGLUTModifiersFunc(glutGetModifiers);
	setupTweakBar();

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_LIGHTING); //Enable lighting
	glEnable(GL_LIGHT0); //Enable light #0
	//glEnable(GL_LIGHT1); //Enable light #1
	//glEnable(GL_NORMALIZE); //Automatically normalize normals
	glShadeModel(GL_SMOOTH); //Enable smooth shading

	glEnable(GL_COLOR_MATERIAL);
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

void idle ()
{  
	currentTime = glutGet(GLUT_ELAPSED_TIME);
	timeInterval = currentTime - previousTime;

	if((timeInterval > 40 && Play!=0)||SpecialEvent!=0)//40ms gives approximately 25 FPS +/-1 if the engine works faster then 25 IPS
	{
		if( DATA_TRANSFER_MUTEX==TAKE_DATA || SpecialEvent!=0)
		{
			switch (SpecialEvent)
			{
				case 0:
				ChangeVectorMode(1);
				case 1:
				totalEnergy = GetTotalEnergy( 	bSx, bSy, bSz, 
							NeighborPairs, AIdxBlock, NIdxBlock, NIdxGridA, NIdxGridB, NIdxGridC, SIdx,
							Jij, Bij, Dij, VDMx, VDMy, VDMz, VKu1, Ku1, VKu2, Ku2, Kc, VHf, Hf, Etot, Mtot, NOS );
				// totalEnergy = 0;
				mtot[0] = Mtot[0]/NOS;
				mtot[1] = Mtot[1]/NOS;
				mtot[2] = Mtot[2]/NOS;
                // int  k=123;
                // mtot[2] = sqrt(bSx[k]*bSx[k]+bSy[k]*bSy[k]+bSz[k]*bSz[k]);//metka
				perSpEnergy = totalEnergy/NOS;
				totalEnergyFerro = GetTotalEnergyFerro( VHf[0], VHf[1], VHf[2], 
							NeighborPairs, AIdxBlock, NIdxBlock, NIdxGridA, NIdxGridB, NIdxGridC, SIdx,
							Jij, Bij, Dij, VDMx, VDMy, VDMz, VKu1, Ku1, VKu2, Ku2, Kc, VHf, Hf, Etot, NOS );
				// totalEnergyFerro = 0;
				totalEnergyFerro = totalEnergyFerro/NOS;	
				perSpEnergyMinusFerro = perSpEnergy - totalEnergyFerro;
				SpecialEvent=0;
			}
			if (timeInterval > 500){//~0.5 second
				FPS = frameCount / (timeInterval * 0.002f);
				previousTime = currentTime;
				frameCount = 0;	
				if (Play!=0){
					IPS = (currentIteration - previousIteration)/ (timeInterval * 0.002f);;
					previousIteration = currentIteration;
				}
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


void TW_CALL CB_SetN_Multisample(const void *value, void *clientData )
{
	// (void)clientData; // unused
    // N_Multisample = 16;//*( int *)value;
    // N_Multisample = *( int *)value;
    // #if !defined(__APPLE__)
	// glutSetOption(GLUT_MULTISAMPLE, N_Multisample);
	// glutInitDisplayMode( GLUT_RGB|GLUT_DOUBLE|GLUT_DEPTH|GLUT_MULTISAMPLE|GLUT_ALPHA);
	// #endif
	// ChangeVectorMode (0);
}


void TW_CALL CB_GetN_Multisample(void *value, void *clientData)
{
    *(int *)value = N_Multisample;
}


void TW_CALL CB_Set_Run( const void *value, void *clientData )
{
	Play = *( int *)value;
    if (Play!=0){
        pthread_mutex_lock(&culc_mutex);
            ENGINE_MUTEX=DO_IT;
            SleepTime=100;
        pthread_mutex_unlock(&culc_mutex);
    }else{
        pthread_mutex_lock(&culc_mutex);
            ENGINE_MUTEX=WAIT;
            SleepTime=3000;
        pthread_mutex_unlock(&culc_mutex);  
    }
    // printf("IchBinHier!\n");

}

void TW_CALL CB_Get_Run(void *value, void *clientData)
{
    *(float *)value = Play;
}


void TW_CALL CB_SetRotX(const void *value, void *clientData )
{
    GetEulerFromQuaternion(q_Rotation, Rot); 
    Rot[0] = *( float *)value;
    GetQuaternionFromEuler(q_Rotation, Rot);
}

void TW_CALL CB_GetRotX(void *value, void *clientData)
{
    *(float *)value = Rot[0];
}


void TW_CALL CB_SetRotY(const void *value, void *clientData )
{
    GetEulerFromQuaternion(q_Rotation, Rot); 
    Rot[1] = *( float *)value;
    GetQuaternionFromEuler(q_Rotation, Rot);
}

void TW_CALL CB_GetRotY(void *value, void *clientData)
{
    *(float *)value = Rot[1];
}

void TW_CALL CB_SetRotZ(const void *value, void *clientData )
{
    GetEulerFromQuaternion(q_Rotation, Rot); 
    Rot[2] = *( float *)value;
    GetQuaternionFromEuler(q_Rotation, Rot);
}

void TW_CALL CB_GetRotZ(void *value, void *clientData)
{
    *(float *)value = Rot[2];
}

void TW_CALL CB_SetScale(const void *value, void *clientData )
{
    Scale = *( float *)value; // copy value to Scale
	ChangeVectorMode (0);
}

//  Callback function called by the tweak bar to get the 'AutoRotate' value
void TW_CALL CB_GetScale(void *value, void *clientData)
{
    *(float *)value = Scale; // just copy Scale to value
}


void TW_CALL CB_SetScaleH(const void *value, void *clientData )
{
    Scale_H = *( float *)value; // copy value to Scale
    ChangeVectorMode (0);
}

//  Callback function called by the tweak bar to get the 'AutoRotate' value
void TW_CALL CB_GetScaleH(void *value, void *clientData)
{
    *(float *)value = Scale_H; // just copy Scale to value
}



void TW_CALL CB_SetVectorMode(const void *value, void *clientData )
{
    WhichVectorMode = *( enVectorMode *)value; // copy value to WhichVectorMode
    ChangeVectorMode (0);
}


void TW_CALL CB_GetVectorMode(void *value, void *clientData)
{
    *(int *)value = WhichVectorMode; // just copy WhichVectorMode to value
}


void TW_CALL CB_SetFaces(const void *value, void *clientData )
{
    arrowFaces = *( int *)value; // copy value to arrowFaces
	ChangeVectorMode (0);
}


void TW_CALL CB_GetFaces(void *value, void *clientData)
{
    *(int *)value = arrowFaces; // just copy arrowFaces to value
}

void TW_CALL CB_SetPivot(const void *value, void *clientData )
{
    Pivot = *( float *)value; // copy value to Pivot
	ChangeVectorMode (0);
}

void TW_CALL CB_GetPivot(void *value, void *clientData)
{
    *(float *)value = Pivot; // just copy Pivot to value
}

void TW_CALL CB_SaveCameraPosition ( void *clientData )//metka has to be fixed beacuse of new quaternion rotation 
{
	CameraPosition[CurrentCameraPositionBank][0]=Rot[0];
	CameraPosition[CurrentCameraPositionBank][1]=Rot[1];
	CameraPosition[CurrentCameraPositionBank][2]=Rot[2];
	CameraPosition[CurrentCameraPositionBank][3]=TransXYZ[0];
	CameraPosition[CurrentCameraPositionBank][4]=TransXYZ[1];
	CameraPosition[CurrentCameraPositionBank][5]=TransXYZ[2];
	CameraPosition[CurrentCameraPositionBank][6]=PerspSet[0];
}

void TW_CALL CB_GetCameraPosition ( void *clientData )//metka has to be fixed beacuse of new quaternion rotation
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
    ColorShift = *( int *)value; // copy value to ColorShift
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
		ChangeVectorMode (1);
}

void TW_CALL CB_GetColorShift(void *value, void *clientData)
{
    *(int *)value = ColorShift; // just copy ColorShift to value
}

void TW_CALL CB_SetInvHue(const void *value, void *clientData )
{
    InvertHue = *( int *)value; // copy value to InvertHue
	ChangeVectorMode (1);
}

void TW_CALL CB_GetInvHue(void *value, void *clientData)
{
    *(int *)value = InvertHue; // just copy InvertHue to value
}

void TW_CALL CB_SetInvVal(const void *value, void *clientData )
{
    InvertValue = *( int *)value; // copy value to InvertValue
	ChangeVectorMode (1);
}

void TW_CALL CB_GetInvVal(void *value, void *clientData)
{
    *(int *)value = InvertValue; // just copy InvertValue to value
}


void TW_CALL CB_SetHfield(const void *value, void *clientData )
{
    Hf = *( float *)value; // copy value to InvertValue
    // metka
    // if (Hf>1.1*fabs(Jij[0])) Hf=1.1*fabs(Jij[0]);
    Bdc[0] = Hf * VHf[0];
    Bdc[1] = Hf * VHf[1];
    Bdc[2] = Hf * VHf[2];
	UpdateVerticesNormalsColors_H(vertexProto_H, normalProto_H, VCNum_H, vertices_H, normals_H, colors_H, Box[0][0]*0.6, Box[1][1]*0.6, Box[2][2]*0.6, VHf[0], VHf[1], VHf[2]);
	UpdateVBO_H(&vboIdV_H, &vboIdN_H, &vboIdC_H, &iboIdI_H, vertices_H, normals_H, colors_H, indices_H);
}

void TW_CALL CB_GetHfield(void *value, void *clientData)
{
    *(float *)value = Hf; // just copy InvertValue to value
}


void TW_CALL CB_SetHfieldTheta(const void *value, void *clientData )
{
    VHtheta = *(float*)value;
    VHf[0]=sin(PI*VHtheta/180)*cos(PI*VHphi/180);
	VHf[1]=sin(PI*VHtheta/180)*sin(PI*VHphi/180);
	VHf[2]=cos(PI*VHtheta/180);
	Bdc[0]=Hf*VHf[0];
	Bdc[1]=Hf*VHf[1];
	Bdc[2]=Hf*VHf[2];
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
	Bdc[0]=Hf*VHf[0];
	Bdc[1]=Hf*VHf[1];
	Bdc[2]=Hf*VHf[2];
	UpdateVerticesNormalsColors_H(vertexProto_H, normalProto_H, VCNum_H, vertices_H, normals_H, colors_H, Box[0][0]*0.6, Box[1][1]*0.6, Box[2][2]*0.6, VHf[0], VHf[1], VHf[2]);
	UpdateVBO_H(&vboIdV_H, &vboIdN_H, &vboIdC_H, &iboIdI_H, vertices_H, normals_H, colors_H, indices_H);
}

void TW_CALL CB_GetHfieldPhi(void *value, void *clientData)
{
    (void)clientData; 
    *(float*)value = VHphi; 
}



// void TW_CALL CB_SetHfieldXYZ(const void *value, void *clientData )
// {
// 	(void)clientData; // unused
//  //    VHtheta = *(float*)value;
//  //    VHf[0]=sin(PI*VHtheta/180)*cos(PI*VHphi/180);
// 	// VHf[1]=sin(PI*VHtheta/180)*sin(PI*VHphi/180);
// 	// VHf[2]=cos(PI*VHtheta/180);
// 	UpdateVerticesNormalsColors_H(vertexProto_H, normalProto_H, VCNum_H, vertices_H, normals_H, colors_H, Box[0][0]*0.8, Box[1][1]*0.8, Box[2][2]*0.8, VHf[0], VHf[1], VHf[2]);
// 	UpdateVBO_H(&vboIdV_H, &vboIdN_H, &vboIdC_H, &iboIdI_H, vertices_H, normals_H, colors_H, indices_H);
// }

// void TW_CALL CB_GetHfieldXYZ(void *value, void *clientData)
// {
//     (void)clientData; 
//     // *(float*)value = VHtheta; //metka check this comment!
// }


void TW_CALL CB_SetNumImages(const void *value, void *clientData )
{
    (void)clientData; // unused
    num_images = *( int *)value; // copy value to Period_dc
    ReallocateMemoryForImages(num_images, NOS);
}

void TW_CALL CB_GetNumImages(void *value, void *clientData)
{
    (void)clientData; // unused
    *(int *)value = num_images; // just copy Period_dc to value
}

void TW_CALL CB_SetACPeriod(const void *value, void *clientData )
{
    (void)clientData; // unused
    Period_dc = *( double *)value; // copy value to Period_dc
    Omega_dc = TPI/Period_dc;
}

void TW_CALL CB_GetACPeriod(void *value, void *clientData)
{
    (void)clientData; // unused
    *(double *)value = Period_dc; // just copy Period_dc to value
}

void TW_CALL CB_SetOmega(const void *value, void *clientData )
{
    (void)clientData; // unused
    Omega_dc = *( double *)value; // copy value to Period_dc
    Period_dc = TPI/Omega_dc;
}

void TW_CALL CB_GetOmega(void *value, void *clientData)
{
    (void)clientData; // unused
    *(double *)value = Omega_dc; // just copy Omega_dc to value
}

void TW_CALL CB_SetInitial( void *clientData )
{
	ChangeInitialState( WhichInitialState );
}

void
UpdateKind(int* Kind,float* Px, float* Py, float* Pz, int NOS, int NOSK)
{
	float dist, dist_max = chSizeG * chSizeG;
	switch(WhichGeometry){
		case CILINDER_G:
			for (int i=0; i<NOS; i++){
				dist = Px[i]*Px[i]+Py[i]*Py[i];
				if (dist>dist_max){
					Kind[i] = 0;
				}else{
					Kind[i] = 1;
					NOSK++;
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
					NOSK++;
				}
			}
		break;

		default:
			for (int i=0; i<NOS; i++){
					Kind[i] = 1;
				}
			NOSK = NOS;
		break;
	}

    for (int n=0; n<NOS; n++)
    {   
        Sx[n]*= Kind[n];
        Sy[n]*= Kind[n];
        Sz[n]*= Kind[n];

        bSx[n]*= Kind[n];
        bSy[n]*= Kind[n];
        bSz[n]*= Kind[n]; 
        
        tSx[n]*= Kind[n];
        tSy[n]*= Kind[n];
        tSz[n]*= Kind[n]; 
        
        t2Sx[n]*= Kind[n];
        t2Sy[n]*= Kind[n];
        t2Sz[n]*= Kind[n];  

        t3Sx[n]*= Kind[n];
        t3Sy[n]*= Kind[n];
        t3Sz[n]*= Kind[n]; 
    }

}

void TW_CALL CB_SetShape( void *clientData )
{
	UpdateKind(Kind, Px, Py, Pz, NOS, NOSK);
	ChangeVectorMode(1);
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
    recordsCounter=0;
 
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

void TW_CALL CB_GetThetaMax1(void *value, void *clientData){
    (void)clientData; // unused
    *(int *)value = theta_max1; 
}

void TW_CALL CB_SetThetaMax1(const void *value, void *clientData ){
	(void)clientData; // unused
	int test= *( int *)value; 
	if (test>=theta_min1 ){
        theta_max1 = test;
        Sz_min1=cos(theta_max1*PI/180.0);
        ChangeVectorMode(0);
	}
}

void TW_CALL CB_GetThetaMax2(void *value, void *clientData){
    (void)clientData; // unused
    *(int *)value = theta_max2; 
}

void TW_CALL CB_SetThetaMax2(const void *value, void *clientData ){
	(void)clientData; // unused
	int test= *( int *)value; 
	if (test>=theta_min2 ){
        theta_max2 = test;
        Sz_min2=cos(theta_max2*PI/180.0);
        ChangeVectorMode(0);
	}
}

void TW_CALL CB_GetThetaMax3(void *value, void *clientData){
    (void)clientData; // unused
    *(int *)value = theta_max3; 
}

void TW_CALL CB_SetThetaMax3(const void *value, void *clientData ){
	(void)clientData; // unused
	int test= *( int *)value; 
	if (test>=theta_min3 ){
        theta_max3 = test;
        Sz_min3=cos(theta_max3*PI/180.0);
        ChangeVectorMode(0);
	}
}

void TW_CALL CB_GetThetaMin1(void *value, void *clientData){
    (void)clientData; // unused
    *(int *)value = theta_min1; 
}

void TW_CALL CB_SetThetaMin1(const void *value, void *clientData ){
	(void)clientData; // unused
	int test= *( int *)value; 
	if (test<=theta_max1 ){
        theta_min1 = test; 
        Sz_max1=cos(theta_min1*PI/180.0);
        ChangeVectorMode(0);
	}
}

void TW_CALL CB_GetThetaMin2(void *value, void *clientData){
    (void)clientData; // unused
    *(int *)value = theta_min2; 
}

void TW_CALL CB_SetThetaMin2(const void *value, void *clientData ){
	(void)clientData; // unused
	int test= *( int *)value; 
	if (test<=theta_max2 ){
        theta_min2 = test; 
        Sz_max2=cos(theta_min2*PI/180.0);
        ChangeVectorMode(0);
	}
}

void TW_CALL CB_GetThetaMin3(void *value, void *clientData){
    (void)clientData; // unused
    *(int *)value = theta_min3; 
}

void TW_CALL CB_SetThetaMin3(const void *value, void *clientData ){
	(void)clientData; // unused
	int test= *( int *)value; 
	if (test<=theta_max3 ){
        theta_min3 = test; 
        Sz_max3=cos(theta_min3*PI/180.0);
        ChangeVectorMode(0);
	}
}

void TW_CALL CB_GetPhiMax1(void *value, void *clientData){
    (void)clientData; // unused
    *(int *)value = phi_max1; 
}


void TW_CALL CB_SetPhiMax1(const void *value, void *clientData ){
	(void)clientData; // unused
	int test= *( int *)value; 
	if (test>=phi_min1 ){
        phi_max1 = test;
        ChangeVectorMode(0);
	}
}

void TW_CALL CB_GetPhiMin1(void *value, void *clientData){
    (void)clientData; // unused
    *(int *)value = phi_min1; 
}


void TW_CALL CB_SetPhiMin1(const void *value, void *clientData ){
	(void)clientData; // unused
	int test= *( int *)value; 
	if (test<=phi_max1 ){
        phi_min1 = test;
        ChangeVectorMode(0);
	}
}

void TW_CALL CB_GetPhiMax2(void *value, void *clientData){
    (void)clientData; // unused
    *(int *)value = phi_max2; 
}


void TW_CALL CB_SetPhiMax2(const void *value, void *clientData ){
	(void)clientData; // unused
	int test= *( int *)value; 
	if (test>=phi_min2 ){
        phi_max2 = test;
        ChangeVectorMode(0);
	}
}

void TW_CALL CB_GetPhiMin2(void *value, void *clientData){
    (void)clientData; // unused
    *(int *)value = phi_min2; 
}


void TW_CALL CB_SetPhiMin2(const void *value, void *clientData ){
	(void)clientData; // unused
	int test= *( int *)value; 
	if (test<=phi_max2 ){
        phi_min2 = test;
        ChangeVectorMode(0);
	}
}

void TW_CALL CB_GetPhiMax3(void *value, void *clientData){
    (void)clientData; // unused
    *(int *)value = phi_max3; 
}


void TW_CALL CB_SetPhiMax3(const void *value, void *clientData ){
	(void)clientData; // unused
	int test= *( int *)value; 
	if (test>=phi_min3 ){
        phi_max3 = test;
        ChangeVectorMode(0);
	}
}

void TW_CALL CB_GetPhiMin3(void *value, void *clientData){
    (void)clientData; // unused
    *(int *)value = phi_min3; 
}

void TW_CALL CB_SetPhiMin3(const void *value, void *clientData ){
	(void)clientData; // unused
	int test= *( int *)value; 
	if (test<=phi_max3 ){
        phi_min3 = test;
        ChangeVectorMode(0);
	}
}
/**********************************************************************/
void TW_CALL CB_GetGreedMaxA(void *value, void *clientData){
    *(int *)value = GreedFilterMaxA; 
}

void TW_CALL CB_SetGreedMaxA(const void *value, void *clientData ){
	int test= *( int *)value; 
	if (test<uABC[0] && test>GreedFilterMinA){
        GreedFilterMaxA = test;
        ChangeVectorMode(0);
	}
}

void TW_CALL CB_GetGreedMinA(void *value, void *clientData){
    *(int *)value = GreedFilterMinA; 
}

void TW_CALL CB_SetGreedMinA(const void *value, void *clientData ){
	int test= *( int *)value; 
	if (test<GreedFilterMaxA && test>=0){
        GreedFilterMinA = test;
        ChangeVectorMode(0);
	}
}
/*************/
void TW_CALL CB_GetGreedMaxB(void *value, void *clientData){
    *(int *)value = GreedFilterMaxB; 
}

void TW_CALL CB_SetGreedMaxB(const void *value, void *clientData ){
	int test= *( int *)value; 
	if (test<uABC[1] && test>GreedFilterMinB){
        GreedFilterMaxB = test;
        ChangeVectorMode(0);
	}
}

void TW_CALL CB_GetGreedMinB(void *value, void *clientData){
    *(int *)value = GreedFilterMinB; 
}

void TW_CALL CB_SetGreedMinB(const void *value, void *clientData ){
	int test= *( int *)value; 
	if (test<GreedFilterMaxB && test>=0){
        GreedFilterMinB = test;
        ChangeVectorMode(0);
	}
}
/*************/
void TW_CALL CB_GetGreedMaxC(void *value, void *clientData){
    *(int *)value = GreedFilterMaxC; 
}

void TW_CALL CB_SetGreedMaxC(const void *value, void *clientData ){
	int test= *( int *)value; 
	if (test<uABC[2] && test>GreedFilterMinC){
        GreedFilterMaxC = test;
        ChangeVectorMode(0);
	}
}

void TW_CALL CB_GetGreedMinC(void *value, void *clientData){
    *(int *)value = GreedFilterMinC; 
}

void TW_CALL CB_SetGreedMinC(const void *value, void *clientData ){
	int test= *( int *)value; 
	if (test<GreedFilterMaxC && test>=0){
        GreedFilterMinC = test;
        ChangeVectorMode(0);
	}
}

void TW_CALL CB_GetGreedFilterInvert(void *value, void *clientData){
    *(bool *)value = GreedFilterInvert; 
}

void TW_CALL CB_SetGreedFilterInvert(const void *value, void *clientData ){
	GreedFilterInvert= *( bool *)value; 
    ChangeVectorMode(0);
}

void TW_CALL CB_GetSpinFilter1(void *value, void *clientData){
    *(bool *)value = SpinFilter1; 
}

void TW_CALL CB_SetSpinFilter1(const void *value, void *clientData ){
	SpinFilter1= *( bool *)value; 
    ChangeVectorMode(0);
}

void TW_CALL CB_GetSpinFilter2(void *value, void *clientData){
    *(bool *)value = SpinFilter2; 
}

void TW_CALL CB_SetSpinFilter2(const void *value, void *clientData ){
	SpinFilter2= *( bool *)value; 
    ChangeVectorMode(0);
}

void TW_CALL CB_GetSpinFilter3(void *value, void *clientData){
    *(bool *)value = SpinFilter3; 
}

void TW_CALL CB_SetSpinFilter3(const void *value, void *clientData ){
	SpinFilter3= *( bool *)value; 
    ChangeVectorMode(0);
}

void TW_CALL CB_GetGreedFilter(void *value, void *clientData){
    *(bool *)value = GreedFilter; 
}

void TW_CALL CB_SetGreedFilter(const void *value, void *clientData ){
	GreedFilter= *( bool *)value; 
    ChangeVectorMode(0);
}


void TW_CALL CB_ResetIterations( void *clientData ){
  ITERATION=0;
  currentIteration=0;
  SpecialEvent=1;
}



void TW_CALL CB_SaveCSV( void *clientData )
{
    FILE * pFile;
    pFile = fopen (outputfilename,"w");
   
    if (pFile!=NULL)
    { 	
        // fputs ("px,py,pz,nx,ny,nz,\n",pFile);
        // for (int i=0;i<NOS;i++)
        // {
        //     snprintf(shortBufer,200,"%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,\n",Px[i],Py[i],Pz[i],bSx[i],bSy[i],bSz[i]);
        //     fputs (shortBufer,pFile);
        // }
        // fclose (pFile);

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

            // for export from x,y,z--> x,z,Nc-y
            // from Filipp's code
            // for (bn = bnini; bn<bnfin; bn++) {
            //     for (cn = cnfin-1; cn>=cnini; cn--) {
            //         for (an = anini; an<anfin; an++) {
            //             n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
            //             n = n*AtomsPerBlock;//index of the first spin in the block
            //             for (atom=0; atom<AtomsPerBlock; atom++){
            //                 N = n + atom;
            //                 snprintf(shortBufer,200,"%.14g,%.14g,%.14g,\n", Sx[N],-Sz[N],Sy[N]);
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
    					n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
    					n = n*AtomsPerBlock;//index of the first spin in the block
    					for (atom=0; atom<AtomsPerBlock; atom++){
    					    N = n + atom;
                            //metka
    						snprintf(shortBufer,200,"%2.5f,%2.5f,%2.5f,%0.15f,%0.15f,%0.15f,\n",Px[N],Py[N],Pz[N],bSx[N],bSy[N],bSz[N]);
                            // snprintf(shortBufer,200,"%.14g,%.14g,%.14g,%.14g,%.14g,%.14g,%.14g,%.14g,\n",
                                // Px[N],Py[N],Pz[N],bSx[N],bSy[N],bSz[N],acos(bSz[N])*R2D,atan2 (bSy[N],bSx[N])*R2D);
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
    					snprintf(shortBufer,200,"%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,\n",tPositin[0],tPositin[1],tPositin[2],tSpin[0]/cN,tSpin[1]/cN,tSpin[2]/cN,modulus/cN);
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
    					snprintf(shortBufer,200,"%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,\n",tPositin[0],tPositin[1],tPositin[2],tSpin[0]/bN,tSpin[1]/bN,tSpin[2]/bN,modulus/bN);
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
    					snprintf(shortBufer,200,"%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,\n",tPositin[0],tPositin[1],tPositin[2],tSpin[0]/aN,tSpin[1]/aN,tSpin[2]/aN,modulus/aN);
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
    //metka dlya schiutyvaniya equilibrium state for dm
    for (int i=0; i<NOS; i++){
        t3Sx[i]=Sx[i];
        t3Sy[i]=Sy[i];
        t3Sz[i]=Sz[i];                       
    }
	ChangeVectorMode(1);
}

void TW_CALL CB_ReadBIN( void *clientData )
{
	

	unsigned short int num = 65535;



  	struct tfshortint {
   		unsigned short int t;
   		unsigned short int f;
  	};
  	int Nx = uABC[0], Ny = uABC[1], Nz = uABC[2];
  	FILE * FilePointer = fopen(inputfilename, "rb");
  	for(int k = 0; k<Nz; k++){
	    for(int j = 0; j<Ny; j++){
			for(int i = 0; i <Nx;i++){
				struct tfshortint my_par_red;
				if (fread(&my_par_red, sizeof(struct tfshortint), 1, FilePointer)){
					double nx,ny,nz;
					unsigned short int p=my_par_red.t, q=my_par_red.f;

					double T = (double)(p+0.5)*PI/num;
					double F = (double)2*(q+0.5)*PI/num;

					nx = sin(T)*cos(F);
					ny = sin(T)*sin(F);
					nz = cos(T);
					int n = (i)+(j)*Nx+k*Nx*Ny;
					Sx[n] = nx; Sy[n] = ny; Sz[n] = nz;
					bSx[n] = nx; bSy[n] = ny; bSz[n] = nz;					
				}
			}
	    }
	}
	fclose (FilePointer);
	ChangeVectorMode(1);
}


void TW_CALL CB_ReadVTK( void *clientData )
{
    Read_VTK(Sx, Sy, Sz, inputfilename);  
    ChangeVectorMode(1);
}


void TW_CALL CB_Save_OVF_b8( void *clientData )
{
    char ovf_filename[64] = "";
    strncpy(ovf_filename, outputfilename, strcspn (outputfilename, "."));
    strcat(ovf_filename, ".ovf");
    if(strncmp(ovf_filename, ".ovf",4)==0){
        printf("Enter the file name. It cannot be empty!");
    }else{
        Save_OVF_b8(bSx, bSy, bSz, ovf_filename);
    }
}

void TW_CALL CB_Save_VTK_b4( void *clientData )
{
    char vtk_filename[64] = "";
    strncpy(vtk_filename, outputfilename, strcspn (outputfilename, "."));
    strcat(vtk_filename, ".vtk");
    if(strncmp(vtk_filename, ".vtk",4)==0){
        printf("Enter the file name. It cannot be empty!");
    }else{
        Save_VTK(bSx, bSy, bSz, 0, vtk_filename);//metka 0->1

        // Save_VTS_b4(bSx, bSy, bSz, Px, Py, Pz, Box, vts_filename);
        // Save_VTS_ascii(bSx, bSy, bSz, Px, Py, Pz, Box, vts_filename);

        // float * Spins_xyz;
        // Spins_xyz = (float *)calloc(NOS*3, sizeof(float));   
        // for (int n=0; n<NOS; n++){
        //         Spins_xyz[3*n+0]=Sx[n];
        //         Spins_xyz[3*n+1]=Sy[n];
        //         Spins_xyz[3*n+2]=Sz[n];         
        // }
        //     save_vtk(vts_filename,"name",3,Spins_xyz,"special_flag",1,Kind,uABC[0],uABC[1],uABC[2],uABC[0],uABC[1],uABC[2],1);
        // 
    }
}

void TW_CALL CB_Save_BIN( void *clientData )
{
    char bin_filename[64] = "";
    strncpy(bin_filename, outputfilename, strcspn (outputfilename, "."));
    strcat(bin_filename, ".bin");
    if(strncmp(bin_filename, ".bin",4)==0){
        printf("Enter the file name. It cannot be empty!");
    }else{
        SaveBin(bSx, bSy, bSz, bin_filename);//metka 0->1
    }
}

void TW_CALL CB_Save_BMP( void *clientData )
{
    char bmp_filename[64] = "";
    strncpy(bmp_filename, outputfilename, strcspn (outputfilename, "."));
    strcat(bmp_filename, ".bmp");
    if(strncmp(bmp_filename, ".bin",4)==0){
        printf("Enter the file name. It cannot be empty!");
    }else{
        SaveBmp(bSx, bSy, bSz, bmp_filename, WhichSliceMode, A_layer_min-1, B_layer_min-1, C_layer_min-1);//metka 0->1
    }
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

			GreedFilterMaxA=uABC[0]-1;
			GreedFilterMaxB=uABC[1]-1;
			GreedFilterMaxC=uABC[2]-1;
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

	// TwAddVarCB(view_bar, "Multisampling", TW_TYPE_INT32, CB_SetN_Multisample, CB_GetN_Multisample, &N_Multisample, " label='Multisamples' min=1 max=32 step=1 help='Multisampling' group='Camera'");
	{
	TwEnumVal		enProjectionsTw[] = { {ORTHO, "Orthogonal"}, {PERSP, "Perspective"} };
	TwType			TW_TYPE_PROJ = TwDefineEnum("ProjectionType", enProjectionsTw, 2);
	TwAddVarRW(view_bar, "Projection", TW_TYPE_PROJ, &WhichProjection, "keyIncr='p' help='Type of 3D projection' group='Camera'");
	}

    TwAddVarRW(view_bar, "ObjRotation", TW_TYPE_QUAT4F, &q_Rotation, " label='Scene rotation' opened=true help='Change the 3D scene orientation.' ");

	TwAddVarRW(view_bar, "CamAng", TW_TYPE_FLOAT, &PerspSet[0], " label='camera angle' min=1 max=120 help='camera angle' group='Camera'");
	TwAddVarRW(view_bar, "PosX", TW_TYPE_FLOAT, &TransXYZ[0], " label='position in x' min=-1000 max=1000 help='camera position along X-axis' group='Camera'");
	TwAddVarRW(view_bar, "PosY", TW_TYPE_FLOAT, &TransXYZ[1], " label='position in y' min=-1000 max=1000 help='camera position along Y-axis' group='Camera'");
	TwAddVarRW(view_bar, "PosZ", TW_TYPE_FLOAT, &TransXYZ[2], " label='position in z' min=-1000 max=1000 help='camera position along Z-axis' group='Camera'");

    TwAddVarCB(view_bar, "RotX", TW_TYPE_FLOAT, CB_SetRotX, CB_GetRotX, &Rot[0], " label='turn around X' help='rotate camera around X-axis' group='Camera'");

    TwAddVarCB(view_bar, "RotY", TW_TYPE_FLOAT, CB_SetRotY, CB_GetRotY, &Rot[1], " label='turn around Y' help='rotate camera around Y-axis' group='Camera'");

    TwAddVarCB(view_bar, "RotZ", TW_TYPE_FLOAT, CB_SetRotZ, CB_GetRotZ, &Rot[2], " label='turn around Z' help='rotate camera around Z-axis' group='Camera'");

	TwAddVarRW(view_bar, "RotSpeed", TW_TYPE_FLOAT, &RotSpeed, " label='rotation speed' min=0 max=10 step=1 help='speed of rotation around any axis' group='Camera'");
	TwAddVarRW(view_bar, "TransSpeed", TW_TYPE_FLOAT, &TransSpeed, " label='translation speed' min=0 max=10 step=1 help='speed of translation along any axis' group='Camera'");

	TwDefine(" View/Camera opened=false ");

	// TwAddVarRW(view_bar, "CamBank", TW_TYPE_INT32, &CurrentCameraPositionBank, " label='Current camera' min=0 max=4 group='CameraRW'");
	// TwAddButton(view_bar, "Read Camera", CB_GetCameraPosition, NULL, "label='read camera pos.' group='CameraRW'");
	// TwAddButton(view_bar, "Write Camera", CB_SaveCameraPosition, NULL, "label='save camera pos.' group='CameraRW'");
	// TwDefine(" View/CameraRW opened=false ");

	/****** Light ******/
	TwAddVarRW(view_bar, "Light_On_Off", TW_TYPE_BOOL32, &Light_On, " label='Light On/Off' key=l help='Reflectins' group='Light'");
	TwAddVarRW(view_bar, "Intensity", TW_TYPE_FLOAT, &g_LightMultiplier, " label='Light intensity' min=0.1 max=4 step=0.02 help='Increase/decrease the light power.' group='Light' ");
	TwAddVarRW(view_bar, "LightDir", TW_TYPE_DIR3F, &g_LightDirection, " label='Light direction' opened=false help='Change the light direction.' group='Light'");
	temp_color[0] = 230;
	temp_color[1] = 230;
	temp_color[2] = 255;
	TwSetParam(view_bar, "LightDir", "arrowcolor", TW_PARAM_INT32, 3, temp_color);


    TwDefine(" View/Light opened=false ");

	{
	TwEnumVal		enColorsTw[] = { {WHITE,"White"}, {BLACK, "Black"}, {RED, "Red"}, {GREEN, "Green"}, {BLUE, "Blue"}, {MANUAL, "Manual"} };
	TwType			TW_TYPE_COLOR = TwDefineEnum("BG_Color", enColorsTw, 6);
	TwAddVarRW(view_bar, "Choose_background", TW_TYPE_COLOR, &WhichBackgroundColor, "help='Background color for 3D scene' group='Background'");
	}
	TwAddVarRW(view_bar, "red", TW_TYPE_FLOAT, &BackgroundColors[5][0], 
	" min=0 max=1 step=0.01 group='Background'");
	TwAddVarRW(view_bar, "green", TW_TYPE_FLOAT, &BackgroundColors[5][1], 
	" min=0 max=1 step=0.01 group='Background'");
	TwAddVarRW(view_bar, "blue", TW_TYPE_FLOAT, &BackgroundColors[5][2], 
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
	TwAddVarCB(view_bar, "Type of vectors", TW_TYPE_VEC_MOD, CB_SetVectorMode, CB_GetVectorMode, &WhichVectorMode, "keyIncr='v' keyDecr='V' help='Type of 3D vectors' group='Appearance' ");
	}

	TwAddVarCB(view_bar, "Pivot", TW_TYPE_FLOAT, CB_SetPivot, CB_GetPivot,  &Pivot, " min=0 max=1 step=0.01 help='Pivot of 3D arrow.' group='Appearance' ");
	TwAddVarCB(view_bar, "Faces", TW_TYPE_INT32, CB_SetFaces, CB_GetFaces,  &arrowFaces, " min=3 max=20 step=1 help='Number of faces for 3D arrow.' group='Appearance' ");
	TwAddVarCB(view_bar, "Scale", TW_TYPE_FLOAT, CB_SetScale, CB_GetScale,  &Scale, " min=0.1 max=10 step=0.01 keyIncr='+' keyDecr='-' help='Scale the vectors.' group='Appearance' ");
    TwAddVarCB(view_bar, "Scale Bext", TW_TYPE_FLOAT, CB_SetScaleH, CB_GetScaleH,  &Scale_H, " min=0.1 max=10 step=0.01 keyIncr='+' keyDecr='-' help='Scale the vectors.' group='Appearance' ");


	TwAddVarRW(view_bar, "Show basis", TW_TYPE_BOOL32, &AxesOn, " key=CTRL+o  group='Appearance'");
	TwAddVarRW(view_bar, "Show box", TW_TYPE_BOOL32, &BoxOn, " key=CTRL+b group='Appearance'");
	TwDefine(" View/Appearance opened=false ");

	{
	TwEnumVal		enSliceModeTw[] = { 
										{A_AXIS, "a-axis"}, 
										{B_AXIS, "b-axis"}, 
										{C_AXIS, "c-axis"},
										{FILTER, "filter"}
									  };
	TwType			TV_TYPE_VEC_MOD = TwDefineEnum("Slicing", enSliceModeTw, 4);
	TwAddVarCB(view_bar, "Slicing mode", TV_TYPE_VEC_MOD, CB_SetSliceMode, CB_GetSliceMode, &WhichSliceMode, "keyIncr='/' keyDecr='?' help='Slising plane perpenticulat to the choosen axis' group='Filters&Slices' ");
	}
	TwAddVarCB(view_bar, "GreedFilter", TW_TYPE_BOOL32, CB_SetGreedFilter, CB_GetGreedFilter, &GreedFilter, 
		" label='Greed filter' true='On' false='Off' group='Filters&Slices' ");
	TwAddVarCB(view_bar, "Spin_filter1", TW_TYPE_BOOL32, CB_SetSpinFilter1, CB_GetSpinFilter1, &SpinFilter1, 
		" label='Spin filter N1' true='On' false='Off' group='Filters&Slices' ");
	TwAddVarCB(view_bar, "Spin_filter2", TW_TYPE_BOOL32, CB_SetSpinFilter2, CB_GetSpinFilter2, &SpinFilter2, 
		" label='Spin filter N2' true='On' false='Off' group='Filters&Slices' ");
	TwAddVarCB(view_bar, "Spin_filter3", TW_TYPE_BOOL32, CB_SetSpinFilter3, CB_GetSpinFilter3, &SpinFilter3, 
		" label='Spin filter N3' true='On' false='Off' group='Filters&Slices' ");

	TwAddVarCB(view_bar, "T_max1", TW_TYPE_INT32, CB_SetThetaMax1, CB_GetThetaMax1, &theta_max1, 
		" group='Spin_filter_1' label='Theta max 1' min=0 max=180  help='max value for polar angle theta'");
	TwAddVarCB(view_bar, "T_min1", TW_TYPE_INT32, CB_SetThetaMin1, CB_GetThetaMin1, &theta_min1, 
		" group='Spin_filter_1' label='Theta min 1' min=0 max=180  help='min value for polar angle theta'");
	TwAddVarCB(view_bar, "F_max1", TW_TYPE_INT32, CB_SetPhiMax1, CB_GetPhiMax1, &phi_max1, 
		" group='Spin_filter_1' label='Phi max 1' min=0 max=360 step=1 help='max value for azimuthal angle phi'");
	TwAddVarCB(view_bar, "F_min1", TW_TYPE_INT32, CB_SetPhiMin1, CB_GetPhiMin1, &phi_min1, 
		" group='Spin_filter_1' label='Phi min 1' min=0 max=360 step=1 help='min value for azimuthal angle phi'");

	TwAddVarCB(view_bar, "T_max2", TW_TYPE_INT32, CB_SetThetaMax2, CB_GetThetaMax2, &theta_max2, 
		" group='Spin_filter_2' label='Theta max 2' min=0 max=180  help='max value for polar angle theta'");
	TwAddVarCB(view_bar, "T_min2", TW_TYPE_INT32, CB_SetThetaMin2, CB_GetThetaMin2, &theta_min2, 
		" group='Spin_filter_2' label='Theta min 2' min=0 max=180  help='min value for polar angle theta'");	
	TwAddVarCB(view_bar, "F_max2", TW_TYPE_INT32, CB_SetPhiMax2, CB_GetPhiMax2, &phi_max2, 
		" group='Spin_filter_2' label='Phi max 2' min=0 max=360 step=1 help='max value for azimuthal angle phi'");
	TwAddVarCB(view_bar, "F_min2", TW_TYPE_INT32, CB_SetPhiMin2, CB_GetPhiMin2, &phi_min2, 
		" group='Spin_filter_2' label='Phi min 2' min=0 max=360 step=1 help='min value for azimuthal angle phi'");

	TwAddVarCB(view_bar, "T_max3", TW_TYPE_INT32, CB_SetThetaMax3, CB_GetThetaMax3, &theta_max3, 
		" group='Spin_filter_3' label='Theta max 3' min=0 max=180  help='max value for polar angle theta'");
	TwAddVarCB(view_bar, "T_min3", TW_TYPE_INT32, CB_SetThetaMin3, CB_GetThetaMin3, &theta_min3, 
		" group='Spin_filter_3' label='Theta min 3' min=0 max=180  help='min value for polar angle theta'");	
	TwAddVarCB(view_bar, "F_max3", TW_TYPE_INT32, CB_SetPhiMax3, CB_GetPhiMax3, &phi_max3, 
		" group='Spin_filter_3' label='Phi max 3' min=0 max=360 step=1 help='max value for azimuthal angle phi'");
	TwAddVarCB(view_bar, "F_min3", TW_TYPE_INT32, CB_SetPhiMin3, CB_GetPhiMin3, &phi_min3, 
		" group='Spin_filter_3' label='Phi min 3' min=0 max=360 step=1 help='min value for azimuthal angle phi'");

	TwDefine(" View/Spin_filter_1 opened=false group='Filters&Slices'");
	TwDefine(" View/Spin_filter_2 opened=false group='Filters&Slices'");
	TwDefine(" View/Spin_filter_3 opened=false group='Filters&Slices'");

	TwAddVarCB(view_bar, "GFmaxA", TW_TYPE_INT32, CB_SetGreedMaxA, CB_GetGreedMaxA, &GreedFilterMaxA, " group='Greed_filter' label='max na' ");
	TwAddVarCB(view_bar, "GFminA", TW_TYPE_INT32, CB_SetGreedMinA, CB_GetGreedMinA, &GreedFilterMinA, " group='Greed_filter' label='min na' ");
	TwAddVarCB(view_bar, "GFmaxB", TW_TYPE_INT32, CB_SetGreedMaxB, CB_GetGreedMaxB, &GreedFilterMaxB, " group='Greed_filter' label='max nb' ");
	TwAddVarCB(view_bar, "GFminB", TW_TYPE_INT32, CB_SetGreedMinB, CB_GetGreedMinB, &GreedFilterMinB, " group='Greed_filter' label='min nb' ");
	TwAddVarCB(view_bar, "GFmaxC", TW_TYPE_INT32, CB_SetGreedMaxC, CB_GetGreedMaxC, &GreedFilterMaxC, " group='Greed_filter' label='max nc' ");
	TwAddVarCB(view_bar, "GFminC", TW_TYPE_INT32, CB_SetGreedMinC, CB_GetGreedMinC, &GreedFilterMinC, " group='Greed_filter' label='min nc' ");

    TwAddVarCB(view_bar, "GreedFilterInvert", TW_TYPE_BOOL32, CB_SetGreedFilterInvert, CB_GetGreedFilterInvert, &GreedFilterInvert, " label='Invert G filter' true='On' false='Off' group='Greed_filter' ");


	TwDefine(" View/Greed_filter opened=false group='Filters&Slices'");
	TwDefine(" View/Filters&Slices opened=false ");


	TwAddVarCB(view_bar, "ColorShift", TW_TYPE_INT32, CB_SetColorShift, CB_GetColorShift, &ColorShift, " label='Rotate hue' min=0 max=360 help='rotate color hue in xy-plane' group='HSV to RGB'");
	TwAddVarCB(view_bar, "InvHue", TW_TYPE_BOOL32, CB_SetInvHue, CB_GetInvHue, &InvertHue, " label='Invert hue' help='invert RGB to RBG color hue' group='HSV to RGB'");
	TwAddVarCB(view_bar, "InvVal", TW_TYPE_BOOL32, CB_SetInvVal, CB_GetInvVal, &InvertValue, " label='Invert value' help='invert black to white' group='HSV to RGB'");
	TwDefine(" View/'HSV to RGB' opened=false ");

/*  Hamiltonian parameters&controls F3 */
	control_bar = TwNewBar("Parameters&Controls");
    TwDefine(" Parameters&Controls iconified=true "); 
    TwDefine(" Parameters&Controls size='220 530' color='100 70 100' alpha=200 "); // change default tweak bar size and color
    TwDefine(" Parameters&Controls help='F3: show/hide Control bar' "); // change default tweak bar size and color

	//TwAddButton(control_bar, "Run", CB_Run, NULL, "key='r' label='RUN simulation' ");
	TwAddVarCB(control_bar, "Run", TW_TYPE_BOOL32, CB_Set_Run, CB_Get_Run, &Play, " keyIncr='r' true='On' false='Off' label='RUN simulation' ");
	TwAddVarRW(control_bar, "Record", TW_TYPE_BOOL32, &Record, "label='Recording' true='On' false='Off' help='Recording <sx>, <sy>, <sz> in each iteration'");
	TwAddVarRW(control_bar, "Rec_Iteration", TW_TYPE_INT32, &rec_iteration, "label='Every iteration' min=1 max=1000 step=1 ");
    TwAddVarRW(control_bar, "Max_Iteration", TW_TYPE_INT32, &Max_Numb_Iteration, "label='Max. Iteration' min=1 max=100000000 step=100 ");
	TwAddButton(control_bar, "Clean the record", CB_CleanSxSySzFile, NULL, "label= 'Clean the record' help='clean the output file with <sx>, <sy>, <sz>' ");
	TwAddButton(control_bar, "Reset iterations", CB_ResetIterations, NULL, "label='Reset iterations' ");

	TwAddSeparator(control_bar, "sep-3", NULL);

	TwAddVarRW(control_bar, "BCinA", TW_TYPE_BOOL32, &Boundary[0], "label='along a' group='Boundary conditions' true='periodic' false='open' help='set boundary conditions along translation vector 'a' '");
	TwAddVarRW(control_bar, "BCinB", TW_TYPE_BOOL32, &Boundary[1], "label='along b' group='Boundary conditions' true='periodic' false='open' help='set boundary conditions along translation vector 'b' '");
	TwAddVarRW(control_bar, "BCinC", TW_TYPE_BOOL32, &Boundary[2], "label='along c' group='Boundary conditions' true='periodic' false='open' help='set boundary conditions along translation vector 'c' '");
    
	{
    TwEnumVal       enIntegrationScheme[] = {{HEUN, "Heun(1) "},
                                             {SIB,  " SIB(2) "},
                                             {RK23, " RK(23) "},
                                             {RK45, " RK(45) "},
                                         	 {RELAX, " RELAX "}};
    TwType          TV_TYPE_INTEGRATION_SCHEME = TwDefineEnum("Solver", enIntegrationScheme, 5);
    TwAddVarRW(control_bar, "Integration scheme", TV_TYPE_INTEGRATION_SCHEME, &WhichIntegrationScheme, "group='LLG' help='Choose the integration scheme'");
    }
    //TwAddVarRW(control_bar, "Preces", TW_TYPE_BOOL32, &Precession, "label='precession' group='LLG' true='On' false='Off' help='On/Off precession'");
    TwAddVarRW(control_bar, "Damping", TW_TYPE_FLOAT, &damping, "label='Damping' min=0 max=100 step=0.000001 group='LLG' ");
	TwAddVarRW(control_bar, "Time_step", TW_TYPE_FLOAT, &t_step, "label='Time step' min=0 max=10.0 step=0.000001   group='LLG' ");
	TwAddVarRW(control_bar, "temperature", TW_TYPE_FLOAT, &Temperature, "label='k_b*T' min=0 max=100 step=0.000001 group='LLG' ");

	TwAddSeparator(control_bar, "sep-4", NULL);
	TwAddVarRW(control_bar, "Xi", TW_TYPE_FLOAT, &Xi, "label='Xi'  step=0.001 help='Non-adiabaticity parameter' ");
	TwAddVarRW(control_bar, "Curr_u", TW_TYPE_FLOAT, &Curr_u, "label='Curr_u'  step=0.001 help='Current in Zhang-Li torque' ");

    TwAddSeparator(control_bar, "sep-2", NULL);
	TwAddVarCB(control_bar, "FieldTheta", TW_TYPE_FLOAT, CB_SetHfieldTheta, CB_GetHfieldTheta, &VHtheta, "label='H theta'  step=0.0001 help='Change the direction of applied field' ");
	TwAddVarCB(control_bar, "FieldPhi", TW_TYPE_FLOAT, CB_SetHfieldPhi, CB_GetHfieldPhi, &VHphi, "label='H phi' step=0.0001  help='Change the direction of applied field' ");
	TwAddVarCB(control_bar, "Field", TW_TYPE_FLOAT, CB_SetHfield, CB_GetHfield, &Hf, "label='Field'  min=0 step=0.0000001 help='Change the strength of applied field' ");
	// TwAddVarCB(control_bar, "FieldDir", TW_TYPE_DIR3F, CB_SetHfieldDir, CB_GetHfieldDir, VHf, 
	// "label='Field direction' opened=true help='Change the direction of applied field' ");
	// temp_color[0] = 55;
	// temp_color[1] = 55;
	// temp_color[2] = 155;
	// TwSetParam(control_bar, "FieldDir", "arrowcolor", TW_PARAM_INT32, 3, temp_color);
	// // TwAddVarRW(control_bar, "Field", TW_TYPE_FLOAT, &Hf, 
	// // "label='Field' help='The value of uniaxial anisotropy' ");
	// TwAddSeparator(control_bar, "control_sep1", NULL);
	// TwAddVarCB(control_bar, "FieldDir", TW_TYPE_DIR3F, CB_SetHfieldXYZ, CB_GetHfieldXYZ, &VHf, 
	// " label='Field direction' opened=false help='Change the applied field direction.' ");

	TwAddSeparator(control_bar, "control_sep2", NULL);


	TwAddVarRW(control_bar, "Kud1Dir", TW_TYPE_DIR3F, &VKu1, 
	"label='Ku1 axis' opened=true help='The axis of the 1-st uniaxial anisotropy' ");
	temp_color[0] = 55;
	temp_color[1] = 155;
	temp_color[2] = 55;
	TwSetParam(control_bar, "Kud1Dir", "arrowcolor", TW_PARAM_INT32, 3, temp_color);
	TwAddVarRW(control_bar, "Ku", TW_TYPE_FLOAT, &Ku1, 
	"label='Ku1' help='The value of uniaxial anisotropy' ");
    //////////////////////////////////////////
    TwAddSeparator(control_bar, "sepbetweenKu1Ku2", NULL);
    //////////////////////////////////////////
    TwAddVarRW(control_bar, "Kud2Dir", TW_TYPE_DIR3F, &VKu2, 
    "label='Ku2 axis' opened=true help='The axis of the 2-nd uniaxial anisotropy' ");
    TwSetParam(control_bar, "Kud2Dir", "arrowcolor", TW_PARAM_INT32, 3, temp_color);
    TwAddVarRW(control_bar, "Ku1", TW_TYPE_FLOAT, &Ku2, 
    "label='Ku2' help='The value of the 2-nd uniaxial anisotropy' ");

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
	snprintf(shortBufer,200,"J%1i",s);
	TwAddVarRW(control_bar, shortBufer, TW_TYPE_FLOAT, &Jij[s], "help='Heisenberg exchange' ");
	TwSetParam(control_bar, shortBufer, "label",  TW_PARAM_CSTRING , 1, shortBufer);
	}
	//////////////////////////////////////////
	TwAddSeparator(control_bar, "sep2", NULL);
	//////////////////////////////////////////	
	for(int s=0; s<ShellNumber; s++)
	{
	snprintf(shortBufer,200,"B%1i",s);
	TwAddVarRW(control_bar, shortBufer, TW_TYPE_FLOAT, &Bij[s], "help='Bi-quadratic exchange' ");
	TwSetParam(control_bar, shortBufer, "label",  TW_PARAM_CSTRING , 1, shortBufer);
	}
	//////////////////////////////////////////
	TwAddSeparator(control_bar, "sep3", NULL);
	//////////////////////////////////////////
	for(int s=0; s<ShellNumber; s++)
	{
	snprintf(shortBufer,200,"D%1i",s);
	TwAddVarRW(control_bar, shortBufer, TW_TYPE_FLOAT, &Dij[s], "help='Dzyaloshinskii-Moriya' ");
	TwSetParam(control_bar, shortBufer, "label",  TW_PARAM_CSTRING , 1, shortBufer);
	}
	//////////////////////////////////////////
	TwAddSeparator(control_bar, "sep4", NULL);
	//////////////////////////////////////////
	TwAddVarRW(control_bar, "ST param.", TW_TYPE_FLOAT, &Cu, "help='Dzyaloshinskii-Moriya' ");
	TwAddVarRW(control_bar, "CurrentDir", TW_TYPE_DIR3F, &VCu, "label='cur. dir.' opened=true help='The polarization direction of electric current' ");

/*  Initial state F4 */
	initial_bar = TwNewBar("Initial_State");
	TwDefine(" Initial_State iconified=true "); 
	TwDefine(" Initial_State size='220 530' color='70 70 100'  alpha=200"); // change default tweak bar size and color
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

	TwAddButton(initial_bar, "Set shape", CB_SetShape, NULL, " label='Set shape' ");

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
	TwAddVarRW(initial_bar, "Choose ini. state", TV_TYPE_INI_STATE, &WhichInitialState, "help='Choose initial spin configuration'");
	}

	TwAddVarRW(initial_bar, "chSize", TW_TYPE_FLOAT,  &chSize, " min=-100000 max=100000 step=0.5 help='characteristic size of modulated state: Skyrmion diameter or spiral period' ");

	TwAddVarRW(initial_bar, "chDir", TW_TYPE_DIR3F,  &chDir, " help='direction of modulations e.g. k-vector of the spin spiral.' ");

	TwAddButton(initial_bar, "Set initial", CB_SetInitial, NULL, "key=I label='insert initial state' ");

	TwAddSeparator(initial_bar, "sep01", NULL);
	TwAddVarRW(initial_bar, "Degrees", TW_TYPE_FLOAT,  &RotateAllSpins, " min=-360 max=360 step=1 help='Rotate all spins about characteristic direction' ");
	TwAddButton(initial_bar, "Rotate spins", CB_RotateAllSpins, NULL, "label='rotate all spins' ");

	TwAddSeparator(initial_bar, "sep1", NULL);
	TwAddButton(initial_bar, "Invert X", CB_InvertX, NULL, "label='invert n_x component' ");
	TwAddButton(initial_bar, "Invert Y", CB_InvertY, NULL, "label='invert n_y component' ");
	TwAddButton(initial_bar, "Invert Z", CB_InvertZ, NULL, "label='invert n_z component' ");

	TwAddSeparator(initial_bar, "sep02", NULL);
	TwAddVarRW(initial_bar, "Input file name:", TW_TYPE_CSSTRING(sizeof(inputfilename)), inputfilename, "");
	TwAddButton(initial_bar, "Read from CSV", CB_ReadCSV, NULL, "label='read from *.csv file' ");
	TwAddButton(initial_bar, "Read from OVF", CB_ReadOVF, NULL, "label='read from *.ovf file' ");
    TwAddButton(initial_bar, "Read from VTK", CB_ReadVTK, NULL, "label='read from *.vtk file' ");
    TwAddButton(initial_bar, "Read from BIN", CB_ReadBIN, NULL, "label='read from *.bin file' ");

	TwAddSeparator(initial_bar, "sep2", NULL);
	TwAddVarRW(initial_bar, "Save slice", TW_TYPE_BOOL32, &save_slice, " label='Save slice' help='Save current slice only' ");
	{
	TwEnumVal		enAverage_mode[] = {{ALONG_A, 	"Along a-axis"	}, 
										{ALONG_B, 	"Along b-axis"	}, 
										{ALONG_C, 	"Along c-axis"	}, 
										{ALONG_0, 	"No averaging"	}
									};

	TwType			TV_TYPE_AVERAGE_MODE = TwDefineEnum("Average mode", enAverage_mode, 4);
	TwAddVarRW(initial_bar, "Averaging mode", TV_TYPE_AVERAGE_MODE, &WhichAverageMode, "help='Choose type of average mode'");}

	TwAddVarRW(initial_bar, "Output file name:", TW_TYPE_CSSTRING(sizeof(outputfilename)), outputfilename, ""); 
	TwAddButton(initial_bar, "Write to CSV", CB_SaveCSV, NULL, "label='write to *.csv file' ");	
	TwAddButton(initial_bar, "Write to OVF", CB_Save_OVF_b8, NULL, "label='write to *.ovf file' ");	
    TwAddButton(initial_bar, "Write to VTK", CB_Save_VTK_b4, NULL, "label='write to *.vtk file' "); 
    TwAddButton(initial_bar, "Write to BIN", CB_Save_BIN, NULL, "label='write to *.bin file' "); 
    TwAddButton(initial_bar, "Write to BMP", CB_Save_BMP, NULL, "label='write to *.bmp file' "); 


    TwAddSeparator(initial_bar, "sep_isoline", NULL);

    




/*  AC field F5 */
	ac_field_bar = TwNewBar("AC_Field");
	TwDefine(" AC_Field iconified=true "); 
	TwDefine(" AC_Field size='220 530' color='100 70 70'  alpha=200"); // change default tweak bar size and color
	TwDefine(" AC_Field help='F5: show/hide AC Field bar' "); // change default tweak bar size and color
	TwAddVarRW(ac_field_bar, "AC field on/off", TW_TYPE_BOOL32, &AC_FIELD_ON, "keyIncr='f' label='AC/DC on/off' true='on' false='off' help='On/off ac field'");
    TwAddVarRW(ac_field_bar, "Mode recording on/off", TW_TYPE_BOOL32, &AC_MODE_REC, "label='Mode Rec. on/off' true='on' false='off' help='On/Off recording dynamical mode'");
    TwAddVarRW(ac_field_bar, "Average Dynamical Mode", TW_TYPE_INT32, &rec_num_mode, "label='Num Mode Average' help='Number over which dynamical mode is averaged'");
    TwAddVarCB(ac_field_bar, "NumImages", TW_TYPE_INT32, CB_SetNumImages, CB_GetNumImages,  &num_images, "min=1 help='period of sin-field or width of gaussian pulse field' ");

	TwAddVarRW(ac_field_bar, "AC field", TW_TYPE_FLOAT, &Hac,	"label='AC field'  help='The value of AC field amplitude' "); 

	TwAddVarRW(ac_field_bar, "ACfieldDir", TW_TYPE_DIR3F, &VHac, "label='Field direction' opened=true help='Change the direction of applied field' ");

	TwAddSeparator(ac_field_bar, "sep0", NULL);

	temp_color[0] = 55;
	temp_color[1] = 55;
	temp_color[2] = 155;

	TwSetParam(ac_field_bar, "ACfieldDir", "arrowcolor", TW_PARAM_INT32, 3, temp_color);
	{
    TwEnumVal       enenACFieldTw[] = { {SIN_FIELD,     "AC Sin(omega*t) "      },
                                        {GAUSSIAN_FIELD,"Gaussian field pulse" },
                                        {SINC_FIELD,"Sinc(omega*[t-t_offset])"},
                                        {CIRCULAR_FIELD,"Circular field (omega*t)"}};
    TwType          TV_TYPE_INI_STATE = TwDefineEnum("AC-field", enenACFieldTw, 4);
	TwAddVarRW(ac_field_bar, "Type of AC field", TV_TYPE_INI_STATE, &WhichACField, "help='Choose type of signal for time dependent magnetic field'");
	}

	TwAddVarRW(ac_field_bar, "t_offset", TW_TYPE_FLOAT,  &t_offset, "min=0 help='offset of time scale.' ");
	TwAddVarRW(ac_field_bar, "pulse width", TW_TYPE_FLOAT,  &GPulseWidth, "min=0 help='width of Gaussian pulse.' ");
	TwAddVarCB(ac_field_bar, "Period/Width", TW_TYPE_DOUBLE, CB_SetACPeriod, CB_GetACPeriod,  &Period_dc, "min=0 help='period of sin-field or width of gaussian pulse field' ");
	TwAddVarCB(ac_field_bar, "Omega=2*pi*P", TW_TYPE_DOUBLE, CB_SetOmega, CB_GetOmega,  &Omega_dc, "min=0 help='period of sin-field or width of gaussian pulse field' ");


/*  Info bar F11 */
	info_bar = TwNewBar("Info");
	TwDefine(" Info refresh=0.5 ");
	TwDefine(" Info iconified = false movable = false alwaysbottom=true resizable=false fontstyle=default fontsize=2"); 
	TwDefine(" Info help='F11: show/hide info-bar' "); // change default tweak bar size and color
	TwDefine(" Info color='10 10 10' alpha=0 "); // change default tweak bar size and color
	TwDefine(" Info help='F11: show/hide info-bar' "); // change default tweak bar size and color
	TwDefine(" Info position = '1 30' size ='170 620' valueswidth=130"); // change default tweak bar size and color
	TwAddVarRO(info_bar, "Run/Stop", TW_TYPE_BOOL32,  &Play, "true='RUNING' false='STOPED' ");
	TwAddVarRO(info_bar, "RECORD", TW_TYPE_BOOL32,  &Record, "true='On' false='Off' ");
	TwAddSeparator(info_bar, "sep+21", NULL);
	//TwAddVarRO(info_bar, "AC field", TW_TYPE_BOOL32,  &AC_FIELD_ON, "true='On' false='Off' help='AC filed on/off'");
	TwAddButton(info_bar, "AC field", NULL, NULL, " ");
	TwAddVarRO(info_bar, " Bx ", TW_TYPE_DOUBLE,  &Bac[0], " ");
	TwAddVarRO(info_bar, " By ", TW_TYPE_DOUBLE,  &Bac[1], " ");
	TwAddVarRO(info_bar, " Bz ", TW_TYPE_DOUBLE,  &Bac[2], " ");

	TwAddSeparator(info_bar, "sep+11", NULL);

	//TwAddVarRO(info_bar, "DC field", TW_TYPE_FLOAT, NULL, " label='DC field'");
	TwAddButton(info_bar, "DC field", NULL, NULL, " ");
	TwAddVarRO(info_bar, " Bx", TW_TYPE_FLOAT,  &Bdc[0], " ");
	TwAddVarRO(info_bar, " By", TW_TYPE_FLOAT,  &Bdc[1], " ");
	TwAddVarRO(info_bar, " Bz", TW_TYPE_FLOAT,  &Bdc[2], " ");

	TwAddSeparator(info_bar, "sep-0", NULL);
	TwAddVarRO(info_bar, "NPB", TW_TYPE_INT32,  &AtomsPerBlock, "help='number of atoms per block' ");
	TwAddVarRO(info_bar, "N_a", TW_TYPE_INT32,  &uABC[0],       "help='translations along a' ");
	TwAddVarRO(info_bar, "N_b", TW_TYPE_INT32,  &uABC[1],       "help='translations along b' ");
	TwAddVarRO(info_bar, "N_c", TW_TYPE_INT32,  &uABC[2],       "help='translations along c' ");	
	TwAddVarRO(info_bar, "NOS", TW_TYPE_INT32,  &NOS,           "help='Number of spins'      ");
	TwAddSeparator(info_bar, "sep", NULL);
	TwAddVarRO(info_bar, "ITR", TW_TYPE_INT32,  &currentIteration, "help='Total number of iterations' ");


	TwAddSeparator(info_bar, "sep0", NULL);

	TwAddVarRO(info_bar, "FPS", TW_TYPE_FLOAT,  &FPS, "help='Frame per second' precision=4");
	TwAddVarRO(info_bar, "IPS", TW_TYPE_FLOAT,  &IPS, "help='Iterations per secon' precision=4");	

	TwAddSeparator(info_bar, "sep01", NULL);


	TwAddVarRO(info_bar, "Etot", TW_TYPE_DOUBLE, &totalEnergy, " precision=10 help='Total energy' ");
	TwAddVarRO(info_bar, "etot", TW_TYPE_DOUBLE, &perSpEnergy, " precision=10 help='Energy density per spin'");
	TwAddVarRO(info_bar, "esat", TW_TYPE_DOUBLE, &totalEnergyFerro, " precision=10 help='Energy density for ferromagnetic state' ");
	TwAddVarRO(info_bar, "de", TW_TYPE_DOUBLE, &perSpEnergyMinusFerro, " precision=10 help='Energy density per spin wrt ferromagnetic state'");

	TwAddSeparator(info_bar, "sep1", NULL);	
	TwAddVarRO(info_bar, "M_x", TW_TYPE_DOUBLE, &Mtot[0], " help='x-component of total moment' precision=10");
	TwAddVarRO(info_bar, "M_y", TW_TYPE_DOUBLE, &Mtot[1], " help='y-component of total moment' precision=10");
	TwAddVarRO(info_bar, "M_z", TW_TYPE_DOUBLE, &Mtot[2], " help='z-component of total moment' precision=10");

	TwAddSeparator(info_bar, "sep2", NULL);
	TwAddVarRO(info_bar, "m_x", TW_TYPE_DOUBLE, &mtot[0], " help='x-component of average moment per spin' precision=10");
	TwAddVarRO(info_bar, "m_y", TW_TYPE_DOUBLE, &mtot[1], " help='y-component of average moment per spin' precision=10");
	TwAddVarRO(info_bar, "m_z", TW_TYPE_DOUBLE, &mtot[2], " help='z-component of average moment per spin' precision=10");
	TwAddSeparator(info_bar, "sep3", NULL);
	TwAddVarRO(info_bar, "Torque", TW_TYPE_DOUBLE, &MAX_TORQUE, " help='maximum torque acting on the spin' precision=10");
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
			TransXYZ[2]+=0.5;	break;

			case 4:		                // scrol down
			TransXYZ[2]-=0.5;	break;

			default:
				b = 0;
				//fprintf( stderr, "Unknown mouse button: %d\n", button );
		}

		// button down sets the bit, up clears the bit:
		if( state == GLUT_DOWN )
		{
			Xmouse = x;
			Ymouse = y;
			ActiveButton |= b;		// set the proper bit
		}else{
			ActiveButton &= ~b;		// clear the proper bit
		}
	}
}


void MouseMotion( int x, int y ) // called when the mouse moves while a button is down
{
	const int 	dx = x - Xmouse;// change in mouse coords
	const int 	dy = y - Ymouse;
	float 	quat1[4];//NSK
	float 	quat2[4];//NSK
if( !TwEventMouseMotionGLUT(x, y) )  // send event to AntTweakBar
    { // event has not been handled by AntTweakBar
      // your code here to handle the event
      if( ( ActiveButton & LEFT ) != 0 )
		{
			// Rot[2] += ( dx )*0.1f;
			// Rot[0] += ( dy )*0.1f;
			// RotAxis[0] = dy;
			// RotAxis[1] = dx;
			if (angle!=0){
			// SetQuaternionFromAxisAngle(axis, angle, quat);
			// MultiplyQuaternions(q_Rotation, quat, q_Rotation);
			SetQuaternionFromAxisAngle(axisY, dx*0.01, quat2);
			SetQuaternionFromAxisAngle(axisX, dy*0.01, quat1);
			MultiplyQuaternions(quat2, quat1, quat1);
			MultiplyQuaternions( quat1, q_Rotation, q_Rotation);
            //metka1 GetEulerFromQuaternion(q_Rotation, Rot);
			}		
		}

		if( ( ActiveButton & MIDDLE ) != 0 )
		{
			TransXYZ[2]+=(dx-dy)*0.05;	
		}

		if( ( ( ActiveButton & RIGHT ) != 0 ) & (true))//WhichProjection == PERSP
		{
			TransXYZ[0]+=dx*0.1;
			TransXYZ[1]-=dy*0.1;
		}

		Xmouse = x;			// new current position
		Ymouse = y;	
    }
}

// the keyboard callback:
void keyboardDown( unsigned char key, int x, int y )
{   
    if( !TwEventKeyboardGLUT(key, x, y) )  // send event to AntTweakBar 
    { // event has not been handled by AntTweakBar
      // your code here to handle the event	
	  // if( false ) fprintf( stderr, "Keyboard: '%c' (0x%0x)\n", c, c );
		switch( key )
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
				dRot[2] = - RotSpeed*0.5;
			break;

			case 'e':
			case 'E':
				dRot[2] = RotSpeed*0.5;
			break;

			case 'w':
			case 'W':
			    dRot[0] = - RotSpeed*0.5;
				break;
			case 's':
			case 'S':
				dRot[0] = RotSpeed*0.5;
				//TransXYZ[2]+=0.5;
				break;
			case 'a':
			case 'A':
				dRot[1] = - RotSpeed*0.5;
				//TransXYZ[0]-=1;	
				break;
			case 'd':
			case 'D':
				dRot[1] = RotSpeed*0.5;
				//TransXYZ[0]+=1;	
				break;

			case 'G':
			case 'g':
				dTransXYZ[2] = -TransSpeed*0.1;	
			break;

			case 'T':
			case 't':
				dTransXYZ[2] = TransSpeed*0.1;
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

			//default:
				//fprintf( stderr, "Don't know what to do with keyboard hit: '%c' (0x%0x)\n", c, c );
		}

	}
}


void keyboardUp( unsigned char key, int x, int y )
{
		switch( key )
		{			

			case 'q':
			case 'Q':
				dRot[2] = 0.0f;
			break;

			case 'e':
			case 'E':
				dRot[2] = 0.0f;
			break;

			case 'w':
			case 'W':
				//Rot[0] -= 0.25;//NSK
			    dRot[0] = 0.0f;
				//TransXYZ[2]-=0.5;
				break;
			case 's':
			case 'S':
				dRot[0] = 0.0f;
				//TransXYZ[2]+=0.5;
				break;
			case 'a':
			case 'A':
				dRot[1] = 0.0f;
				//TransXYZ[0]-=1;	
				break;
			case 'd':
			case 'D':
				dRot[1] = 0.0f;
				//TransXYZ[0]+=1;	
				break;

			case 'G':
			case 'g':
				dTransXYZ[2] = 0.0f;
			break;


			case 'T':
			case 't':
				dTransXYZ[2] = 0.0f;	
			break;
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
			case  GLUT_KEY_F4:
				TwGetParam(control_bar, NULL, "iconified", TW_PARAM_INT32, 1, &isiconified);
				if (isiconified){
					TwDefine(" Parameters&Controls iconified=false ");				
				}else{
					TwDefine(" Parameters&Controls iconified=true ");
				}
				break;
			case  GLUT_KEY_F6:
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
			//default:
				//fprintf( stderr, "Don't know what to do with special key: '%c' (0x%0x)\n", key, key );
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

		// default:
		// 	fprintf( stderr, "Don't know what to do with Button ID %d\n", id );
	}
}

void ChangeInitialState ( int id )
{
	InitSpinComponents( Px, Py, Pz, Sx, Sy, Sz, id );
	for (int i=0;i<NOS;i++) { bSx[i]=Sx[i]; bSy[i]=Sy[i]; bSz[i]=Sz[i];}
	ChangeVectorMode ( 1 );
	SpecialEvent=1;
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
		HueMap[0]=HueMapRGB[0]+ColorShift;
		HueMap[1]=HueMapRGB[1]+ColorShift;
		HueMap[2]=HueMapRGB[2]+ColorShift;
		HueMap[3]=HueMapRGB[3]+ColorShift;
		HueMap[4]=HueMapRGB[4]+ColorShift;
		HueMap[5]=HueMapRGB[5]+ColorShift;
		// InitRGB(RHue, GHue, BHue, HueMap);
		case 1:
		ChangeVectorMode(1);
	}
}

void ReallocateArrayDrawing()
{
	free(vertexProto); free(normalProto); free(indicesProto);
	free(vertices); free(normals); free(colors); free(indices);

	int NOS_L=0;
	int NOB_L=0;
  
	switch( WhichSliceMode)  // see also UpdateIndices()
	{
		case A_AXIS:
		NOS_L=NOS_AL * (1+A_layer_max - A_layer_min);
		if (A_layer_max - A_layer_min<2){
			NOB_L=NOB_AL * (1+A_layer_max - A_layer_min);
		}else{
			NOB_L=2*NOB_AL+2*(A_layer_max - A_layer_min-1)*(uABC[1]-1+uABC[2]-1);
            if (uABC[2]==1 || uABC[1]==1){
                NOB_L=NOB_AL * (1+A_layer_max - A_layer_min);
            } 
		}
		break;
		case B_AXIS:
		NOS_L=NOS_BL * (1+B_layer_max - B_layer_min);
		if (B_layer_max - B_layer_min<2){
			NOB_L=NOB_BL * (1+B_layer_max - B_layer_min);
		}else{
			NOB_L=2*NOB_BL + 2*(B_layer_max - B_layer_min-1)*(uABC[0]-1+uABC[2]-1);
            if (uABC[0]==1 || uABC[2]==1){
                NOB_L=NOB_BL * (1+B_layer_max - B_layer_min);
            } 
		}
		break;
		case C_AXIS:
		NOS_L=NOS_CL * (1+C_layer_max - C_layer_min);
		if (C_layer_max - C_layer_min<2){
			NOB_L=NOB_CL * (1+C_layer_max - C_layer_min);
		}else{
			NOB_L=2*NOB_CL + 2*(C_layer_max - C_layer_min-1)*(uABC[0]-1+uABC[1]-1);
            if (uABC[0]==1 || uABC[1]==1){
                NOB_L=NOB_CL * (1+C_layer_max - C_layer_min);
            } 
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

void ReallocateArrayDrawing_AC_phase()
{
    free(vertices_AC_phase); free(colors_AC_phase); free(indices_AC_phase);         
    IdNum_AC_phase = 2; // number of indixes
    VCNum_AC_phase = 3*2; // metka
    vertices_AC_phase  = (float  *)malloc(VCNum_AC_phase * sizeof( float  ));
    colors_AC_phase    = (float  *)malloc(VCNum_AC_phase * sizeof( float  ));
    indices_AC_phase   = (GLuint *)malloc(IdNum_AC_phase * sizeof( GLuint ));             
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
    // NumberElementInSlice(NOS_L, NOB_L);
    switch( WhichSliceMode)
    {
        case A_AXIS:
        NOS_L=NOS_AL * (1+A_layer_max - A_layer_min);
        if (A_layer_max - A_layer_min<2){
            NOB_L=NOB_AL * (1+A_layer_max - A_layer_min);
        }else{
            NOB_L=2*NOB_AL+2*(A_layer_max - A_layer_min-1)*(uABC[1]-1+uABC[2]-1);
            if (uABC[2]==1 || uABC[1]==1){
                NOB_L=NOB_AL * (1+A_layer_max - A_layer_min);
            } 
        }
        break;
        case B_AXIS:
        NOS_L=NOS_BL * (1+B_layer_max - B_layer_min);
        if (B_layer_max - B_layer_min<2){
            NOB_L=NOB_BL * (1+B_layer_max - B_layer_min);
        }else{
            NOB_L=2*NOB_BL + 2*(B_layer_max - B_layer_min-1)*(uABC[0]-1+uABC[2]-1);
            if (uABC[0]==1 || uABC[2]==1){
                NOB_L=NOB_BL * (1+B_layer_max - B_layer_min);
            } 
        }
        break;
        case C_AXIS:
        NOS_L=NOS_CL * (1+C_layer_max - C_layer_min);
        if (C_layer_max - C_layer_min<2){
            NOB_L=NOB_CL * (1+C_layer_max - C_layer_min);
        }else{
            NOB_L=2*NOB_CL + 2*(C_layer_max - C_layer_min-1)*(uABC[0]-1+uABC[1]-1);
            if (uABC[0]==1 || uABC[1]==1){
                NOB_L=NOB_CL * (1+C_layer_max - C_layer_min);
            } 
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
	bool F;//factor

	switch( WhichSliceMode)//see also ReallocateArrayDrawing()
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
			for (int an = 0; an<uABC[0]; an++) {
			for (int bn = 0; bn<uABC[1]; bn++) {
			for (int cn = 0; cn<uABC[2]; cn++) {
				n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
				n = n*AtomsPerBlock;//index of the first spin in the block
				if (GreedFilter){
					F=(an>=GreedFilterMinA && an<=GreedFilterMaxA) &&
						  (bn>=GreedFilterMinB && bn<=GreedFilterMaxB) &&
						  (cn>=GreedFilterMinC && cn<=GreedFilterMaxC) ;
					if (GreedFilterInvert) F=!F;					
				}else{F=true;}


				if (F)
				{
					for (int atom=0; atom<AtomsPerBlock; atom++)//atom is the index of the atom in block
					{
					    N = n + atom;
						//slow version is commented but easy to read: 
						S[0] = Sx[N];
						S[1] = Sy[N];
						S[2] = Sz[N];
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
						F =	(SpinFilter1 && ((S[2]>=Sz_min1 && S[2]<=Sz_max1) && (phi>=phi_min1 && phi<=phi_max1)) ) ||
							(SpinFilter2 && ((S[2]>=Sz_min2 && S[2]<=Sz_max2) && (phi>=phi_min2 && phi<=phi_max2)) ) ||
							(SpinFilter3 && ((S[2]>=Sz_min3 && S[2]<=Sz_max3) && (phi>=phi_min3 && phi<=phi_max3)) ) ;													

						if (F)
						{
							HSVtoRGB(S, RGB, InvertValue, InvertHue);//metka
                            vlength = Unitf(S,S);
							N_filter++;
					        // HSVtoRGB(S, RGB, InvertValue, InvertHue);
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

								Vout[i+0] =Kind[N]*( (-S[1]*A + Vinp[3*k+0]*S[2] + S[0]*Vinp[3*k+2]			  )*Scale*vlength + Px[N]);
								Vout[i+1] =Kind[N]*( ( S[0]*A + Vinp[3*k+1]*S[2] + S[1]*Vinp[3*k+2]			  )*Scale*vlength + Py[N]);
								Vout[i+2] =Kind[N]*( ( Vinp[3*k+2]*S[2] - (S[0]*Vinp[3*k+0]+S[1]*Vinp[3*k+1]) )*Scale*vlength + Pz[N]);	

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
                    HSVtoRGB(S, RGB, InvertValue, InvertHue);//metka
					vlength = Unitf(S,S);
			        // HSVtoRGB( S, RGB, InvertValue, InvertHue);
					j++;
					if (S[2]==-1){
					for (int k=0; k<Kinp/3; k++){// k runs over vertices 
							i = j*Kinp + 3*k;
							Vout[i+0] = Kind[N]*((-Vinp[3*k+0])*Scale*vlength + Px[N]);
							Vout[i+1] = Kind[N]*(( Vinp[3*k+1])*Scale*vlength + Py[N]);
							Vout[i+2] = Kind[N]*((-Vinp[3*k+2])*Scale*vlength + Pz[N]);	

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

							Vout[i+0] =Kind[N]*((-S[1]*A + Vinp[3*k+0]*S[2] + S[0]*Vinp[3*k+2]			 )*Scale*vlength + Px[N]);
							Vout[i+1] =Kind[N]*(( S[0]*A + Vinp[3*k+1]*S[2] + S[1]*Vinp[3*k+2]			 )*Scale*vlength + Py[N]);
							Vout[i+2] =Kind[N]*(( Vinp[3*k+2]*S[2] - (S[0]*Vinp[3*k+0]+S[1]*Vinp[3*k+1]) )*Scale*vlength + Pz[N]);

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
				if (GreedFilter){
						F =	(an>=GreedFilterMinA && an<=GreedFilterMaxA) &&
							(bn>=GreedFilterMinB && bn<=GreedFilterMaxB) &&
							(cn>=GreedFilterMinC && cn<=GreedFilterMaxC) ;
						if (GreedFilterInvert) F=!F;					
					}else{F=true;}

				if (F)
				{
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
					//if ( (S[2]>=Sz_min1 && S[2]<=Sz_max1) && (phi>=phi_min1 && phi<=phi_max1) )
					F =	(SpinFilter1 && ((S[2]>=Sz_min1 && S[2]<=Sz_max1) && (phi>=phi_min1 && phi<=phi_max1)) ) ||
						(SpinFilter2 && ((S[2]>=Sz_min2 && S[2]<=Sz_max2) && (phi>=phi_min2 && phi<=phi_max2)) ) ||
						(SpinFilter3 && ((S[2]>=Sz_min3 && S[2]<=Sz_max3) && (phi>=phi_min3 && phi<=phi_max3)) ) ;

					if (F)
					{
						N_filter++;
                        HSVtoRGB(S, RGB, InvertValue, InvertHue);//metka
					    (void)Unitf(S,S);
					    // HSVtoRGB( S, RGB, InvertValue, InvertHue);

						j++;
						for (int k=0; k<Kinp/3; k++) // k runs over vertices of the box 
						{
							i = j*Kinp + 3*k;	// vertex index
                            N = n*AtomsPerBlock;//first index of the atom in the blok defines visible/invisible
							Vout[i+0] = (Vinp[3*k+0] + BPx[n])*Kind[N];
							Vout[i+1] = (Vinp[3*k+1] + BPy[n])*Kind[N];
							Vout[i+2] = (Vinp[3*k+2] + BPz[n])*Kind[N];	

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

					n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
					for (int atom=0; atom<AtomsPerBlock; atom++)//atom is the index of the atom in block
					{
					    N = atom + n*AtomsPerBlock;//n*AtomsPerBlock=index of the first spin in the block;
		 
						S[0]+= Sx[N];
						S[1]+= Sy[N];
						S[2]+= Sz[N];
				    }
                    HSVtoRGB(S, RGB, InvertValue, InvertHue);//metka
				    (void)Unitf(S,S);
				    // HSVtoRGB( S, RGB, InvertValue, InvertHue);

					j++;
					for (int k=0; k<Kinp/3; k++) // k runs over vertices of the box 
					{
						i = j*Kinp + 3*k;	// vertex index
                        N=n*AtomsPerBlock;//first index of the atom in the blok defines visible/invisible
						Vout[i+0] = (Vinp[3*k+0] + BPx[n])*Kind[N];
						Vout[i+1] = (Vinp[3*k+1] + BPy[n])*Kind[N];
						Vout[i+2] = (Vinp[3*k+2] + BPz[n])*Kind[N];	

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
				if (GreedFilter){
					F =	(an>=GreedFilterMinA && an<=GreedFilterMaxA) &&
						(bn>=GreedFilterMinB && bn<=GreedFilterMaxB) &&
						(cn>=GreedFilterMinC && cn<=GreedFilterMaxC) ;
					if (GreedFilterInvert) F=!F;					
				}else{F=true;}

				if (F)
				{
					n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
					n = n*AtomsPerBlock;//index of the first spin in the block
					for (int atom=0; atom<AtomsPerBlock; atom++)//atom is the index of the atom in block
					{
						S[0] = Sx[n+atom];
						S[1] = Sy[n+atom];
						S[2] = Sz[n+atom];
                        HSVtoRGB( S, RGB, InvertValue, InvertHue);
						vlength = Unitf(S,S);//metka
						int phi = atan2int( S[1], S[0] );// return integer angle phi 0 - 360
						//if (S[2]>=Sz_min1 && S[2]<=Sz_max1)
						//if ( (S[2]>=Sz_min1 && S[2]<=Sz_max1) && (phi>=phi_min1 && phi<=phi_max1) )
						F =	(SpinFilter1 && ((S[2]>=Sz_min1 && S[2]<=Sz_max1) && (phi>=phi_min1 && phi<=phi_max1)) ) ||
							(SpinFilter2 && ((S[2]>=Sz_min2 && S[2]<=Sz_max2) && (phi>=phi_min2 && phi<=phi_max2)) ) ||
							(SpinFilter3 && ((S[2]>=Sz_min3 && S[2]<=Sz_max3) && (phi>=phi_min3 && phi<=phi_max3)) ) ;
						if (F)
						{
							N_filter++;
					        j++;
							i = j*Kinp;			// index of first cane vertex 
							Vout[i+0] = Px[n+atom];	// new x-component of vertex + translation
							Vout[i+1] = Py[n+atom];	// new y-component of vertex + translation
							Vout[i+2] = Pz[n+atom];	// new z-component of vertex + translation
							Cout[i+0] = RGB[0]*Kind[n+atom];	// x-component of vertex normal
							Cout[i+1] = RGB[1]*Kind[n+atom];	// y-component of vertex normal
							Cout[i+2] = RGB[2]*Kind[n+atom];	// z-component of vertex normal
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
				n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
				n = n*AtomsPerBlock;//index of the first spin in the block
				for (int atom=0; atom<AtomsPerBlock; atom++)//atom is the index of the atom in block
				{
					S[0] = Sx[n+atom];
					S[1] = Sy[n+atom];
					S[2] = Sz[n+atom];
                    HSVtoRGB(S, RGB, InvertValue, InvertHue);//metka
					vlength = Unitf(S,S);
			        // HSVtoRGB( S, RGB, InvertValue, InvertHue);
			        j++;
					i = j*Kinp;			// index of first cane vertex 
					int Factor = Kind[n+atom];
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
		if (WhichSliceMode==FILTER)
		{
			N_filter=0;
			for (int an = 0; an<uABC[0]; an++) // an runs over slice along A_AXIS
			{
			for (int bn = 0; bn<uABC[1]; bn++) // bn runs over slice along B_AXIS
			{
			for (int cn = cnini; cn<uABC[2]; cn++) // cn runs over slice along C_AXIS
			{
				if (GreedFilter){
					F =	(an>=GreedFilterMinA && an<=GreedFilterMaxA) &&
						(bn>=GreedFilterMinB && bn<=GreedFilterMaxB) &&
						(cn>=GreedFilterMinC && cn<=GreedFilterMaxC) ;
					if (GreedFilterInvert) F=!F;					
				}else{F=true;}

				if (F)
				{
					n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
					n = n*AtomsPerBlock;//index of the first spin in the block
					for (int atom=0; atom<AtomsPerBlock; atom++)//atom is the index of the atom in block
					{
						S[0] = Sx[n+atom];
						S[1] = Sy[n+atom];
						S[2] = Sz[n+atom];
						int phi = atan2int( S[1], S[0] );// return integer angle phi 0 - 360
						//if (S[2]>=Sz_min1 && S[2]<=Sz_max1)
						//if ( (S[2]>=Sz_min1 && S[2]<=Sz_max1) && (phi>=phi_min1 && phi<=phi_max1) )
						F =	(SpinFilter1 && ((S[2]>=Sz_min1 && S[2]<=Sz_max1) && (phi>=phi_min1 && phi<=phi_max1)) ) ||
							(SpinFilter2 && ((S[2]>=Sz_min2 && S[2]<=Sz_max2) && (phi>=phi_min2 && phi<=phi_max2)) ) ||
							(SpinFilter3 && ((S[2]>=Sz_min3 && S[2]<=Sz_max3) && (phi>=phi_min3 && phi<=phi_max3)) ) ;	
						if (F)
						{
                            HSVtoRGB( S, RGB, InvertValue, InvertHue);
                            vlength = Unitf(S,S);//metka
							N_filter++;
					        j++;
							//i = (n-nini)*Kinp;							// index of ferst cane vertex 
							i = j*Kinp;
							Vout[i+0] = Kind[n+atom]*( S[0]*(1-Pivot)*Scale*vlength + Px[n+atom]);	// new x-component of vertex + translation
							Vout[i+1] = Kind[n+atom]*( S[1]*(1-Pivot)*Scale*vlength + Py[n+atom]);	// new y-component of vertex + translation
							Vout[i+2] = Kind[n+atom]*( S[2]*(1-Pivot)*Scale*vlength + Pz[n+atom]);	// new z-component of vertex + translation
							//i = n*Kinp/3*4;		// colors contains 4 floats
							Cout[i+0] = RGB[0];					// x-component of vertex normal
							Cout[i+1] = RGB[1];					// y-component of vertex normal
							Cout[i+2] = RGB[2];					// z-component of vertex normal
							//Cout[i+3] = 1.f;
							//printf( "|V1=%f,%f,%f \n",Vout[i+0],Vout[i+1],Vout[i+2]);
							//i = (n-nini)*Kinp + 3*1;					// index of ferst cane vertex 
							i = j*Kinp+ 3*1;
							Vout[i+0] = Kind[n+atom]*(-S[0]*(Pivot)*Scale*vlength + Px[n+atom]);		// new x-component of vertex + translation
							Vout[i+1] = Kind[n+atom]*(-S[1]*(Pivot)*Scale*vlength + Py[n+atom]);		// new y-component of vertex + translation
							Vout[i+2] = Kind[n+atom]*(-S[2]*(Pivot)*Scale*vlength + Pz[n+atom]);		// new z-component of vertex + translation
							//i = n*Kinp/3*4+4;			// colors contains 4 floats
							Cout[i+0] = RGB[0];//*Kind[n+atom];					// x-component of vertex normal
							Cout[i+1] = RGB[1];//*Kind[n+atom];					// y-component of vertex normal
							Cout[i+2] = RGB[2];//*Kind[n+atom];					// z-component of vertex normal
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
				n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
				n = n*AtomsPerBlock;//index of the first spin in the block
				for (int atom=0; atom<AtomsPerBlock; atom++)//atom is the index of the atom in block
				{
					S[0] = Sx[n+atom];
					S[1] = Sy[n+atom];
					S[2] = Sz[n+atom];
                    HSVtoRGB( S, RGB, InvertValue, InvertHue);//metka
					vlength = Unitf(S,S);
			        // HSVtoRGB( S, RGB, InvertValue, InvertHue);
			        j++;
					//i = (n-nini)*Kinp;							// index of ferst cane vertex 
					i = j*Kinp;//intMAX
					int Factor = Kind[n+atom];
					if (Factor==0) Factor=intMAX;
					Vout[i+0] = Factor*( S[0]*(1-Pivot)*Scale*vlength + Px[n+atom]);	// new x-component of vertex + translation
					Vout[i+1] = Factor*( S[1]*(1-Pivot)*Scale*vlength + Py[n+atom]);	// new y-component of vertex + translation
					Vout[i+2] = Factor*( S[2]*(1-Pivot)*Scale*vlength + Pz[n+atom]);	// new z-component of vertex + translation
					//i = n*Kinp/3*4;		// colors contains 4 floats
					Cout[i+0] = RGB[0];					// x-component of vertex normal
					Cout[i+1] = RGB[1];					// y-component of vertex normal
					Cout[i+2] = RGB[2];					// z-component of vertex normal
					//Cout[i+3] = 1.f;
					//printf( "|V1=%f,%f,%f \n",Vout[i+0],Vout[i+1],Vout[i+2]);
					//i = (n-nini)*Kinp + 3*1;					// index of ferst cane vertex 
					i = j*Kinp+ 3*1;
					Vout[i+0] = Factor*(-S[0]*(Pivot)*Scale*vlength + Px[n+atom]);		// new x-component of vertex + translation
					Vout[i+1] = Factor*(-S[1]*(Pivot)*Scale*vlength + Py[n+atom]);		// new y-component of vertex + translation
					Vout[i+2] = Factor*(-S[2]*(Pivot)*Scale*vlength + Pz[n+atom]);		// new z-component of vertex + translation
					//i = n*Kinp/3*4+4;			// colors contains 4 floats
					Cout[i+0] = Kind[n+atom]*RGB[0];					// x-component of vertex normal
					Cout[i+1] = Kind[n+atom]*RGB[1];					// y-component of vertex normal
					Cout[i+2] = Kind[n+atom]*RGB[2];					// z-component of vertex normal
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
							float Vx, float Vy, float Vz)
{
	int i;
    float Sx=Vx;
    float Sy=Vy;
    float Sz=Vz;
	float U,A;
	if (Sz==-1){
		for (int k=0; k<Kinp/3; k++){// k runs over vertices 
			i = 3*k;	// vertex index
			Vout[i+0] = (-Vinp[i+0])*Hf*Scale_H + Px;
			Vout[i+1] = ( Vinp[i+1])*Hf*Scale_H + Py;
			Vout[i+2] = (-Vinp[i+2])*Hf*Scale_H + Pz;	
		
			Nout[i+0] = -Ninp[i+0];
			Nout[i+1] =  Ninp[i+1];
			Nout[i+2] = -Ninp[i+2];
			if (WhichBackgroundColor == BLACK){
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
			if (WhichBackgroundColor == BLACK){
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

void UpdateVerticesNormalsColors_AC_phase(float * vertices, float * colors, GLuint * indices)
{
    int i=-1;
    int j=-1;
    vertices[++i]=0; colors[i]=1;
    vertices[++i]=0; colors[i]=1;
    vertices[++i]=0; colors[i]=1;
    indices[++j]=0;

    vertices[++i]=200; colors[i]=1;
    vertices[++i]=100; colors[i]=1;
    vertices[++i]=0;   colors[i]=1;
    indices[++j]=1;

    // float   d = WireWidth;
    // float   cube[3][3] = {
    //             {   1.0f, 0.0f, 0.0f }, // a
    //             {   0.0f, 1.0f, 0.0f }, // b
    //             {   0.0f, 0.0f, 1.0f }};// c
    // float tr[3] = {-d/2, -d/2, -d/2};
    // float length = 20*d;
    // float colorR[3]={1, 0, 0};
    // float colorG[3]={0, 1, 0};
    // float colorB[3]={0, 0, 1};
    // float color0[3]={0.1,0.1,0.1};
    
    // parallelepiped( cube, tr, length, d, d, colorR, 0, vertices, normals, colors, indices );//X
    // parallelepiped( cube, tr, d, length, d, colorG, 1, vertices, normals, colors, indices );//Y
    // parallelepiped( cube, tr, d, d, length, colorB, 2, vertices, normals, colors, indices );//Z
    // d *= 1.5; tr[0]=-d/2; tr[1]=-d/2; tr[2]=-d/2;

    // parallelepiped( cube, tr, d, d, d, color0, 3, vertices, normals, colors, indices );//O

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

void CreateNewVBO_AC_phase( ){
    glGenBuffers(1, &vboIdV_AC_phase);
    glGenBuffers(1, &vboIdC_AC_phase);
    glGenBuffers(1, &iboIdI_AC_phase);
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

void UpdateVBO_AC_phase(GLuint * V, GLuint * C, GLuint * I, float * ver, float * col, GLuint * ind)
{//metka   
    int VCNum=3*2;
    glBindBuffer(GL_ARRAY_BUFFER, *V);
    glBufferData(GL_ARRAY_BUFFER, VCNum * sizeof(float), 0, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0, VCNum * sizeof(float), ver); 
    // for cane and point modes we do not need normals  
    glBindBuffer(GL_ARRAY_BUFFER, *C);
    glBufferData(GL_ARRAY_BUFFER, VCNum * sizeof(float), 0, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0, VCNum * sizeof(float), col);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, *I);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, IdNum * sizeof(GLuint), 0, GL_DYNAMIC_DRAW);
    glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, IdNum * sizeof(GLuint), ind);
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

		case CANE: 
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

void drawVBO_AC_phase()
{
    glDisable(GL_LIGHTING);
    glBindBuffer(GL_ARRAY_BUFFER, vboIdC_AC_phase);      glColorPointer(3, GL_FLOAT, 0, (void*)0);
    glBindBuffer(GL_ARRAY_BUFFER, vboIdV_AC_phase);      glVertexPointer(3, GL_FLOAT, 0, (void*)0);  

    glEnableClientState(GL_COLOR_ARRAY);        // enable color arrays
    glEnableClientState(GL_VERTEX_ARRAY);       // enable vertex arrays 

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, iboIdI_AC_phase);
    
    glLineWidth(5.0f);
    //glDrawElements(GL_LINES, iNum, GL_UNSIGNED_INT, (void*)(iStart*sizeof(GLuint)));// draw a stick VCNum
    glDrawElements(GL_LINES, IdNum, GL_UNSIGNED_INT, (void*)0);// draw a stick VCNum

    glDisableClientState(GL_VERTEX_ARRAY);      // disable vertex arrays
    glDisableClientState(GL_COLOR_ARRAY);       // disable normal arrays

    glBindBuffer(GL_ARRAY_BUFFER,           0); // disable vertex arrays
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,   0); // disable normal arrays

    glBindBuffer(GL_ARRAY_BUFFER, vboIdC_AC_phase);      glColorPointer(3, GL_FLOAT, 6*sizeof(float), (void*)0);
    glBindBuffer(GL_ARRAY_BUFFER, vboIdV_AC_phase);      glVertexPointer(3, GL_FLOAT, 6*sizeof(float), (void*)0);

    glEnableClientState(GL_COLOR_ARRAY);        // enable normal arrays
    glEnableClientState(GL_VERTEX_ARRAY);       // enable vertex arrays     

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, iboIdI_AC_phase);
    
    glPointSize(10.5f*Scale);
    glDrawElements(GL_POINTS, IdNum_AC_phase/2, GL_UNSIGNED_INT, (void*)0);// draw a ball (point at the end of the stick)

    glDisableClientState(GL_VERTEX_ARRAY);      // disable vertex arrays
    glDisableClientState(GL_COLOR_ARRAY);       // disable normal arrays

    glBindBuffer(GL_ARRAY_BUFFER,           0); // disable vertex arrays
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,   0); // disable normal arrays
    glEnable(GL_LIGHTING);
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
