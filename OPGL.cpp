GLuint 			GLUT_window;
const char *	WINDOWTITLE = { "Magnoom v1.0" };
int				window_width	= 1024;
int				window_height	= 768;
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

float			PerspSet[4]		= {60.0, asp_rat, 0.1, 5000}; // {Setings: Field of view vertical, apect ratio, zNear, zFar}  

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

typedef enum	{ARROW1, CONE1, CANE, POINT, PLANE} enVectorMode; // which mode
enVectorMode	WhichVectorMode	= CANE;	// CANE by default 
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
float			Rot[3];	// rotation angles in degrees
float			TransXYZ[3];	// set by glui translation widgets
float			Scale = 1.f;	// scaling factors for arrows [0.1-2] 
float			Pivot = 0.55f;

int				Alayer = 1;	// which layer (along c tr. vect. ) to show max=ABC[2]
int				Blayer = 1;	// which layer (along c tr. vect. ) to show max=ABC[2]
int				Clayer = 1;	// which layer (along c tr. vect. ) to show max=ABC[2]
typedef enum	{A_AXIS, B_AXIS, C_AXIS, ABC_AXIS} enSliceMode; // which mode
enSliceMode	    WhichSliceMode	= C_AXIS;	// CANE by default 

float			chSize = 71; // characteristic size of initial state in units of "a"
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
float g_LightDirection[] = { 0.0f, 0.0f, -1.0f };




////////////////////////////////// Parameters for visualisation ///////////////////////////////////

GLfloat*	vertexProto   = NULL; // array of vertexes for prototipe arrow or cane
GLfloat*	normalProto   = NULL; // array of normals for prototipe arrow (not used in vector "cane mode")
GLuint*		indicesProto  = NULL; // array of indices for prototipe arrow or cane

GLfloat*	vertices	= NULL; // array of vertexes for tatal vector field 
GLfloat*	normals		= NULL; // array of normals for tatal vector field
GLfloat*	colors		= NULL; // array of colors 
GLuint*		indices		= NULL; // array of indices for tatal vector field

int			arrowFaces	= 6; // number of arrow faces, default number

GLuint		vboIdV;   // ID of VBO for vertex arrays
GLuint		vboIdN;   // ID of VBO for normal arrays
GLuint		vboIdC;   // ID of VBO for color arrays
GLuint		iboIdI;   // ID of IBO for index arrays

int			ElNumProto;   // number of triangles per arrow
int			IdNumProto;   // number of indixes per arrow
int			VCNumProto;   // number of vertex and normals component per arrow

int			ElNum; // total number of triangles for the whole vector field = ElNumProto * Number of spins
int			IdNum; // total number of indixes for the whole vector field = IdNumProto * Number of spins
int			VCNum; // total number of component of vertices for whole vector field = VCNumProto * Number of spins

int			Play=0;

int			DataTransfer=1;

void			ChangeVectorMode ( int );
void			ChangeColorMap ( int );
void			ChangeScale ( int );
void			ChangeInitialState ( int );
void			Buttons( int );
void			Keyboard( unsigned char, int, int );
void			KeyboardAdd( int, int, int );
void			MouseButton( int, int, int, int );
void			MouseMotion( int, int );
float			ElapsedSeconds( );
// color control functions
void			HSVtoRGB( float[3], float [3] , int, int);
void			InitRGB(float* , float* , float* , int*);
// VBO array preparing functions
void			ReallocateArrayDrawing(int, int);
void			UpdatePrototypeVerNorInd(float *, float *, GLuint * , int, int);
void			CreateNewVBO();
void			UpdateVBO(	GLuint * , GLuint * , GLuint * , GLuint * , float * , float * , float * , GLuint * );
void			UpdateSpinComponents(float * , float * , float * , int);
void			UpdateSpinPositions(float[][3], int[3], float[][3], int, float[3][3], float*, float*, float*);
void			InitSpinComponents(float * , float * , float * , double * , double * , double * , int N);
void			UpdateIndices(GLuint * , int, GLuint *, int, int);
void			UpdateVerticesNormalsColors(float *, float *, int, float *, float *, float *, int, float * , float * , float *, double * , double * , double *, int);
// drawing functions
void			GetBox(float[][3], int[3], float[4][3]);
void			drawVBO( int, int );
void			idle ();
void			setupTweakBar();
// return the number of seconds since the start of the program:
float
ElapsedSeconds( )	
{
	int ms = glutGet( GLUT_ELAPSED_TIME );	// get # of milliseconds since the start of the program
	return (float)ms / 1000.;				// convert it to seconds:
}

void
Resize( int window_width, int window_height) // called when user resizes the window
{
	asp_rat     = (float)((   ((double)window_width)/((double)window_height)   ));
	asp_rat_inv = (float)((   ((double)window_height)/((double)window_width)   ));

	glutSetWindow( GLUT_window );
	glViewport(0, 0, window_width, window_height);
	    // Send the new window size to AntTweakBar
    TwWindowSize(window_width, window_height);
}

void
Reset( )
{
	Rot[0] = Rot[1] = Rot[2] = TransXYZ[0] = TransXYZ[1] = TransXYZ[2] = 0.;
	PerspSet[0] =60.f;
}

void
Xup( )
{
	Rot[1] = 0.;
	Rot[0] = Rot[2] = 270.;
	TransXYZ[0] = TransXYZ[1] = 0.;
}

void
Yup( )
{
	Rot[2] = 180;
	Rot[0] = 270.;
	TransXYZ[0] = TransXYZ[1] = 0.;
}

void
Zup( )
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

void GetBox(float abc[][3], int ABC[3], float box[3][3])
{
	//origine of the box is (0,0,0)
	//three vectors, b[1]+b[2]+b[2] define main diagonal of the box:
	box[0][0] = (ABC[0])*abc[0][0];
	box[0][1] = (ABC[0])*abc[0][1];
	box[0][2] = (ABC[0])*abc[0][2];

	box[1][0] = (ABC[1])*abc[1][0];
	box[1][1] = (ABC[1])*abc[1][1];
	box[1][2] = (ABC[1])*abc[1][2];

	box[2][0] = (ABC[2])*abc[2][0];
	box[2][1] = (ABC[2])*abc[2][1];
	box[2][2] = (ABC[2])*abc[2][2];
}

void
Parallelepiped( float abc[][3], float tr[3], float scale1, float scale2, float scale3)
{
	float p0[3]={0,0,0};
	float p1[3]={abc[0][0],abc[0][1],abc[0][2]};
	(void)Unitf(p1,p1);
	p1[0]*=scale1;p1[1]*=scale1;p1[2]*=scale1;
	float p2[3]={abc[1][0],abc[1][1],abc[1][2]};
	(void)Unitf(p2,p2 );
	p2[0]*=scale2;p2[1]*=scale2;p2[2]*=scale2;
	float p3[3]={abc[2][0],abc[2][1],abc[2][2]};
	(void)Unitf(p3,p3);
	p3[0]*=scale3;p3[1]*=scale3;p3[2]*=scale3;
	float p12[3]={p1[0]+p2[0],p1[1]+p2[1],p1[2]+p2[2]};
	float p13[3]={p1[0]+p3[0],p1[1]+p3[1],p1[2]+p3[2]};
	float p23[3]={p2[0]+p3[0],p2[1]+p3[1],p2[2]+p3[2]};
	float p123[3]={p1[0]+p23[0],p1[1]+p23[1],p1[2]+p23[2]};
	float normal1[3];
	float normal2[3];
	float normal3[3];
	Enorm( p0, p2, p3, normal1);
	Enorm( p0, p1, p3, normal2);
	Enorm( p0, p1, p2, normal3);
	p0[0]+=tr[0]; p0[1]+=tr[1]; p0[2]+=tr[2];
	p1[0]+=tr[0]; p1[1]+=tr[1]; p1[2]+=tr[2];
	p2[0]+=tr[0]; p2[1]+=tr[1]; p2[2]+=tr[2];
	p3[0]+=tr[0]; p3[1]+=tr[1]; p3[2]+=tr[2];
	p12[0]+=tr[0];p12[1]+=tr[1];p12[2]+=tr[2];
	p13[0]+=tr[0];p13[1]+=tr[1];p13[2]+=tr[2];
	p23[0]+=tr[0];p23[1]+=tr[1];p23[2]+=tr[2];
	p123[0]+=tr[0];p123[1]+=tr[1];p123[2]+=tr[2];
			glNormal3fv(normal3);
			//top 
				glVertex3fv( p3);  
				glVertex3fv( p123); 
				glVertex3fv( p13); 
				glVertex3fv( p3);  
				glVertex3fv( p23); 
				glVertex3fv( p123);
			//bottom
				glVertex3fv( p0);  
				glVertex3fv( p1); 
				glVertex3fv( p12); 
				glVertex3fv( p0);  
				glVertex3fv( p12); 
				glVertex3fv( p2);	
			//front
			glNormal3fv(normal2);
				glVertex3fv( p0);  
				glVertex3fv( p3); 
				glVertex3fv( p13); 
				glVertex3fv( p0);  
				glVertex3fv( p13); 
				glVertex3fv( p1);
			//back
				glVertex3fv( p2);  
				glVertex3fv( p12); 
				glVertex3fv( p123); 
				glVertex3fv( p2);  
				glVertex3fv( p123); 
				glVertex3fv( p23);	
			//righ
			glNormal3fv(normal1);
				glVertex3fv( p1);  
				glVertex3fv( p13); 
				glVertex3fv( p123); 
				glVertex3fv( p1);  
				glVertex3fv( p123); 
				glVertex3fv( p12);
			//left
				glVertex3fv( p0);  
				glVertex3fv( p2); 
				glVertex3fv( p23); 
				glVertex3fv( p0);  
				glVertex3fv( p23); 
				glVertex3fv( p3);	
}

void
InitLists(float abc[][3], int ABC[3])
{
	float d=0.06f;
	float Tr[3] = {	
					-(abc[0][0]*ABC[0]+abc[1][0]*ABC[1]+abc[2][0]*ABC[2])/2.f,
					-(abc[0][1]*ABC[0]+abc[1][1]*ABC[1]+abc[2][1]*ABC[2])/2.f,
					-(abc[0][2]*ABC[0]+abc[1][2]*ABC[1]+abc[2][2]*ABC[2])/2.f
				  };
	Tr[0]-= d/2;
	Tr[1]-= d/2;
	Tr[2]-= d/2;
	float tr[3] = {0., 0., 0.};
	float length1 = ABC[0]*sqrt(abc[0][0]*abc[0][0]+abc[0][1]*abc[0][1]+abc[0][2]*abc[0][2]);
	float length2 = ABC[1]*sqrt(abc[1][0]*abc[1][0]+abc[1][1]*abc[1][1]+abc[1][2]*abc[1][2]);
	float length3 = ABC[2]*sqrt(abc[2][0]*abc[2][0]+abc[2][1]*abc[2][1]+abc[2][2]*abc[2][2]);
	// create the box:
	BoxList = glGenLists( 1 );
	glNewList( BoxList, GL_COMPILE );
		glBegin( GL_TRIANGLES );
			glColor3f(0.7,0.7,0.7);
			tr[0]=Tr[0]; tr[1]=Tr[1]; tr[2]=Tr[2];
			Parallelepiped( abc, tr, length1, d, d);//(0,0,0)-->(1,0,0)
			Parallelepiped( abc, tr, d, length2, d);//(0,0,0)-->(0,1,0)
			Parallelepiped( abc, tr, d, d, length3);//(0,0,0)-->(0,0,1)

			tr[0]=Tr[0]+abc[1][0]*ABC[1]; 
			tr[1]=Tr[1]+abc[1][1]*ABC[1]; 
			tr[2]=Tr[2]+abc[1][2]*ABC[1];
			Parallelepiped( abc, tr, length1, d, d);//(0,1,0)-->(1,1,0)
			tr[0]=Tr[0]+abc[2][0]*ABC[2]; 
			tr[1]=Tr[1]+abc[2][1]*ABC[2]; 
			tr[2]=Tr[2]+abc[2][2]*ABC[2];
			Parallelepiped( abc, tr, length1, d, d);//(0,0,1)-->(0,1,1)
			tr[0]+=abc[1][0]*ABC[1]; 
			tr[1]+=abc[1][1]*ABC[1]; 
			tr[2]+=abc[1][2]*ABC[1];
			Parallelepiped( abc, tr, length1, d, d);//(1,0,1)-->(1,1,1)

			tr[0]=Tr[0]+abc[0][0]*ABC[0]; 
			tr[1]=Tr[1]+abc[0][1]*ABC[0]; 
			tr[2]=Tr[2]+abc[0][2]*ABC[0];
			Parallelepiped( abc, tr, d, length2, d);//(1,0,0)-->(1,1,0)
			tr[0]=Tr[0]+abc[2][0]*ABC[2]; 
			tr[1]=Tr[1]+abc[2][1]*ABC[2]; 
			tr[2]=Tr[2]+abc[2][2]*ABC[2];
			Parallelepiped( abc, tr, d, length2, d);//(0,0,1)-->(0,1,1)
			tr[0]+=abc[0][0]*ABC[0]; 
			tr[1]+=abc[0][1]*ABC[0]; 
			tr[2]+=abc[0][2]*ABC[0];
			Parallelepiped( abc, tr, d, length2, d);//(0,1,1)-->(1,1,1)

			tr[0]=Tr[0]+abc[0][0]*ABC[0]; 
			tr[1]=Tr[1]+abc[0][1]*ABC[0]; 
			tr[2]=Tr[2]+abc[0][2]*ABC[0];
			Parallelepiped( abc, tr, d, d, length3);//(1,0,0)-->(1,0,1)
			tr[0]=Tr[0]+abc[1][0]*ABC[1]; 
			tr[1]=Tr[1]+abc[1][1]*ABC[1]; 
			tr[2]=Tr[2]+abc[1][2]*ABC[1];
			Parallelepiped( abc, tr, d, d, length3);//(0,1,0)-->(0,1,1)
			tr[0]+=abc[0][0]*ABC[0]; 
			tr[1]+=abc[0][1]*ABC[0]; 
			tr[2]+=abc[0][2]*ABC[0];
			Parallelepiped( abc, tr, d, d, length3+d);//(1,1,0)-->(1,1,1)
		glEnd( );//GL_TRIANGLES
	glEndList( );

	// create the boundary signal doted lines:
	BoundaryListA = glGenLists( 1 );
	glNewList( BoundaryListA, GL_COMPILE );
		glBegin( GL_TRIANGLES );
			glColor3f(0.7,0.7,0.7);
			tr[0]=Tr[0]+abc[0][0]*(ABC[0]+6*d); 
			tr[1]=Tr[1]+abc[0][1]*(ABC[0]+6*d);
			tr[2]=Tr[2]+abc[0][2]*(ABC[0]+6*d);
			Parallelepiped( abc, tr, 5*d, d, d);
			tr[0]=Tr[0]+abc[0][0]*(ABC[0]+16*d); 
			tr[1]=Tr[1]+abc[0][1]*(ABC[0]+16*d);
			tr[2]=Tr[2]+abc[0][2]*(ABC[0]+16*d);
			Parallelepiped( abc, tr, 5*d, d, d);
			tr[0]=Tr[0]+abc[0][0]*(-10*d); 
			tr[1]=Tr[1]+abc[0][1]*(-10*d);
			tr[2]=Tr[2]+abc[0][2]*(-10*d);
			Parallelepiped( abc, tr, 5*d, d, d);
			tr[0]=Tr[0]+abc[0][0]*(-20*d); 
			tr[1]=Tr[1]+abc[0][1]*(-20*d);
			tr[2]=Tr[2]+abc[0][2]*(-20*d);
			Parallelepiped( abc, tr, 5*d, d, d);

			tr[0]=Tr[0]+abc[1][0]*ABC[1]+abc[0][0]*(ABC[0]+6*d);
			tr[1]=Tr[1]+abc[1][1]*ABC[1]+abc[0][1]*(ABC[0]+6*d);
			tr[2]=Tr[2]+abc[1][2]*ABC[1]+abc[0][2]*(ABC[0]+6*d);
			Parallelepiped( abc, tr, 5*d, d, d);
			tr[0]=Tr[0]+abc[1][0]*ABC[1]+abc[0][0]*(ABC[0]+16*d);
			tr[1]=Tr[1]+abc[1][1]*ABC[1]+abc[0][1]*(ABC[0]+16*d);
			tr[2]=Tr[2]+abc[1][2]*ABC[1]+abc[0][2]*(ABC[0]+16*d);
			Parallelepiped( abc, tr, 5*d, d, d);
			tr[0]=Tr[0]+abc[1][0]*ABC[1]+abc[0][0]*(-10*d);
			tr[1]=Tr[1]+abc[1][1]*ABC[1]+abc[0][1]*(-10*d);
			tr[2]=Tr[2]+abc[1][2]*ABC[1]+abc[0][2]*(-10*d);
			Parallelepiped( abc, tr, 5*d, d, d);
			tr[0]=Tr[0]+abc[1][0]*ABC[1]+abc[0][0]*(-20*d);
			tr[1]=Tr[1]+abc[1][1]*ABC[1]+abc[0][1]*(-20*d);
			tr[2]=Tr[2]+abc[1][2]*ABC[1]+abc[0][2]*(-20*d);
			Parallelepiped( abc, tr, 5*d, d, d);

			tr[0]=Tr[0]+abc[2][0]*ABC[2]+abc[0][0]*(ABC[0]+6*d);
			tr[1]=Tr[1]+abc[2][1]*ABC[2]+abc[0][1]*(ABC[0]+6*d);
			tr[2]=Tr[2]+abc[2][2]*ABC[2]+abc[0][2]*(ABC[0]+6*d);
			Parallelepiped( abc, tr, 5*d, d, d);
			tr[0]=Tr[0]+abc[2][0]*ABC[2]+abc[0][0]*(ABC[0]+16*d);
			tr[1]=Tr[1]+abc[2][1]*ABC[2]+abc[0][1]*(ABC[0]+16*d);
			tr[2]=Tr[2]+abc[2][2]*ABC[2]+abc[0][2]*(ABC[0]+16*d);
			Parallelepiped( abc, tr, 5*d, d, d);
			tr[0]=Tr[0]+abc[2][0]*ABC[2]+abc[0][0]*(-10*d);
			tr[1]=Tr[1]+abc[2][1]*ABC[2]+abc[0][1]*(-10*d);
			tr[2]=Tr[2]+abc[2][2]*ABC[2]+abc[0][2]*(-10*d);
			Parallelepiped( abc, tr, 5*d, d, d);
			tr[0]=Tr[0]+abc[2][0]*ABC[2]+abc[0][0]*(-20*d);
			tr[1]=Tr[1]+abc[2][1]*ABC[2]+abc[0][1]*(-20*d);
			tr[2]=Tr[2]+abc[2][2]*ABC[2]+abc[0][2]*(-20*d);
			Parallelepiped( abc, tr, 5*d, d, d);		

			tr[0]=Tr[0]+abc[1][0]*ABC[1]+abc[2][0]*ABC[2]+abc[0][0]*(ABC[0]+6*d);
			tr[1]=Tr[1]+abc[1][1]*ABC[1]+abc[2][1]*ABC[2]+abc[0][1]*(ABC[0]+6*d);
			tr[2]=Tr[2]+abc[1][2]*ABC[1]+abc[2][2]*ABC[2]+abc[0][2]*(ABC[0]+6*d);
			Parallelepiped( abc, tr, 5*d, d, d);
			tr[0]=Tr[0]+abc[1][0]*ABC[1]+abc[2][0]*ABC[2]+abc[0][0]*(ABC[0]+16*d);
			tr[1]=Tr[1]+abc[1][1]*ABC[1]+abc[2][1]*ABC[2]+abc[0][1]*(ABC[0]+16*d);
			tr[2]=Tr[2]+abc[1][2]*ABC[1]+abc[2][2]*ABC[2]+abc[0][2]*(ABC[0]+16*d);
			Parallelepiped( abc, tr, 5*d, d, d);
			tr[0]=Tr[0]+abc[1][0]*ABC[1]+abc[2][0]*ABC[2]+abc[0][0]*(-10*d);
			tr[1]=Tr[1]+abc[1][1]*ABC[1]+abc[2][1]*ABC[2]+abc[0][1]*(-10*d);
			tr[2]=Tr[2]+abc[1][2]*ABC[1]+abc[2][2]*ABC[2]+abc[0][2]*(-10*d);
			Parallelepiped( abc, tr, 5*d, d, d);
			tr[0]=Tr[0]+abc[1][0]*ABC[1]+abc[2][0]*ABC[2]+abc[0][0]*(-20*d);
			tr[1]=Tr[1]+abc[1][1]*ABC[1]+abc[2][1]*ABC[2]+abc[0][1]*(-20*d);
			tr[2]=Tr[2]+abc[1][2]*ABC[1]+abc[2][2]*ABC[2]+abc[0][2]*(-20*d);
			Parallelepiped( abc, tr, 5*d, d, d);	
		glEnd( );//GL_TRIANGLES
	glEndList( );

	BoundaryListB = glGenLists( 1 );
	glNewList( BoundaryListB, GL_COMPILE );
		glBegin( GL_TRIANGLES );
			glColor3f(0.7,0.7,0.7);
			tr[0]=Tr[0]+abc[1][0]*(ABC[1]+6*d); 
			tr[1]=Tr[1]+abc[1][1]*(ABC[1]+6*d);
			tr[2]=Tr[2]+abc[1][2]*(ABC[1]+6*d);
			Parallelepiped( abc, tr, d, 5*d, d);
			tr[0]=Tr[0]+abc[1][0]*(ABC[1]+16*d); 
			tr[1]=Tr[1]+abc[1][1]*(ABC[1]+16*d);
			tr[2]=Tr[2]+abc[1][2]*(ABC[1]+16*d);
			Parallelepiped( abc, tr, d, 5*d, d);
			tr[0]=Tr[0]+abc[1][0]*(-10*d); 
			tr[1]=Tr[1]+abc[1][1]*(-10*d);
			tr[2]=Tr[2]+abc[1][2]*(-10*d);
			Parallelepiped( abc, tr, d, 5*d, d);
			tr[0]=Tr[0]+abc[1][0]*(-20*d); 
			tr[1]=Tr[1]+abc[1][1]*(-20*d);
			tr[2]=Tr[2]+abc[1][2]*(-20*d);
			Parallelepiped( abc, tr, d, 5*d, d);

			tr[0]=Tr[0]+abc[0][0]*ABC[0]+abc[1][0]*(ABC[1]+6*d);
			tr[1]=Tr[1]+abc[0][1]*ABC[0]+abc[1][1]*(ABC[1]+6*d);
			tr[2]=Tr[2]+abc[0][2]*ABC[0]+abc[1][2]*(ABC[1]+6*d);
			Parallelepiped( abc, tr, d, 5*d, d);
			tr[0]=Tr[0]+abc[0][0]*ABC[0]+abc[1][0]*(ABC[1]+16*d);
			tr[1]=Tr[1]+abc[0][1]*ABC[0]+abc[1][1]*(ABC[1]+16*d);
			tr[2]=Tr[2]+abc[0][2]*ABC[0]+abc[1][2]*(ABC[1]+16*d);
			Parallelepiped( abc, tr, d, 5*d, d);
			tr[0]=Tr[0]+abc[0][0]*ABC[0]+abc[1][0]*(-10*d);
			tr[1]=Tr[1]+abc[0][1]*ABC[0]+abc[1][1]*(-10*d);
			tr[2]=Tr[2]+abc[0][2]*ABC[0]+abc[1][2]*(-10*d);
			Parallelepiped( abc, tr, d, 5*d, d);
			tr[0]=Tr[0]+abc[0][0]*ABC[0]+abc[1][0]*(-20*d);
			tr[1]=Tr[1]+abc[0][1]*ABC[0]+abc[1][1]*(-20*d);
			tr[2]=Tr[2]+abc[0][2]*ABC[0]+abc[1][2]*(-20*d);
			Parallelepiped( abc, tr, d, 5*d, d);		

			tr[0]=Tr[0]+abc[2][0]*ABC[2]+abc[1][0]*(ABC[1]+6*d);
			tr[1]=Tr[1]+abc[2][1]*ABC[2]+abc[1][1]*(ABC[1]+6*d);
			tr[2]=Tr[2]+abc[2][2]*ABC[2]+abc[1][2]*(ABC[1]+6*d);
			Parallelepiped( abc, tr, d, 5*d, d);
			tr[0]=Tr[0]+abc[2][0]*ABC[2]+abc[1][0]*(ABC[1]+16*d);
			tr[1]=Tr[1]+abc[2][1]*ABC[2]+abc[1][1]*(ABC[1]+16*d);
			tr[2]=Tr[2]+abc[2][2]*ABC[2]+abc[1][2]*(ABC[1]+16*d);
			Parallelepiped( abc, tr, d, 5*d, d);
			tr[0]=Tr[0]+abc[2][0]*ABC[2]+abc[1][0]*(-10*d);
			tr[1]=Tr[1]+abc[2][1]*ABC[2]+abc[1][1]*(-10*d);
			tr[2]=Tr[2]+abc[2][2]*ABC[2]+abc[1][2]*(-10*d);
			Parallelepiped( abc, tr, d, 5*d, d);
			tr[0]=Tr[0]+abc[2][0]*ABC[2]+abc[1][0]*(-20*d);
			tr[1]=Tr[1]+abc[2][1]*ABC[2]+abc[1][1]*(-20*d);
			tr[2]=Tr[2]+abc[2][2]*ABC[2]+abc[1][2]*(-20*d);
			Parallelepiped( abc, tr, d, 5*d, d);

			tr[0]=Tr[0]+abc[0][0]*ABC[0]+abc[2][0]*ABC[2]+abc[1][0]*(ABC[1]+6*d);
			tr[1]=Tr[1]+abc[0][1]*ABC[0]+abc[2][1]*ABC[2]+abc[1][1]*(ABC[1]+6*d);
			tr[2]=Tr[2]+abc[0][2]*ABC[0]+abc[2][2]*ABC[2]+abc[1][2]*(ABC[1]+6*d);
			Parallelepiped( abc, tr, d, 5*d, d);
			tr[0]=Tr[0]+abc[0][0]*ABC[0]+abc[2][0]*ABC[2]+abc[1][0]*(ABC[1]+16*d);
			tr[1]=Tr[1]+abc[0][1]*ABC[0]+abc[2][1]*ABC[2]+abc[1][1]*(ABC[1]+16*d);
			tr[2]=Tr[2]+abc[0][2]*ABC[0]+abc[2][2]*ABC[2]+abc[1][2]*(ABC[1]+16*d);
			Parallelepiped( abc, tr, d, 5*d, d);
			tr[0]=Tr[0]+abc[0][0]*ABC[0]+abc[2][0]*ABC[2]+abc[1][0]*(-10*d);
			tr[1]=Tr[1]+abc[0][1]*ABC[0]+abc[2][1]*ABC[2]+abc[1][1]*(-10*d);
			tr[2]=Tr[2]+abc[0][2]*ABC[0]+abc[2][2]*ABC[2]+abc[1][2]*(-10*d);
			Parallelepiped( abc, tr, d, 5*d, d);
			tr[0]=Tr[0]+abc[0][0]*ABC[0]+abc[2][0]*ABC[2]+abc[1][0]*(-20*d);
			tr[1]=Tr[1]+abc[0][1]*ABC[0]+abc[2][1]*ABC[2]+abc[1][1]*(-20*d);
			tr[2]=Tr[2]+abc[0][2]*ABC[0]+abc[2][2]*ABC[2]+abc[1][2]*(-20*d);
			Parallelepiped( abc, tr, d, 5*d, d);
		glEnd( );//GL_TRIANGLES
	glEndList( );

	BoundaryListC = glGenLists( 1 );
	glNewList( BoundaryListC, GL_COMPILE );
		glBegin( GL_TRIANGLES );
			glColor3f(0.7,0.7,0.7);
			tr[0]=Tr[0]+abc[2][0]*(ABC[2]+6*d); 
			tr[1]=Tr[1]+abc[2][1]*(ABC[2]+6*d);
			tr[2]=Tr[2]+abc[2][2]*(ABC[2]+6*d);
			Parallelepiped( abc, tr, d, d, 5*d);
			tr[0]=Tr[0]+abc[2][0]*(ABC[2]+16*d); 
			tr[1]=Tr[1]+abc[2][1]*(ABC[2]+16*d);
			tr[2]=Tr[2]+abc[2][2]*(ABC[2]+16*d);
			Parallelepiped( abc, tr, d, d, 5*d);
			tr[0]=Tr[0]+abc[2][0]*(-10*d); 
			tr[1]=Tr[1]+abc[2][1]*(-10*d);
			tr[2]=Tr[2]+abc[2][2]*(-10*d);
			Parallelepiped( abc, tr, d, d, 5*d);
			tr[0]=Tr[0]+abc[2][0]*(-20*d); 
			tr[1]=Tr[1]+abc[2][1]*(-20*d);
			tr[2]=Tr[2]+abc[2][2]*(-20*d);
			Parallelepiped( abc, tr, d, d, 5*d);

			tr[0]=Tr[0]+abc[0][0]*ABC[0]+abc[2][0]*(ABC[2]+6*d);
			tr[1]=Tr[1]+abc[0][1]*ABC[0]+abc[2][1]*(ABC[2]+6*d);
			tr[2]=Tr[2]+abc[0][2]*ABC[0]+abc[2][2]*(ABC[2]+6*d);
			Parallelepiped( abc, tr, d, d, 5*d);
			tr[0]=Tr[0]+abc[0][0]*ABC[0]+abc[2][0]*(ABC[2]+16*d);
			tr[1]=Tr[1]+abc[0][1]*ABC[0]+abc[2][1]*(ABC[2]+16*d);
			tr[2]=Tr[2]+abc[0][2]*ABC[0]+abc[2][2]*(ABC[2]+16*d);
			Parallelepiped( abc, tr, d, d, 5*d);
			tr[0]=Tr[0]+abc[0][0]*ABC[0]+abc[2][0]*(-10*d);
			tr[1]=Tr[1]+abc[0][1]*ABC[0]+abc[2][1]*(-10*d);
			tr[2]=Tr[2]+abc[0][2]*ABC[0]+abc[2][2]*(-10*d);
			Parallelepiped( abc, tr, d, d, 5*d);
			tr[0]=Tr[0]+abc[0][0]*ABC[0]+abc[2][0]*(-20*d);
			tr[1]=Tr[1]+abc[0][1]*ABC[0]+abc[2][1]*(-20*d);
			tr[2]=Tr[2]+abc[0][2]*ABC[0]+abc[2][2]*(-20*d);
			Parallelepiped( abc, tr, d, d, 5*d);		

			tr[0]=Tr[0]+abc[1][0]*ABC[1]+abc[2][0]*(ABC[2]+6*d);
			tr[1]=Tr[1]+abc[1][1]*ABC[1]+abc[2][1]*(ABC[2]+6*d);
			tr[2]=Tr[2]+abc[1][2]*ABC[1]+abc[2][2]*(ABC[2]+6*d);
			Parallelepiped( abc, tr, d, d, 5*d);
			tr[0]=Tr[0]+abc[1][0]*ABC[1]+abc[2][0]*(ABC[2]+16*d);
			tr[1]=Tr[1]+abc[1][1]*ABC[1]+abc[2][1]*(ABC[2]+16*d);
			tr[2]=Tr[2]+abc[1][2]*ABC[1]+abc[2][2]*(ABC[2]+16*d);
			Parallelepiped( abc, tr, d, d, 5*d);
			tr[0]=Tr[0]+abc[1][0]*ABC[1]+abc[2][0]*(-10*d);
			tr[1]=Tr[1]+abc[1][1]*ABC[1]+abc[2][1]*(-10*d);
			tr[2]=Tr[2]+abc[1][2]*ABC[1]+abc[2][2]*(-10*d);
			Parallelepiped( abc, tr, d, d, 5*d);
			tr[0]=Tr[0]+abc[1][0]*ABC[1]+abc[2][0]*(-20*d);
			tr[1]=Tr[1]+abc[1][1]*ABC[1]+abc[2][1]*(-20*d);
			tr[2]=Tr[2]+abc[1][2]*ABC[1]+abc[2][2]*(-20*d);
			Parallelepiped( abc, tr, d, d, 5*d);

			tr[0]=Tr[0]+abc[0][0]*ABC[0]+abc[1][0]*ABC[1]+abc[2][0]*(ABC[2]+6*d);
			tr[1]=Tr[1]+abc[0][1]*ABC[0]+abc[1][1]*ABC[1]+abc[2][1]*(ABC[2]+6*d);
			tr[2]=Tr[2]+abc[0][2]*ABC[0]+abc[1][2]*ABC[1]+abc[2][2]*(ABC[2]+6*d);
			Parallelepiped( abc, tr, d, d, 5*d);
			tr[0]=Tr[0]+abc[0][0]*ABC[0]+abc[1][0]*ABC[1]+abc[2][0]*(ABC[2]+16*d);
			tr[1]=Tr[1]+abc[0][1]*ABC[0]+abc[1][1]*ABC[1]+abc[2][1]*(ABC[2]+16*d);
			tr[2]=Tr[2]+abc[0][2]*ABC[0]+abc[1][2]*ABC[1]+abc[2][2]*(ABC[2]+16*d);
			Parallelepiped( abc, tr, d, d, 5*d);
			tr[0]=Tr[0]+abc[0][0]*ABC[0]+abc[1][0]*ABC[1]+abc[2][0]*(-10*d);
			tr[1]=Tr[1]+abc[0][1]*ABC[0]+abc[1][1]*ABC[1]+abc[2][1]*(-10*d);
			tr[2]=Tr[2]+abc[0][2]*ABC[0]+abc[1][2]*ABC[1]+abc[2][2]*(-10*d);
			Parallelepiped( abc, tr, d, d, 5*d);
			tr[0]=Tr[0]+abc[0][0]*ABC[0]+abc[1][0]*ABC[1]+abc[2][0]*(-20*d);
			tr[1]=Tr[1]+abc[0][1]*ABC[0]+abc[1][1]*ABC[1]+abc[2][1]*(-20*d);
			tr[2]=Tr[2]+abc[0][2]*ABC[0]+abc[1][2]*ABC[1]+abc[2][2]*(-20*d);
			Parallelepiped( abc, tr, d, d, 5*d);
		glEnd( );//GL_TRIANGLES
	glEndList( );

	// create the axes:
	float xyz[][3]={ {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0} };
	tr[0] = -d/2;
	tr[1] = -d/2;
	tr[2] = -d/2;
	AxesList = glGenLists( 1 );
	glNewList( AxesList, GL_COMPILE );
		glBegin( GL_TRIANGLES );
			glColor3f(1,0,0);
			Parallelepiped( xyz, tr, 2, d, d);//(0,0,0)-->(1,0,0)
			glColor3f(0,1,0);
			Parallelepiped( xyz, tr, d, 2, d);//(0,0,0)-->(0,1,0)
			glColor3f(0,0,1);
			Parallelepiped( xyz, tr, d, d, 2);//(0,0,0)-->(0,0,1)
			tr[0]+= -d/10;
			tr[1]+= -d/10;
			tr[2]+= -d/10;
			glColor3f(0.2,0.2,0.2);		
			Parallelepiped( xyz, tr, d+d/5, d+d/5, d+d/5);//(0,0,0)-->(1,0,0)
		glEnd();
	glEndList( );
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

	drawVBO(VCNum, WhichVectorMode); // Draw VBO for spins
	
	// possibly draw the box and periodic boundary condition :
	if( BoxOn != 0 ) 
		{	glCallList( BoxList);
			if(Boundary[0]!=0) glCallList( BoundaryListA );
			if(Boundary[1]!=0) glCallList( BoundaryListB );
			if(Boundary[2]!=0) glCallList( BoundaryListC );
		}
	// possibly draw the axes:
	if( AxesOn != 0 ) glCallList( AxesList );
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
	#ifndef __APPLE__
	glewInit();
	#endif

// initialize element of the drawing scene which remain unchanged e.g. coordinate basis, domain boundary etc.
	InitLists(abc, ABC);
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
	glCullFace(GL_FRONT);
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

	if(timeInterval > 40)//40ms gives approximately 25 FPS +/-1 if the engine works faster then 25 IPS
	{
		if( FLAG_SHOW==TAKE_DATA )
		{
			if (DataTransfer==1)
			{
			UpdateVerticesNormalsColors(vertexProto, normalProto, VCNumProto, vertices, normals, colors, VCNum, Px, Py, Pz, bSx, bSy, bSz, WhichVectorMode);
			UpdateVBO( &vboIdV, &vboIdN, &vboIdC, &iboIdI, vertices, normals, colors, indices);	
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
			if (timeInterval > 2000)//~2 seconds
			{
				FPS = frameCount / (timeInterval * 0.002f);
				previousTime = currentTime;
				frameCount = 0;	
				IPS = (currentIteration - previousIteration)/ (timeInterval * 0.002f);;
				previousIteration = currentIteration;
			}
			EnterCriticalSection( &show_mutex);
				FLAG_SHOW = READY; // meaning that OpenGL is ready to take new data from engine
			LeaveCriticalSection( &show_mutex);
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
		EnterCriticalSection(&culc_mutex);
			FLAG_CALC=DO_IT;
		LeaveCriticalSection(&culc_mutex);
	}else{
		Play=0; 
		TwDefine(" Parameters&Controls/Run  label='RUN simulation' ");
		EnterCriticalSection(&culc_mutex);
			FLAG_CALC=WAIT;
		LeaveCriticalSection(&culc_mutex);	
	}
}

void TW_CALL CB_SetScale(const void *value, void *clientData )
{
	(void)clientData; // unused
    Scale = *( float *)value; // copy value to Scale
if(WhichVectorMode==ARROW1||WhichVectorMode==CONE1) 
   {UpdatePrototypeVerNorInd(vertexProto, normalProto, indicesProto, arrowFaces, WhichVectorMode);}
	UpdateVerticesNormalsColors(vertexProto, normalProto, VCNumProto, vertices, normals, colors, VCNum, Px, Py, Pz, Sx, Sy, Sz, WhichVectorMode);
	UpdateVBO(&vboIdV, &vboIdN, &vboIdC, &iboIdI,vertices, normals, colors, indices);
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
    WhichVectorMode = *( enVectorMode *)value; // copy value to Scale
		ReallocateArrayDrawing(arrowFaces, WhichVectorMode);
		UpdatePrototypeVerNorInd(vertexProto, normalProto, indicesProto, arrowFaces, WhichVectorMode);
		UpdateIndices(indicesProto , IdNumProto, indices, IdNum, VCNumProto); 
		UpdateVerticesNormalsColors(vertexProto, normalProto, VCNumProto, vertices, normals, colors, VCNum, Px, Py, Pz, Sx, Sy, Sz, WhichVectorMode);
		UpdateVBO(&vboIdV, &vboIdN, &vboIdC, &iboIdI,vertices, normals, colors, indices);
}

//  Callback function called by the tweak bar to get the 'AutoRotate' value
void TW_CALL CB_GetVectorMode(void *value, void *clientData)
{
    (void)clientData; // unused
    *(int *)value = WhichVectorMode; // just copy Scale to value
}

void TW_CALL CB_SetSliceMode(const void *value, void *clientData )
{
	(void)clientData; // unused
    WhichSliceMode = *( enSliceMode *)value; // copy value to Scale
		ReallocateArrayDrawing(arrowFaces, WhichVectorMode);
		UpdatePrototypeVerNorInd(vertexProto, normalProto, indicesProto, arrowFaces, WhichVectorMode);
		UpdateIndices(indicesProto , IdNumProto, indices, IdNum, VCNumProto); 
		UpdateVerticesNormalsColors(vertexProto, normalProto, VCNumProto, vertices, normals, colors, VCNum, Px, Py, Pz, Sx, Sy, Sz, WhichVectorMode);
		UpdateVBO(&vboIdV, &vboIdN, &vboIdC, &iboIdI,vertices, normals, colors, indices);
}

//  Callback function called by the tweak bar to get the 'AutoRotate' value
void TW_CALL CB_GetSliceMode(void *value, void *clientData)
{
    (void)clientData; // unused
    *(int *)value = WhichSliceMode; // just copy Scale to value
}

void TW_CALL CB_SetFaces(const void *value, void *clientData )
{
	(void)clientData; // unused
    arrowFaces = *( int *)value; // copy value to Scale
		ReallocateArrayDrawing(arrowFaces, WhichVectorMode);
		UpdatePrototypeVerNorInd(vertexProto, normalProto, indicesProto, arrowFaces, WhichVectorMode);
		UpdateIndices(indicesProto , IdNumProto, indices, IdNum, VCNumProto); 
		UpdateVerticesNormalsColors(vertexProto, normalProto, VCNumProto, vertices, normals, colors, VCNum, Px, Py, Pz, Sx, Sy, Sz, WhichVectorMode);
		UpdateVBO(&vboIdV, &vboIdN, &vboIdC, &iboIdI,vertices, normals, colors, indices);
}

//  Callback function called by the tweak bar to get the 'AutoRotate' value
void TW_CALL CB_GetFaces(void *value, void *clientData)
{
    (void)clientData; // unused
    *(int *)value = arrowFaces; // just copy Scale to value
}

void TW_CALL CB_SetPivot(const void *value, void *clientData )
{
	(void)clientData; // unused
    Pivot = *( float *)value; // copy value to Scale
		ReallocateArrayDrawing(arrowFaces, WhichVectorMode);
		UpdatePrototypeVerNorInd(vertexProto, normalProto, indicesProto, arrowFaces, WhichVectorMode);
		UpdateIndices(indicesProto , IdNumProto, indices, IdNum, VCNumProto); 
		UpdateVerticesNormalsColors(vertexProto, normalProto, VCNumProto, vertices, normals, colors, VCNum, Px, Py, Pz, Sx, Sy, Sz, WhichVectorMode);
		UpdateVBO(&vboIdV, &vboIdN, &vboIdC, &iboIdI,vertices, normals, colors, indices);
}

//  Callback function called by the tweak bar to get the 'AutoRotate' value
void TW_CALL CB_GetPivot(void *value, void *clientData)
{
    (void)clientData; // unused
    *(float *)value = Pivot; // just copy Scale to value
}


void TW_CALL CB_SetColorShift(const void *value, void *clientData )
{
	(void)clientData; // unused
    ColorShift = *( int *)value; // copy value to Scale
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
		UpdateVerticesNormalsColors(vertexProto, normalProto, VCNumProto, vertices, normals, colors, VCNum, Px, Py, Pz, Sx, Sy, Sz, WhichVectorMode);
		UpdateVBO(&vboIdV, &vboIdN, &vboIdC, &iboIdI, vertices, normals, colors, indices);
}

void TW_CALL CB_GetColorShift(void *value, void *clientData)
{
    (void)clientData; // unused
    *(int *)value = ColorShift; // just copy Scale to value
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
		UpdateVerticesNormalsColors(vertexProto, normalProto, VCNumProto, vertices, normals, colors, VCNum, Px, Py, Pz, Sx, Sy, Sz, WhichVectorMode);
		UpdateVBO(&vboIdV, &vboIdN, &vboIdC, &iboIdI, vertices, normals, colors, indices);
}

void TW_CALL CB_GetColorScheme(void *value, void *clientData)
{
    (void)clientData; // unused
    *(enColorScheme *)value = WhichColorScheme; // just copy Scale to value
}


void TW_CALL CB_SetInvHue(const void *value, void *clientData )
{
	(void)clientData; // unused
    InvertHue = *( int *)value; // copy value to Scale
		UpdateVerticesNormalsColors(vertexProto, normalProto, VCNumProto, vertices, normals, colors, VCNum, Px, Py, Pz, Sx, Sy, Sz, WhichVectorMode);
		UpdateVBO(&vboIdV, &vboIdN, &vboIdC, &iboIdI, vertices, normals, colors, indices);
}

void TW_CALL CB_GetInvHue(void *value, void *clientData)
{
    (void)clientData; // unused
    *(int *)value = InvertHue; // just copy Scale to value
}

void TW_CALL CB_SetInvVal(const void *value, void *clientData )
{
	(void)clientData; // unused
    InvertValue = *( int *)value; // copy value to Scale
		UpdateVerticesNormalsColors(vertexProto, normalProto, VCNumProto, vertices, normals, colors, VCNum, Px, Py, Pz, Sx, Sy, Sz, WhichVectorMode);
		UpdateVBO(&vboIdV, &vboIdN, &vboIdC, &iboIdI, vertices, normals, colors, indices);
}

void TW_CALL CB_GetInvVal(void *value, void *clientData)
{
    (void)clientData; // unused
    *(int *)value = InvertValue; // just copy InvertValue to value
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
    *(float *)value = Period_dc; // just copy Scale to value
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
    *(float *)value = Omega_dc; // just copy Scale to value
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
			Sx[i] = tmp[0];
			Sy[i] = tmp[1];
			Sz[i] = tmp[2];
		}
		ChangeVectorMode(1);		
	}
}

void TW_CALL CB_InvertX( void *clientData )
{
	for (int i=0; i<NOS; i++) {Sx[i] = -Sx[i];}
	ChangeVectorMode(1);
}

void TW_CALL CB_InvertY( void *clientData )
{
	for (int i=0; i<NOS; i++) {Sy[i] = -Sy[i];}
	ChangeVectorMode(1);
}

void TW_CALL CB_InvertZ( void *clientData )
{
	for (int i=0; i<NOS; i++) {Sz[i] = -Sz[i];}
	ChangeVectorMode(1);
}

void TW_CALL CB_CleanSxSySzFile( void *clientData )
{
  fclose (outFile);//outFile is a global variable - pointer FILE* see also CALC_THREAD in ENGINE.cpp
 	outFile = fopen ("sxsysz.csv","w");
	if (outFile!=NULL) {fputs ("iter,sx,sy,sz,e_tot,\n",outFile);}
 }


void TW_CALL CB_ResetIterations( void *clientData )
{
  ITERATION=0;
 }


void TW_CALL CB_SaveCSV( void *clientData )
{
FILE * pFile;
  pFile = fopen (outputfilename,"w");
  if (pFile!=NULL)
  {
  	fputs ("px,py,pz,nx,ny,nz,\n",pFile);
  	for (int i=0;i<NOS;i++)
  	{
	snprintf(shortBufer,80,"%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,%2.5f,\n",Px[i],Py[i],Pz[i],Sx[i],Sy[i],Sz[i]);
    fputs (shortBufer,pFile);  		
  	}
    fclose (pFile);
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
				Sx[i]=sx;
				Sy[i]=sy;
				Sz[i]=sz;
		 	}
		 	i++;
		}while(c != EOF); 
    fclose(pFile);           
    }
		UpdateVerticesNormalsColors(vertexProto, normalProto, VCNumProto, vertices, normals, colors, VCNum, Px, Py, Pz, Sx, Sy, Sz, WhichVectorMode);
		UpdateVBO(&vboIdV, &vboIdN, &vboIdC, &iboIdI, vertices, normals, colors, indices);
}

void setupTweakBar()
{
	TwDefine(" GLOBAL iconpos=topleft "); // icons go to top-left corner of the window
	TwDefine(" GLOBAL iconalign=horizontal "); // icons will be aligned horizontally
	TwDefine(" GLOBAL contained=true "); // bars cannot move outside of the window
	TwDefine(" GLOBAL iconmargin='1 1' ");
////////////////////////////////////////////// Help Bar F1 /////////////////////////////////////////////////////////////
	help_bar = TwGetBarByIndex(0);
	TwDefine(" TW_HELP size='440 510' color='70 100 100'");
	TwDefine(" TW_HELP help='F1: show/hide (this) Help bar' "); // change default tweak bar size and color
//////////////////////////////////////////////// View F2 ///////////////////////////////////////////////////////////////
    view_bar = TwNewBar("View");
    TwDefine(" View iconified=true "); 
    TwDefine(" View size='220 510' color='224 216 96' alpha=200 "); // change default tweak bar size and color
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
	TwEnumVal		enVectorModeTw[] = { {ARROW1, "Arrows"}, {CONE1, "Cones"}, {CANE, "Canes"}, {POINT, "Points"} };
	TwType			TV_TYPE_VEC_MOD = TwDefineEnum("Type_of_vectors", enVectorModeTw, 4);
	TwAddVarCB(view_bar, "Type of vectors", TV_TYPE_VEC_MOD, CB_SetVectorMode, CB_GetVectorMode, &WhichVectorMode, "help='Type of 3D vectors' ");
	}

	{
	TwEnumVal		enSliceModeTw[] = { {A_AXIS, "a-axis"}, {B_AXIS, "b-axis"}, {C_AXIS, "c-axis"}, {ABC_AXIS, "a+b+c"}};
	TwType			TV_TYPE_VEC_MOD = TwDefineEnum("Slicing", enSliceModeTw, 4);
	TwAddVarCB(view_bar, "Slicing plane", TV_TYPE_VEC_MOD, CB_SetSliceMode, CB_GetSliceMode, &WhichSliceMode, "help='Slising plane perpenticulat to the choosen axis' ");
	}

	TwAddVarCB(view_bar, "Pivot", TW_TYPE_FLOAT, CB_SetPivot, CB_GetPivot,  &Pivot, " min=0 max=1 step=0.01 help='Pivot of 3D arrow.' ");
	TwAddVarCB(view_bar, "Faces", TW_TYPE_INT32, CB_SetFaces, CB_GetFaces,  &arrowFaces, " min=3 max=20 step=1 help='Number of faces for 3D arrow.' ");
	TwAddVarCB(view_bar, "Scale", TW_TYPE_FLOAT, CB_SetScale, CB_GetScale,  &Scale, " min=0.1 max=2.5 step=0.01 keyIncr='+' keyDecr='-' help='Scale the vectors.' ");

	TwAddVarRW(view_bar, "Show basis", TW_TYPE_BOOL32, &AxesOn, " key=ALT+b ");
	TwAddVarRW(view_bar, "Show box", TW_TYPE_BOOL32, &BoxOn, " key=CTRL+b ");

	TwAddVarRW(view_bar, "Intensity", TW_TYPE_FLOAT, &g_LightMultiplier, 
	" label='Light intensity' min=0.1 max=4 step=0.02 help='Increase/decrease the light power.' group='Light' ");
	TwAddVarRW(view_bar, "LightDir", TW_TYPE_DIR3F, &g_LightDirection, 
	" label='Light direction' opened=true help='Change the light direction.' group='Light'");
	temp_color[0] = 230;
	temp_color[1] = 230;
	temp_color[2] = 255;
	TwSetParam(view_bar, "LightDir", "arrowcolor", TW_PARAM_INT32, 3, temp_color);

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

///////////////////////////////////////// Parameters&Controls F3 ///////////////////////////////////////////////////////
	control_bar = TwNewBar("Parameters&Controls");
    TwDefine(" Parameters&Controls iconified=true "); 
    TwDefine(" Parameters&Controls size='220 510' color='224 96 216' alpha=200 "); // change default tweak bar size and color
    TwDefine(" Parameters&Controls help='F3: show/hide Control bar' "); // change default tweak bar size and color

	TwAddButton(control_bar, "Run", CB_Run, NULL, "key=SPACE label='RUN simulation' ");
	TwAddVarRW(control_bar, "Record", TW_TYPE_BOOL32, &Record, 
	"keyIncr='r' label='Recording' true='Rec.' false='Pause' help='Recording <sx>, <sy>, <sz> in each iteration'");
	TwAddVarRW(control_bar, "Rec_Iteration", TW_TYPE_INT32, &rec_iteration, 
	"label='Every iteration' min=1 max=1000 step=1 ");
	TwAddButton(control_bar, "Clean the record", CB_CleanSxSySzFile, NULL, "label= 'Clean the record' help='clean the output file with <sx>, <sy>, <sz>' ");
	TwAddButton(control_bar, "Reset iterations", CB_ResetIterations, NULL, "label='Reset iterations' ");

	TwAddSeparator(control_bar, "sep", NULL);

	TwAddVarRW(control_bar, "BCinA", TW_TYPE_BOOL32, &Boundary[0], 
	"label='along a' group='Boundary conditions' true='periodic' false='open' help='set boundary conditions along translation vector 'a' '");
	TwAddVarRW(control_bar, "BCinB", TW_TYPE_BOOL32, &Boundary[1], 
	"label='along b' group='Boundary conditions' true='periodic' false='open' help='set boundary conditions along translation vector 'b' '");
	TwAddVarRW(control_bar, "BCinC", TW_TYPE_BOOL32, &Boundary[2], 
	"label='along c' group='Boundary conditions' true='periodic' false='open' help='set boundary conditions along translation vector 'c' '");
    
	TwAddVarRW(control_bar, "Damping", TW_TYPE_FLOAT, &damping, 
	"label='Damping' min=0 max=10 group='LLG' ");
	TwAddVarRW(control_bar, "Time_step", TW_TYPE_FLOAT, &t_step, 
	"label='Time step' min=0 max=0.1 group='LLG' ");
	TwAddVarRW(control_bar, "temperature", TW_TYPE_FLOAT, &Temperature, 
	"label='k_b*T' min=0 max=100 step=0.01 group='LLG' ");

	TwAddVarRW(control_bar, "FieldDir", TW_TYPE_DIR3F, &VHf, 
	"label='Field direction' opened=true help='Change the direction of applied field' ");
	temp_color[0] = 55;
	temp_color[1] = 55;
	temp_color[2] = 155;
	TwSetParam(control_bar, "FieldDir", "arrowcolor", TW_PARAM_INT32, 3, temp_color);
	TwAddVarRW(control_bar, "Field", TW_TYPE_FLOAT, &Hf, 
	"label='Field' help='The value of uniaxial anisotropy' ");

	TwAddVarRW(control_bar, "KudDir", TW_TYPE_DIR3F, &VKu, 
	"label='Ku' opened=true help='The direction of uniaxial anisotropy' ");
	temp_color[0] = 55;
	temp_color[1] = 155;
	temp_color[2] = 55;
	TwSetParam(control_bar, "KudDir", "arrowcolor", TW_PARAM_INT32, 3, temp_color);
	TwAddVarRW(control_bar, "Ku", TW_TYPE_FLOAT, &Ku, 
	"label='Ku' help='The value of uniaxial anisotropy' ");
TwAddSeparator(control_bar, "sep0", NULL);

	TwAddVarRW(control_bar, "Kcub", TW_TYPE_FLOAT, &Kc, 
	"label='Kc' help='The value of cubic anisotropy' ");

TwAddSeparator(control_bar, "sep1", NULL);

	for(int s=0; s<ShellNumber; s++)
	{
	snprintf(shortBufer,50,"J%1i",s);
	TwAddVarRW(control_bar, shortBufer, TW_TYPE_FLOAT, &Jij[s], 
	"help='Heisenberg exchange' ");
	TwSetParam(control_bar, shortBufer, "label",  TW_PARAM_CSTRING , 1, shortBufer);
	}

TwAddSeparator(control_bar, "sep2", NULL);
	
	for(int s=0; s<ShellNumber; s++)
	{
	snprintf(shortBufer,50,"B%1i",s);
	TwAddVarRW(control_bar, shortBufer, TW_TYPE_FLOAT, &Bij[s], 
	"help='Bi-quadratic exchange' ");
	TwSetParam(control_bar, shortBufer, "label",  TW_PARAM_CSTRING , 1, shortBufer);
	}

TwAddSeparator(control_bar, "sep3", NULL);

	for(int s=0; s<ShellNumber; s++)
	{
	snprintf(shortBufer,50,"D%1i",s);
	TwAddVarRW(control_bar, shortBufer, TW_TYPE_FLOAT, &Dij[s], 
	"help='Dzyaloshinskii-Moriya' ");
	TwSetParam(control_bar, shortBufer, "label",  TW_PARAM_CSTRING , 1, shortBufer);
	}

TwAddSeparator(control_bar, "sep4", NULL);
	TwAddVarRW(control_bar, "ST param.", TW_TYPE_FLOAT, &Cu, 
	"help='Dzyaloshinskii-Moriya' ");
	TwAddVarRW(control_bar, "CurrentDir", TW_TYPE_DIR3F, &VCu, 
	"label='cur. dir.' opened=true help='The polarization direction of electric current' ");

///////////////////////////////////////// Initial state F4 ///////////////////////////////////////////////////////
	initial_bar = TwNewBar("Initial_State");
	TwDefine(" Initial_State iconified=true "); 
	TwDefine(" Initial_State size='220 510' color='180 180 254'  alpha=200"); // change default tweak bar size and color
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

	TwAddButton(initial_bar, "Set initial", CB_SetInitial, NULL, "key=INSERT label='insert initial state' ");

	TwAddSeparator(initial_bar, "sep01", NULL);
	TwAddVarRW(initial_bar, "Degrees", TW_TYPE_FLOAT,  &RotateAllSpins, 
	" min=-360 max=360 step=1 help='Rotate all spins about characteristic direction' ");
	TwAddButton(initial_bar, "Rotate spins", CB_RotateAllSpins, NULL, "label='rotate all spins' ");
	TwAddSeparator(initial_bar, "sep02", NULL);

	TwAddSeparator(initial_bar, "sep1", NULL);
	TwAddButton(initial_bar, "Invert X", CB_InvertX, NULL, "label='invert n_x component' ");
	TwAddButton(initial_bar, "Invert Y", CB_InvertY, NULL, "label='invert n_y component' ");
	TwAddButton(initial_bar, "Invert Z", CB_InvertZ, NULL, "label='invert n_z component' ");
	TwAddSeparator(initial_bar, "sep2", NULL);

	TwAddVarRW(initial_bar, "Output file name:", TW_TYPE_CSSTRING(sizeof(outputfilename)), outputfilename, ""); 
	TwAddButton(initial_bar, "Write to CSV", CB_SaveCSV, NULL, "label='save state as comma separated values' ");
	TwAddVarRW(initial_bar, "Input file name:", TW_TYPE_CSSTRING(sizeof(inputfilename)), inputfilename, "");
	TwAddButton(initial_bar, "Read from CSV", CB_ReadCSV, NULL, "label='read state from JSINI.csv' ");	

///////////////////////////////////////// AC field F5 ///////////////////////////////////////////////////////
	ac_field_bar = TwNewBar("AC_Field");
	TwDefine(" AC_Field iconified=true "); 
	TwDefine(" AC_Field size='300 510' color='180 100 140'  alpha=200"); // change default tweak bar size and color
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


///////////////////////////////////////// Info bar F12 ///////////////////////////////////////////////////////
	info_bar = TwNewBar("Info");
	TwDefine(" Info refresh=0.5 ");
	TwDefine(" Info iconified = false movable = false alwaysbottom=true resizable=false fontstyle=fixed fontsize=2"); 
	TwDefine(" Info help='F12: show/hide info-bar' "); // change default tweak bar size and color
	TwDefine(" Info color='10 10 10' alpha=0 "); // change default tweak bar size and color
	TwDefine(" Info help='F12: show/hide info-bar' "); // change default tweak bar size and color
	TwDefine(" Info position = '1 30' size ='170 400' valueswidth=120"); // change default tweak bar size and color
	TwAddVarRO(info_bar, "R/S", TW_TYPE_BOOL32,  &Play, "true='RUNING' false='STOPED' ");
	TwAddVarRO(info_bar, "Rec.", TW_TYPE_BOOL32,  &Record, "true='On' false='Off' ");
	TwAddVarRO(info_bar, "ACF", TW_TYPE_BOOL32,  &AC_FIELD_ON, "true='On' false='Off' help='AC filed on/off'");
	TwAddSeparator(info_bar, "sep", NULL);
	TwAddVarRO(info_bar, "NOS", TW_TYPE_INT32,  &NOS, "help='Number of spins' ");
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
}


void
MouseButton( int button, int state, int x, int y ) // called when the mouse button transitions down or up
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


void 
MouseMotion( int x, int y ) // called when the mouse moves while a button is down
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
void
Keyboard( unsigned char c, int x, int y )
{
if( !TwEventKeyboardGLUT(c, x, y) )  // send event to AntTweakBar 
  { // event has not been handled by AntTweakBar
    // your code here to handle the event	
	if( false )
		fprintf( stderr, "Keyboard: '%c' (0x%0x)\n", c, c );
		switch( c )
		{
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
				Rot[0] += 0.25;
				//TransXYZ[1]+=1;
				break;
			case 's':
			case 'S':
				Rot[0] -= 0.25;
				//TransXYZ[1]-=1;
				break;
			case 'a':
			case 'A':
				Rot[1] -= 0.25;
				//TransXYZ[0]-=1;	
				break;
			case 'd':
			case 'D':
				Rot[1] += 0.25;
				//TransXYZ[0]+=1;	
				break;
			case 'o':
			case 'O':
				
				break;

			case 'v':
			case 'V':
				WhichVectorMode=ARROW1;
				ChangeVectorMode(0);
			break;

			case 'c':
			case 'C':
				WhichVectorMode=CONE1;
				ChangeVectorMode(0);
			break;

			case 'b':
			case 'B':
				WhichVectorMode=CANE;
				ChangeVectorMode(0);
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

void
KeyboardAdd( int key, int x, int y )
{
int isiconified;
if( !TwEventSpecialGLUT(key, x, y) )  // send event to AntTweakBar TwEventSpecialGLUT
    { // event has not been handled by AntTweakBar
      // your code here to handle the event	
		switch( key )
		{
			case GLUT_KEY_UP:
				switch (WhichSliceMode)
				{case A_AXIS:
					if (Alayer<ABC[0]) Alayer+=1;
				 break;
				 case B_AXIS:
				    if (Blayer<ABC[1]) Blayer+=1;
				 break;
				 case C_AXIS:
				 	if (Clayer<ABC[2]) Clayer+=1;
				 break;
				 case ABC_AXIS:
				 	//if (Clayer<ABC[2]) Clayer+=1;
				 break;				 
				}
				ChangeVectorMode(1);
				break;
			case GLUT_KEY_DOWN:
				switch (WhichSliceMode)
				{case A_AXIS:
					if (Alayer>1) Alayer-=1;
				 break;
				 case B_AXIS:
				    if (Blayer>1) Blayer-=1;
				 break;
				 case C_AXIS:
				 	if (Clayer>1) Clayer-=1;
				 break;
				 case ABC_AXIS:
				 	//if (Clayer<ABC[2]) Clayer+=1;
				 break;	
				}
				ChangeVectorMode(1);
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
			case  GLUT_KEY_F12:
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




// glui buttons callback:
void 
Buttons( int id )
{
	switch( id )
	{
		case XUP:

			Xup( );
			//glui_window->sync_live( );
			glutSetWindow( GLUT_window );
			glutPostRedisplay( );
			break;

		case YUP:

			Yup( );
			//glui_window->sync_live( );
			glutSetWindow( GLUT_window );
			glutPostRedisplay( );
			break;

		case ZUP:

			Zup( );
			//glui_window->sync_live( );
			glutSetWindow( GLUT_window );
			glutPostRedisplay( );
			break;

		case RESET:

			Reset( );
			//glui_window->sync_live( );
			glutSetWindow( GLUT_window );
			glutPostRedisplay( );
			break;

		case QUIT:

			//glui_window->close( );					// gracefully close the glui window
			glutSetWindow( GLUT_window );	//
			glFinish( );					// gracefully close out the graphics
			glutDestroyWindow( GLUT_window );// gracefully close the graphics window
			//glutLeaveMainLoop();
			exit( 0 );						// gracefully exit the program
			break;

		// case PLAY:
		// 	if (Play==0)
		// 	{
		// 		Play=1;
		// 		TwDefine(" Parameters&Controls/Run  label='STOP simulation' ");
		// 		EnterCriticalSection(&culc_mutex);
		// 			FLAG_CALC=DO_IT;
		// 		LeaveCriticalSection(&culc_mutex);
		// 	}else{
		// 		Play=0; 
		// 		TwDefine(" Parameters&Controls/Run  label='RUN simulation' ");
		// 		EnterCriticalSection(&culc_mutex);
		// 			FLAG_CALC=WAIT;
		// 		LeaveCriticalSection(&culc_mutex);	
		// 	}
		// 	//glui_window->sync_live( );
		// 	break;

		default:
			fprintf( stderr, "Don't know what to do with Button ID %d\n", id );
	}

}

void
ChangeInitialState ( int id )
{
	InitSpinComponents( Px, Py, Pz, Sx, Sy, Sz, id );
	ChangeVectorMode ( 1 );
}

void
ChangeVectorMode ( int id )
{
	switch ( id )
	{
		case 0: // change of mode e.g. from point to arrow
		ReallocateArrayDrawing(arrowFaces, WhichVectorMode);
		// Fill array for prototype (arrow or cane) array 
		UpdatePrototypeVerNorInd(vertexProto, normalProto, indicesProto, arrowFaces, WhichVectorMode);
		// Fill big array for indecies for all ARROW1 or cans 
		UpdateIndices(indicesProto , IdNumProto, indices, IdNum, VCNumProto); 
		case 1:	// change layer for drawing
		UpdateVerticesNormalsColors(vertexProto, normalProto, VCNumProto, vertices, normals, colors, VCNum, Px, Py, Pz, Sx, Sy, Sz, WhichVectorMode);
		UpdateVBO(&vboIdV, &vboIdN, &vboIdC, &iboIdI,vertices, normals, colors, indices);
	}
}


void
ChangeColorMap ( int id )
{
	switch (id)
	{
		case 0:
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
		case 1:
		UpdateVerticesNormalsColors(vertexProto, normalProto, VCNumProto, vertices, normals, colors, 
													VCNum, Px, Py, Pz, Sx, Sy, Sz, WhichVectorMode);
		UpdateVBO(&vboIdV, &vboIdN, &vboIdC, &iboIdI, vertices, normals, colors, indices);
	}
}

void
HSVtoRGB(float Vec[3], float RGB[3], int inV, int inH )
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

void
InitRGB(float* R, float* G, float* B, int *hueMap)
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



void UpdatePrototypeVerNorInd(float * V, float * N, GLuint * I, int faces, int mode)//faces = arrowFaces
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

	case POINT:
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

void  
UpdateIndices(GLuint * Iinp , int Kinp, GLuint * Iout, int Kout, int VerN)
{	// VerN number of vertesiec components per one prototype arrow
	int i,j,NOS_L=NOS_CL;
	switch( WhichSliceMode)
	{
		case A_AXIS:
		NOS_L=NOS_AL;
		break;
		case B_AXIS:
		NOS_L=NOS_BL;
		break;
		case C_AXIS:
		NOS_L=NOS_CL;
		break;
		case ABC_AXIS:
		//NOS_L=NOS_AL+NOS_BL+NOS_CL
		break;	
	}

	for (int n = 0; n<NOS_L; n++)//NOS_L number of spins per slic layer.
	{
		i = n*(VerN/3); // shift in vertex array 
		j = n*Kinp; // shift in index array (for arrow and cone vertex num<index num!)
		for (int k=0; k<Kinp; k++) Iout[j+k] = Iinp[k]+i;		
	}
}



void 
UpdateVerticesNormalsColors (float * Vinp, float * Ninp, int Kinp, 
							float * Vout, float * Nout, float * Cout, int Kout, 
							float * Px, float * Py, float * Pz,
							double * Sx, double * Sy, double * Sz, int mode)
{
	//float tmpV1[3], tmpV2[3], tmpV3[3], RGB[3], U, A;
	float tmpV2[3], RGB[3], U, A;
	int i,j,n;
	int anini=0;
	int anfin=ABC[0];
	int bnini=0;
	int bnfin=ABC[1];
	int cnini=0;
	int cnfin=ABC[2];

	switch( WhichSliceMode)
	{
		case A_AXIS:
			anini=(Alayer-1);
	        anfin=Alayer;
		break;
		case B_AXIS:
			bnini=(Blayer-1);
	        bnfin=Blayer;
		break;
		case C_AXIS:
			cnini=(Clayer-1);
	        cnfin=Clayer;
		break;
		case ABC_AXIS:
            //
		break;

	}
	j=-1;
	switch (mode)
	{
	case ARROW1:
	case CONE1:
		for (int an = anini; an<anfin; an++) // n runs over spins (ARROW1)
		{
		for (int bn = bnini; bn<bnfin; bn++) // n runs over spins (ARROW1)
		{
		for (int cn = cnini; cn<cnfin; cn++) // n runs over spins (ARROW1)
		{
			n = an+bn*ABC[0]+cn*ABC[0]*ABC[1];
			//slow version is commented but easy to read: 
			tmpV2[0] = Sx[n];
			tmpV2[1] = Sy[n];
			tmpV2[2] = Sz[n];
	        HSVtoRGB( tmpV2, RGB, InvertValue, InvertHue);
			j++;
			for (int k=0; k<Kinp/3; k++) // k runs over vertices of the arrow/cone 
			{
				i = j*Kinp + 3*k;	// vertex index
				//slow version is commented but easy to read: 
				// tmpV1[0] = Vinp[3*k+0];
				// tmpV1[1] = Vinp[3*k+1];
				// tmpV1[2] = Vinp[3*k+2];
				// NewBasisCartesian(tmpV1, tmpV2, tmpV3); // to find arrow vector components w.r.t basis based on Sx,Sy,Sz
				// Vout[i+0] = tmpV3[0]*Scale + Px[n];	// new x-component of vertex + translation
				// Vout[i+1] = tmpV3[1]*Scale + Py[n];	// new y-component of vertex + translation
				// Vout[i+2] = tmpV3[2]*Scale + Pz[n];	// new z-component of vertex + translation
				U = Sx[n]*Sx[n]+Sy[n]*Sy[n]+(1e-37f); 
				
				A = (-Sy[n]*Vinp[3*k+0] + Sx[n]*Vinp[3*k+1])*(1. - Sz[n])/U; 

				Vout[i+0] =(-Sy[n]*A + Vinp[3*k+0]*Sz[n] + Sx[n]*Vinp[3*k+2]			)*Scale + Px[n];
				Vout[i+1] =( Sx[n]*A + Vinp[3*k+1]*Sz[n] + Sy[n]*Vinp[3*k+2]			)*Scale + Py[n];
				Vout[i+2] =( Vinp[3*k+2]*Sz[n] - (Sx[n]*Vinp[3*k+0]+Sy[n]*Vinp[3*k+1])	)*Scale + Pz[n];	

				//slow version is commented but easy to read:
				// tmpV1[0] = Ninp[3*k+0];
				// tmpV1[1] = Ninp[3*k+1];
				// tmpV1[2] = Ninp[3*k+2];
				//NewBasisCartesian(tmpV1, tmpV2, tmpV3);	// to find vertices normals w.r.t basis based on Sx,Sy,Sz
				// Nout[i+0] = tmpV3[0];		// x-component of vertex normal
				// Nout[i+1] = tmpV3[1];		// y-component of vertex normal
				// Nout[i+2] = tmpV3[2];		// z-component of vertex normal

				A = (-Sy[n]*Ninp[3*k+0] + Sx[n]*Ninp[3*k+1])*(1. - Sz[n])/U; 

				Nout[i+0] =-Sy[n]*A + Ninp[3*k+0]*Sz[n] + Sx[n]*Ninp[3*k+2];
				Nout[i+1] = Sx[n]*A + Ninp[3*k+1]*Sz[n] + Sy[n]*Ninp[3*k+2];
				Nout[i+2] = Ninp[3*k+2]*Sz[n] - (Sx[n]*Ninp[3*k+0]+Sy[n]*Ninp[3*k+1]);

				Cout[i+0] = RGB[0];			// x-component of vertex normal
				Cout[i+1] = RGB[1];			// y-component of vertex normal
				Cout[i+2] = RGB[2];			// z-component of vertex normal
			}
		}
		}
		}
	break;

	case POINT:
		for (int an = anini; an<anfin; an++) // n runs over spins (ARROW1)
		{
		for (int bn = bnini; bn<bnfin; bn++) // n runs over spins (ARROW1)
		{
		for (int cn = cnini; cn<cnfin; cn++) // n runs over spins (ARROW1)
		{
			n = an+bn*ABC[0]+cn*ABC[0]*ABC[1];
			tmpV2[0] = Sx[n];
			tmpV2[1] = Sy[n];
			tmpV2[2] = Sz[n];
	        HSVtoRGB( tmpV2, RGB, InvertValue, InvertHue);
	        j++;
			i = j*Kinp;			// index of first cane vertex 
			Vout[i+0] = Px[n];	// new x-component of vertex + translation
			Vout[i+1] = Py[n];	// new y-component of vertex + translation
			Vout[i+2] = Pz[n];	// new z-component of vertex + translation
			Cout[i+0] = RGB[0];	// x-component of vertex normal
			Cout[i+1] = RGB[1];	// y-component of vertex normal
			Cout[i+2] = RGB[2];	// z-component of vertex normal	
		}
		}
		}

	break;

	case CANE:
	default:
		for (int an = anini; an<anfin; an++) // n runs over spins (ARROW1)
		{
		for (int bn = bnini; bn<bnfin; bn++) // n runs over spins (ARROW1)
		{
		for (int cn = cnini; cn<cnfin; cn++) // n runs over spins (ARROW1)
		{
			n = an+bn*ABC[0]+cn*ABC[0]*ABC[1];
			
			tmpV2[0] = Sx[n];
			tmpV2[1] = Sy[n];
			tmpV2[2] = Sz[n];
	        HSVtoRGB( tmpV2, RGB, InvertValue, InvertHue);
	        j++;
			//i = (n-nini)*Kinp;							// index of ferst cane vertex 
			i = j*Kinp;
			Vout[i+0] = tmpV2[0]*(1-Pivot)*Scale + Px[n];	// new x-component of vertex + translation
			Vout[i+1] = tmpV2[1]*(1-Pivot)*Scale + Py[n];	// new y-component of vertex + translation
			Vout[i+2] = tmpV2[2]*(1-Pivot)*Scale + Pz[n];	// new z-component of vertex + translation
			//i = n*Kinp/3*4;		// colors contains 4 floats
			Cout[i+0] = RGB[0];					// x-component of vertex normal
			Cout[i+1] = RGB[1];					// y-component of vertex normal
			Cout[i+2] = RGB[2];					// z-component of vertex normal
			//Cout[i+3] = 1.f;
			//printf( "|V1=%f,%f,%f \n",Vout[i+0],Vout[i+1],Vout[i+2]);
			//i = (n-nini)*Kinp + 3*1;					// index of ferst cane vertex 
			i = j*Kinp+ 3*1;
			Vout[i+0] = -tmpV2[0]*(Pivot)*Scale + Px[n];		// new x-component of vertex + translation
			Vout[i+1] = -tmpV2[1]*(Pivot)*Scale + Py[n];		// new y-component of vertex + translation
			Vout[i+2] = -tmpV2[2]*(Pivot)*Scale + Pz[n];		// new z-component of vertex + translation
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



void UpdateSpinPositions(float abc[][3], int ABC[3], float BD[][3], int NBD, float box[][3], float * Px, float * Py, float * Pz)
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
	int i=-1; // spin index
	for( int L=0;L<ABC[2];L++)// translation of basic domain along vector 'c' L times
	{
		for(int K=0;K<ABC[1];K++)// translation of basic domain along vector 'b' K times
		{
			for(int J=0;J<ABC[0];J++) // translation of basic domain along vector 'a' J times
			{	
				for(int I=0; I < NBD; I++) // runs over atoms in basic domain 
				{	
					i++; 
					Px[i] = BD[I][0] + abc[0][0]*J + abc[1][0]*K + abc[2][0]*L-Tr[0]; 
					Py[i] = BD[I][1] + abc[0][1]*J + abc[1][1]*K + abc[2][1]*L-Tr[1];
					Pz[i] = BD[I][2] + abc[0][2]*J + abc[1][2]*K + abc[2][2]*L-Tr[2];
				}
			}
		}
	}	
}
///////////////////////////////////////////////////////////////////////////////////////////////////

void
ReallocateArrayDrawing(int faces, int mode)
{
	free(vertexProto); free(normalProto); free(indicesProto);
	free(vertices); free(normals); free(colors); free(indices);
	int NOS_L=NOS_CL;
	switch( WhichSliceMode)
	{
		case A_AXIS:
		NOS_L=NOS_AL;
		break;
		case B_AXIS:
		NOS_L=NOS_BL;
		break;
		case C_AXIS:
		NOS_L=NOS_CL;
		break;
		case ABC_AXIS:
		//
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
			// arrowFaces - number of arrow faces
			ElNumProto = 2*arrowFaces-2; // number of triangles per arrow
			IdNumProto = 3*ElNumProto; // number of indixes per arrow
			VCNumProto = 3*((arrowFaces)+3*arrowFaces);// number of vertecies per arrow			
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

		case POINT:		
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
			vertices	= (float  *)calloc(VCNum, sizeof( float  ));
			normals 	= NULL;
			colors 		= (float  *)calloc(VCNum, sizeof( float  ));
			indices		= (GLuint *)calloc(IdNum, sizeof( GLuint ));	
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

void 
CreateNewVBO( )
{
	glGenBuffers(1, &vboIdV);
	glGenBuffers(1, &vboIdN);
	glGenBuffers(1, &vboIdC);
	glGenBuffers(1, &iboIdI);
}

void UpdateVBO(GLuint * V, GLuint * N, GLuint * C, GLuint * I, float * ver, float * nor, float * col, GLuint * ind)
{	//ver, nor, col and ind pointer to arrays of vertxcies components, norlamls, colors and indecies 
	switch (WhichVectorMode)
	{
		case ARROW1:
		case CONE1:			
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
				glBufferData(GL_ELEMENT_ARRAY_BUFFER, IdNum * sizeof(GLuint), NULL, GL_DYNAMIC_DRAW);
				glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, IdNum * sizeof(GLuint), ind);
		break;
		case POINT:
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

void 
drawVBO(int VCNum, int mode)
{
	switch (mode)
	{
		case ARROW1:
		case CONE1:
			glBindBuffer(GL_ARRAY_BUFFER, vboIdC);		glColorPointer(3, GL_FLOAT, 0, (void*)0);
			glBindBuffer(GL_ARRAY_BUFFER, vboIdN);		glNormalPointer(GL_FLOAT, 0, (void*)0);
			glBindBuffer(GL_ARRAY_BUFFER, vboIdV);		glVertexPointer(3, GL_FLOAT, 0, (void*)0);	

			glEnableClientState(GL_COLOR_ARRAY);		// enable color arrays
			glEnableClientState(GL_NORMAL_ARRAY);		// enable normal arrays
			glEnableClientState(GL_VERTEX_ARRAY);		// enable vertex arrays	

			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, iboIdI);
			glDrawElements(GL_TRIANGLES, VCNum, GL_UNSIGNED_INT, (void*)(0));

			glDisableClientState(GL_VERTEX_ARRAY);		// disable vertex arrays
			glDisableClientState(GL_NORMAL_ARRAY);		// disable normal arrays
			glDisableClientState(GL_COLOR_ARRAY);		// disable color arrays

			glBindBuffer(GL_ARRAY_BUFFER,			0);	// disable vertex arrays
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,	0);	// disable normal arrays

		break;

		case POINT:
		    glDisable(GL_LIGHTING);
			glBindBuffer(GL_ARRAY_BUFFER, vboIdC);		glColorPointer(3, GL_FLOAT, 0, (void*)0);
			glBindBuffer(GL_ARRAY_BUFFER, vboIdV);		glVertexPointer(3, GL_FLOAT, 0, (void*)0);	

			glEnableClientState(GL_COLOR_ARRAY);		// enable normal arrays
			glEnableClientState(GL_VERTEX_ARRAY);		// enable vertex arrays	

			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, iboIdI);
			
			glPointSize(10.f*Scale);
			//Note, the range specifies the range of index values in the index region being rendered from.
			//glDrawElements(GL_POINTS, iNum, GL_UNSIGNED_INT, (void*)(iStart*sizeof(GLuint)));
			glDrawElements(GL_POINTS, VCNum/3, GL_UNSIGNED_INT, (void*)0);//draw whole VBO

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
			glDrawElements(GL_LINES, VCNum/3, GL_UNSIGNED_INT, (void*)0);// draw a stick VCNum
			glPointSize(4.0f*Scale);
			//glDrawElements(GL_POINTS, iNum, GL_UNSIGNED_INT, (void*)(iStart*sizeof(GLuint)));// draw a ball (point at the end of the stick)
			glDrawElements(GL_POINTS, VCNum/3, GL_UNSIGNED_INT, (void*)0);// draw a ball (point at the end of the stick)

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
			glDrawElements(GL_POINTS, VCNum/3/2, GL_UNSIGNED_INT, (void*)0);// draw a ball (point at the end of the stick)

			glDisableClientState(GL_VERTEX_ARRAY);		// disable vertex arrays
			glDisableClientState(GL_COLOR_ARRAY);		// disable normal arrays

			glBindBuffer(GL_ARRAY_BUFFER,			0);	// disable vertex arrays
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,	0);	// disable normal arrays
			glEnable(GL_LIGHTING);
	}
}


