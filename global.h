// Constants:
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
// Global vatiables:

//  Cartesian components of spins:
double*         Sx; // array of spin x-component 
double*         Sy; // array of spin y-component
double*         Sz; // array of spin z-component

double*         tSx; // temporal array of spin x-component 
double*         tSy; // temporal array of spin y-component
double*         tSz; // temporal array of spin z-component

double*         t2Sx; // temporal array of spin x-component  RK45
double*         t2Sy; // temporal array of spin y-component  RK45
double*         t2Sz; // temporal array of spin z-component  RK45

double*         bSx; // buffer array for data transfer for visulaization
double*         bSy; // buffer array for data transfer for visulaization
double*         bSz; // buffer array for data transfer for visulaization

float*          IsoLineX;
float*          IsoLineY;
float*          IsoLineZ;

float*          RNx; // array of x-components for random vector 
float*          RNy; // array of y-components for random vector 
float*          RNz; // array of z-components for random vector 

float*          Px; // x spin position array
float*          Py; // y spin position array
float*          Pz; // z spin position array 

float*          BPx; // x Block position array
float*          BPy; // y Block position array
float*          BPz; // z Block position array 

int*            Kind; 

double*         Heffx; // x component position array
double*         Heffy; // y component position array
double*         Heffz; // z component position array

int             Boundary[3] = {1, 1, 1};// boundary conditions along a, b, c translation vectors

double*         Etot; // array of total energy per spin 
double*         Etot0; // array of total energy per spin 
double          totalEnergy;
double          perSpEnergy;
double          totalEnergyFerro;
double          perSpEnergyMinusFerro;
double          Mtot[3]; // Total magnetic moment
double          mtot[3]; // Total magnetic moment
double          Max_torque[THREADS_NUMBER];
double          MAX_TORQUE=0;

double          BigDataBank[10][1000];
int             recordsCounter=0;

//Color scheme variables
int             HueMapRGB[6]={0,60,120,180,240,300};// initial (equidistant RGB) hue map for the color sphere
int             HueMapRYGB[6]={0,90,180,225,270,315};// initial (equidistant RYGB) hue map for the color sphere
int             HueMap[6]={0,60,120,180,240,300};// current color huemap
float*          RHue;
float*          GHue;
float*          BHue;

//Allocate arrays for neighbours map
int             NeighborPairs;
int*            AIdxBlock;// index of the atom within the block
int*            NIdxBlock;// index of the neighbour within the block
int*            NIdxGridA;// index of the relative position of the block of the neighbour in the greed along tr. vect. a
int*            NIdxGridB;// index of the relative position of the block of the neighbour in the greed along tr. vect. b
int*            NIdxGridC;// index of the relative position of the block of the neighbour in the greed along tr. vect. c
int*            SIdx     ;// index of the shell corresponding to this pair
float           Box[3][3];// box for the simulated domain 
float*          Jexc;//array of isotropic Heisenberg exchange constants for each pair
float*          Bexc;//array of bi-quadratic exchange constants for each pair
float*          Dexc;//array of DMI constants (module of DMI vectors) for each pair

float*          VDMx;//array of Dzyaloshinskii vectors x-component (normalized)
float*          VDMy;//array of Dzyaloshinskii vectors y-component (normalized)
float*          VDMz;//array of Dzyaloshinskii vectors z-component (normalized)
//Heisenberg exchange
float           Jij[]={     // Jij[shell]
                1.0,    // first shell
                1.0/6.0, // second shell
                -1.0/3.0,  // third shell
                0.0,    // fourth shell
                0.0,
                0.0
                };
//bi-quadratic exchange
float           Bij[]={     // Bij[shell]
                0.0,    // first shell
                0.0,    // second shell
                0.0,    // third shell
                0.0,    // fourth shell
                0.0,
                0.0
                };
//Dzyaloshinskii-Moriya Interaction
float           Dij[]={ // Dij[shell] abs value for DMI vector 
                0.125,//12pi350.20943951023932,// 2pi/30
                0.0,//0.483321946706122,//0.41887902,//0.0369138485,    // first shell
                0.0,//0.1,  // second shell
                0.0,//0.085,    // third shell
                0.0,//0.024,    // fourth shell
                0.0,
                0.0
                };
//Magnetocrystalline anisotropy:
float           VKu1[]={    0.0 , 0.0, 1.0 }; // uniaxial anisotropy vector
float           VKu2[]={    0.0 , 0.0, 1.0 }; // easy plane anisotropy vector
float           Ku1 = 0.0;//uniaxial anisotropy constant
float           Ku2 = 0.0;//uniaxial anisotropy constant
float           Kc = 0;//cubic anisotropy constant 
//DC applied H-field:
float           VHf[]={ 0.0 , 0.0, 1.0 };
float           VHtheta=0;//[111] T=54.74, F=45
float           VHphi=0;
// float*      VHf=(float *)calloc(3, sizeof(float));
float           Hf=0.0;//0.01632;
float           Bdc[]={ Hf*VHf[0] , Hf*VHf[1], Hf*VHf[2] };
//AC applied H-field:
float           VHac[]={ 1.0 , 0.0, 0.0 };
float           Hac=0.001;//001632;
float           Bac[]={ 0.0 , 0.0, 0.0 };//Bac is components of external field for info panel only.
float           Period_dc=244.994;
float           Omega_dc=TPI/Period_dc;
float           GPulseWidth=20.0;
float           t_offset=80;
float           HacTime=0.0;//time (iteration) dependent value of ac field
enum            enACField{SIN_FIELD, GAUSSIAN_FIELD} ; // which mode
enACField       WhichACField = SIN_FIELD;   // RND by default 

enum            Average_mode{ALONG_A,ALONG_B, ALONG_C, ALONG_0};
int             WhichAverageMode = ALONG_0;//ALONG_0 means do no average.
int             save_slice=0;
//Current polarization direction
float           VCu[]={ 0.0 , 0.0, 1.0 };
//Spin-torque parameter ~ density of ingected current
float           Cu = 0.0;

//on/off precession term
int             Precession=1;
//damping parameter:
float           damping=1.0;
//timestep
float           t_step=0.3;//1.0/16.0;
//temperature
float           Temperature=0.0;

//GUI control
int             Play=0;
int             SpecialEvent=1;
int             DataTransfer=1;//metka possibly to delete
unsigned int    ITERATION=0; 
unsigned int    Max_Numb_Iteration=1000000;
int             Record=0;// record <sx>, <sy>, <sz> into fole sxsysz.csv
int             AC_FIELD_ON=0;//ON/OFF AC field signal.

// FPS & IPS
int             currentTime=0;
int             previousTime=0;
int             frameCount=0;
int             timeInterval=0;
int             previousIteration=0; 
int             currentIteration=0;
float           FPS, IPS;
FILE*           outFile;//table output file
unsigned int    rec_iteration=1;//each rec_iteration one puts into sxsysz.csv file
char            BuferString[800];//for output file table// metka 80 -> 800 test LLG
double          outputEtotal;
double          outputMtotal[3];
int             SleepTime=1000;