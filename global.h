#include <math.h>
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

// Slicing parameters
int             A_layer_min = 1;    // which layer (along a tr. vect. ) to show max=uABC[0]
int             B_layer_min = 1;    // which layer (along b tr. vect. ) to show max=uABC[1]
int             C_layer_min = 1;    // which layer (along c tr. vect. ) to show max=uABC[2]
int             A_layer_max = 1;    // which layer (along a tr. vect. ) to show max=uABC[0]
int             B_layer_max = 1;    // which layer (along b tr. vect. ) to show max=uABC[1]
int             C_layer_max = 1;    // which layer (along c tr. vect. ) to show max=uABC[2]

typedef enum    {A_AXIS, B_AXIS, C_AXIS, FILTER} enSliceMode; // which mode
enSliceMode     WhichSliceMode  = C_AXIS; 



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

double*         t3Sx; // temporal array of spin x-component  RK45
double*         t3Sy; // temporal array of spin y-component  RK45
double*         t3Sz; // temporal array of spin z-component  RK45

double**        Image_x; // Snepshots of the system
double**        Image_y; // Snepshots of the system
double**        Image_z; // Snepshots of the system

double**        dImage_x; // Snepshots of dm
double**        dImage_y; // Snepshots of dm
double**        dImage_z; // Snepshots of dm

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

bool*           Proj;
double          ALPHA;

int             Boundary[3] = {1, 1, 0};// boundary conditions along a, b, c translation vectors

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

double          coef[THREADS_NUMBER][6];

double          BigDataBank[10][1500];
int             recordsCounter=0;
int             WhereAmI = 0;
int             NumOfPoints = 10000;
int             TempPoints = 10;
int             TempInt = 1;

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
                0.0,
                0.0,//-1.0/16.0, // second shell
                0.0,  // third shell
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
                0.251327412, //LD=64//0.196349541,//LD=32//,//TPI/32,//(TPI/12.0),//12pi350.20943951023932,// 2pi/30
                0.0,//-(TPI/12.0)/8,//0.1,  // second shell
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
float           VHtheta=0.0;//[111] T=54.74, F=45
float           VHphi= 0;
// float*      VHf=(float *)calloc(3, sizeof(float));
float           Hf= 0.0;//0.0134936;//0.006149226;//0.006264886;//0.00607212;//0.02428848;////0.65*Dij[0]*Dij[0]/Jij[0];
                //Theta=3 degree Hf= 0.016342396686686
float           Bdc[]={ Hf*VHf[0] , Hf*VHf[1], Hf*VHf[2] };
//AC applied H-field:
float           VHac[]={ 0.0 , 0.0, 1.0 };
float           Hac=Hf*sin(1*PI/180);//0.1;//001632;
double          Bac[]={ 0.0 , 0.0, 0.0 };//Bac is components of external field for info panel only.
double          Omega_dc=0.005;//0.029481249;//0.014141413;
double          Period_dc=TPI/Omega_dc;
float           GPulseWidth=20.0;
float           t_offset=0.0;
float           HacTime=0.0;//time (iteration) dependent value of ac field
enum            enACField{SIN_FIELD, GAUSSIAN_FIELD, SINC_FIELD, CIRCULAR_FIELD} ; // which mode
enACField       WhichACField = SIN_FIELD;  

enum            IntegrationScheme{HEUN,SIB,RK23,RK45,RELAX};
int             WhichIntegrationScheme = SIB;

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
float           damping=0.01;//0.01;
//timesteps
float           t_step=0.1;//0.05;//1.0/16.0;
//temperature
float           Temperature = 0;


//Zhang-Li
float           Xi=0.0;//0.05;//nonadiabatisity parameter
float           Curr_u = 0.0;//-0.05;// effective current (u) in x direction

//GUI control
int             Play=0;
int             SpecialEvent=1;
int             DataTransfer=1;//metka possibly to delete
unsigned int    ITERATION=0; 
unsigned int    Max_Numb_Iteration=100000;//530000;
int             Record=0;// record <sx>, <sy>, <sz> into the file 
int             AC_FIELD_ON=0;//ON/OFF AC field signal.
int             AC_MODE_REC=0;

// FPS & IPS
int             currentTime=0;
int             previousTime=0;
int             frameCount=0;
int             timeInterval=0;
int             previousIteration=0; 
int             currentIteration=0;
float           FPS, IPS;
FILE*           outFile;//table output file
unsigned int    rec_iteration=100;//each rec_iteration one puts into sxsysz.csv file
int             rec_num_mode=10;//rec_num_mode is the number of snapshots over which the dynalical mode at this phase is averaged
int             current_rec_num_mode=0;
int             num_images=32;// phase at which one should record snapshot of the mode.
char            BuferString[800];//for output file table// metka 80 -> 800 test LLG
double          outputEtotal;
double          outputMtotal[3];
int             SleepTime=1000;

char            shortBufer[200];
char            inputfilename[64] = "input.csv";
char            outputfilename[64] = "output.csv";


////////////////////////////////////////////// GEOMETRY ////////////////////////////////////////////
//translation vectors cubic:
//float         abc[3][3] = {
//                              {       1.0f, 0.0f, 0.0f }, // a
//                              {       0.0f, 1.0f, 0.0f }, // b
//                              {       0.0f, 0.0f, 1.0f }};// c
                                
//atom positions in basic domain: {x0,y0,z0, x1,y1,z1,...}
//not more then 100 atoms per basic domain
//FCC
//float Block[4*3] = 
//              {       0.0f,   0.0f,   0.0f,   
//                      0.5f,   0.5f,   0.0f,
//                      0.5f,   0.0f,   0.5f,
//                      0.0f,   0.5f,   0.5f    
//              };
//BCC
/*
float Block[2*3] = 
    { 0.0f, 0.0f, 0.0f, 
      0.5f, 0.5f, 0.5f, 
    };
 */




//Simple Cubic 1 [001]
// float                abc[3][3] = {
//                              {       1.0f, 0.0f, 0.0f }, // a
//                              {       0.0f, 1.0f, 0.0f }, // b
//                              {       0.0f, 0.0f, 1.0f }};// c
// float                        Block[][3] = { 
//                              {0.5f, 0.5f, 0.5f},  
//                              };

/*float         abc[3][3] = {
                                {       3.0f, 0.0f, 0.0f }, // a
                                {       0.0f, 1.0f, 0.0f }, // b
                                {       0.0f, 0.0f, 10.0f }};// c
float                   Block[][3] = { 
                                {0.0f, 0.5f, 0.0f}, 
                                {1.5f, 0.5f, 0.0f}, 
                                {3.0f, 0.5f, 0.0f}, 
                                {0.75f, 0.0f, 0.0f},
                                {2.25f, 0.0f, 0.0f},
                                {0.75f, 1.0f, 0.0f},
                                {2.25f, 1.0f, 0.0f}  
                                };*/

// float                abc[3][3] = {
//                              {       1.0f, 0.0f, 0.0f }, // a
//                              {       0.0f, 1.0f, 0.0f }, // b
//                              {       0.0f, 0.0f, 1.0f }};// c
// float                        Block[][3] = { 
//                              {0.5f, 0.5f, 0.5f},  
//                              };
// Simple Cubic 2 [011]
// float                abc[3][3] = {
//                              {       1.0f, 0.0f, 0.0f }, // a
//                              {       0.0f, sqrt(2.f), 0.0f }, // b
//                              {       0.0f, 0.0f, sqrt(2.f) }};// c                           
// float                Block[][3] = {          
//                              {0,                   0,        0},     
//                              {0.5,   sqrt(2.f)/2,    0}
//                              };      

//Simple Cubic 3 [111]
// float                abc[3][3] = {
//                              {       sqrt(2.f), 0.0f, 0.0f }, // a
//                              {       0.0f, sqrt(3.f)*sqrt(2.f), 0.0f }, // b
//                              {       0.0f, 0.0f, sqrt(3.f) }};// c                           
// float                Block[][3] = {          
//                              {0,                   0,                0},     
//                              {sqrt(2.f)/2, sqrt(6.f)/2,      0},
//                              {sqrt(2.f)/2, sqrt(6.f)/6,  sqrt(3.f)/3.f},
//                              {0, 2*sqrt(6.f)/3,  sqrt(3.f)/3.f},
//                              {0, 1*sqrt(6.f)/3, -sqrt(3.f)/3.f},
//                              {sqrt(2.f)/2, 5*sqrt(6.f)/6, -sqrt(3.f)/3.f}
//                              };

//FCC 2
// float                abc[3][3] = {
//                              {       1.0f, 0.0f, 0.0f }, // a
//                              {       0.0f, sqrt(2.f), 0.0f }, // b
//                              {       0.0f, 0.0f, sqrt(2.f) }};// c                           
// float                Block[][3] = {          
//                              {0,                   0,                  0},   
//                              {0,     sqrt(2.f)/2,                      0},
//                              {0.5,   sqrt(2.f)/4,sqrt(2.f)/4},
//                              {0.5, 3*sqrt(2.f)/4, sqrt(2.f)/4},
//                              {0,     sqrt(2.f)/2,    sqrt(2.f)/2}
//                              };      
//FCC 3
// float                abc[3][3] = {
//                              {sqrt(2.f), 0.0f, 0.0f }, // a
//                              {sqrt(2.f)/2, sqrt(6.f)/2, 0.0f }, // b
//                              {sqrt(2.f)/2, sqrt(6.f)/6, sqrt(3.f)/sqrt(2.f) }};// c                          
// float                Block[][3] = {          
//                              {0,                   0,                  0},   
//                              };


//B20(1)
// float                abc[3][3] = {
//                              {       1.0f, 0.0f, 0.0f }, // a
//                              {       0.0f, 1.0f, 0.0f }, // b
//                              {       0.0f, 0.0f, 1.0f }};// c
// #define uB20         0.136506//0.138f//0.133
// float                Block[][3] = {          
//                              {   uB20,                  uB20,            uB20},      
//                              {0.5f+uB20,     0.5f-uB20,      1.0f-uB20},
//                              {1.0f-uB20,     0.5f+uB20,      0.5f-uB20},     
//                              {0.5f-uB20,     1.0f-uB20,      0.5f+uB20}      
//                              };
//B20(2)

float           abc[3][3] = {
                                {       1.0f, 0.0f, 0.0f }, // a
                                {       0.0f, 1.0f, 0.0f }, // b
                                {       0.0f, 0.0f, 1.0f }};// c

float           Block[][3] = {          
                                {   0.5,    0.5,     0.5}        
                                };
// float           abc[3][3] = {
//                                 {       sqrt(3), 0.0f, 0.0f }, // a
//                                 {       0.0f, 1.0f, 0.0f }, // b
//                                 {       0.0f, 0.0f, 1.0f }};// c

// float           Block[][3] = {          
//                                 {   0.0,    0.0,     0.0},
//                                 {   sqrt(3)/2,    0.5,     0.0}        
//                                 };

// #define uB20         0.138f
// float                Block[][3] = {          
//                              {   0.,            0.,      0.},        
//                              {0.5f,  0.5f-2*uB20,    1.0f-2*uB20},
//                              {1.0f-2*uB20,   0.5f,   0.5f-2*uB20},   
//                              {0.5f-2*uB20,   1.0f-2*uB20,    0.5f}   
//                              };

//B20(Maria)
// float                abc[3][3] = {
//                              {       1.0f, 0.0f, 0.0f }, // a
//                              {       0.0f, 1.0f, 0.0f }, // b
//                              {       0.0f, 0.0f, 1.0f }};// c
// //#define uB20               0.138f//MnSi
// // float             Block[][3] = {          
// //                           {       0.0f,        0.0f,        0.0f},//r1= u, u, u   
// //                           {       0.5f, 0.5f-2*uB20,     -2*uB20},//r2= 0.5+u, 0.5-u,-u
// //                           {    -2*uB20,        0.5f, 0.5f-2*uB20},//r3= -u, 0.5+u, 0.5 -u 
// //                           {0.5f-2*uB20,     -2*uB20,        0.5f} //r4=0.5-u, -u, 0.5+u
// //                           };
// #define uB20         0.135f//FeGe
// float                Block[][3] = {          
//                              {       0.0f,        0.0f,        0.0f},//r1= u, u, u   
//                              {       0.5f, 0.5f-2*uB20,     -2*uB20},//r2= 0.5+u, 0.5-u,-u
//                              {    -2*uB20,        0.5f, 0.5f-2*uB20},//r3= -u, 0.5+u, 0.5 -u 
//                              {0.5f-2*uB20,     -2*uB20,        0.5f} //r4=0.5-u, -u, 0.5+u
//                              };

//number of translations for the basic domain along a,b, and c verctors respectively 
int             uABC[3] = {10,10,10};//Grid dimensionality along translation vectors a, b, c; uABC[i]>=1 
int             ShellNumber = 1;
int             AtomsPerBlock = sizeof(Block)/sizeof(float)/3;
float*          RadiusOfShell = (float *)calloc(ShellNumber , sizeof(float));  
int*            NeighborsPerAtom = (int *)calloc(AtomsPerBlock, sizeof(int));
// total number of neighbour pairs per whole map of neighbours
int             NOS=AtomsPerBlock*uABC[0]*uABC[1]*uABC[2]; // number of spins
int             NOS_AL=AtomsPerBlock*uABC[1]*uABC[2]; // number of spins per A layer
int             NOS_BL=AtomsPerBlock*uABC[0]*uABC[2]; // number of spins per B layer
int             NOS_CL=AtomsPerBlock*uABC[0]*uABC[1]; // number of spins per C layer

int             NOSK = 0;

double          iNOS = 1.0/NOS;

int             NOB=uABC[0]*uABC[1]*uABC[2]; // number of Blocks
int             NOB_AL=uABC[1]*uABC[2]; // number of blocks per A layer
int             NOB_BL=uABC[0]*uABC[2]; // number of blocks per B layer
int             NOB_CL=uABC[0]*uABC[1]; // number of blocks per C layer

/////////////////////////////////////////////Functions//////////////////////////////

//Color map control parameters
int             ColorShift=0; // by default the red collor corresponds to phi=0
int             InvertHue=0; // by default the sequence is red-green-blue for phi=0,120,240(-120) degrees respectively
int             InvertValue=0; //n_z=+1 (white), -1 (black). For InvertHue=1 vice versa 





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


void HSVtoRGB(float Vec[3], float rgb[3], int inV, int inH ){
    // int H = inH*359+(1-2*inH)*((int) atan2int (Vec[1]+0.01,Vec[0])); //it's fast atan function [int deg], see MATH.cpp
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


    rgb[0] = color_function(H+120+ColorShift)*dV+minV;//rad
    rgb[1] = color_function(H+000+ColorShift)*dV+minV;//green
    rgb[2] = color_function(H-120+ColorShift)*dV+minV;//blue
    
    // rgb[0] = colorGreen(H+120+ColorShift)*dV*S+minV+dV*(1-S)*0.5;//rad
    // rgb[1] = colorGreen(H    +ColorShift)*dV*S+minV+dV*(1-S)*0.5;//green
    // rgb[2] = colorGreen(H-120+ColorShift)*dV*S+minV+dV*(1-S)*0.5;//blue
}




void Save_OVF_b8(double* Sx, double* Sy, double* Sz, char ovf_filename[64]){
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

         temp0  = abc[0][0]*abc[0][0];
         temp0 += abc[0][1]*abc[0][1];
         temp0 += abc[0][2]*abc[0][2];
         temp1 = sqrt(temp0);
         snprintf(shortBufer,80,"# xmax: %.6g\n",uABC[0]*temp1*a_lattice);
         fputs (shortBufer,pFile);

         temp0  = abc[1][0]*abc[1][0];
         temp0 += abc[1][1]*abc[1][1];
         temp0 += abc[1][2]*abc[1][2];
         temp2 = sqrt(temp0);
         snprintf(shortBufer,80,"# ymax: %.6g\n",uABC[1]*temp2*a_lattice);
         fputs (shortBufer,pFile);

         temp0  = abc[2][0]*abc[2][0];
         temp0 += abc[2][1]*abc[2][1];
         temp0 += abc[2][2]*abc[2][2];
         temp3 = sqrt(temp0);            
         snprintf(shortBufer,80,"# zmax: %.6g\n",uABC[2]*temp3*a_lattice);
         fputs (shortBufer,pFile);
         fputs ("# valuedim: 3\n",pFile);
         fputs ("# valuelabels: m_x m_y m_z\n",pFile);
         fputs ("# valueunits: 1 1 1\n",pFile);
         fputs ("# Desc: Total simulation time:  0  s\n",pFile);

         snprintf(shortBufer,80,"# xbase: %.6g\n",temp1*0.5*a_lattice);
         fputs (shortBufer,pFile);           

         snprintf(shortBufer,80,"# ybase: %.6g\n",temp2*0.5*a_lattice);
         fputs (shortBufer,pFile);

         snprintf(shortBufer,80,"# zbase: %.6g\n",temp3*0.5*a_lattice);
         fputs (shortBufer,pFile);

         snprintf(shortBufer,80,"# xnodes: %d\n",uABC[0]);
         fputs (shortBufer,pFile);
         snprintf(shortBufer,80,"# ynodes: %d\n",uABC[1]);
         fputs (shortBufer,pFile);
         snprintf(shortBufer,80,"# znodes: %d\n",uABC[2]);
         fputs (shortBufer,pFile);

         snprintf(shortBufer,80,"# xstepsize:  %.6g\n",temp1*a_lattice);
         fputs (shortBufer,pFile);           

         snprintf(shortBufer,80,"# ystepsize: %.6g\n",temp2*a_lattice);
         fputs (shortBufer,pFile);

         snprintf(shortBufer,80,"# zstepsize: %.6g\n",temp3*a_lattice);
         fputs (shortBufer,pFile);           

         fputs ("# End: Header\n",pFile);
         fputs ("# Begin: Data Binary 8\n",pFile);
         double Temp1[]= {123456789012345.0};
         fwrite (Temp1, sizeof(double), 1, pFile);
         for (int cn = 0; cn<uABC[2]; cn++)
            {
             for (int bn = 0; bn<uABC[1]; bn++)
                {
                 for (int an = 0; an<uABC[0]; an++)
                    {
                     int n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
                     n = n*AtomsPerBlock;//index of the first spin in the block
                     for (int atom=0; atom<AtomsPerBlock; atom++)
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


void SaveBin(double* Sx, double* Sy, double* Sz, char bin_filename[64]){
    unsigned short int num = 65535;
    struct tfshortint {
        unsigned short int t;
        unsigned short int f;
    };
    FILE * pFile = fopen (bin_filename,"wb");
    for(int k = 0; k < uABC[2]; k++){
        for(int j = 0; j < uABC[1]; j++){
            for(int i = 0; i < uABC[0]; i++){
            int n = i+j*uABC[0]+k*uABC[0]*uABC[1];
        
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


void SaveBmp(double* Sx, double* Sy, double* Sz, char bmp_filename[64], enSliceMode WhichSliceMode, int x1, int y1, int z1){
  if (WhichSliceMode!=FILTER)
  {
    int scale = 10;
    int N1 = uABC[0], N2 = uABC[1];
    switch (WhichSliceMode){
      case A_AXIS:
          N1 = uABC[1]; N2 = uABC[2];
          break;
      case B_AXIS:
          N1 = uABC[0]; N2 = uABC[2];
          break;
      case C_AXIS:
          N1 = uABC[0]; N2 = uABC[1];
          break; 
      default: N1 = uABC[0]; N2 = uABC[1];
          break;      
    }
        
    bitmap_image image(scale*N1,scale*N2);
    // set background to black
    image.set_all_channels(0, 0, 0);
    // image_drawer draw(image);
    float rgb[3],vec[3]; 
    for(int i = 0; i < N1; i++){
      for(int j = 0; j < N2; j++){
        int n;
        switch (WhichSliceMode){
          case A_AXIS:
              n = x1 + i*uABC[0] + j*uABC[0]*uABC[1];
              break;
          case B_AXIS:
              n = N1-i-1 + y1*uABC[0] + j*uABC[0]*uABC[1];
              break;
          case C_AXIS:
              n = i + j*uABC[0] + z1*uABC[0]*uABC[1];
              break; 
          default: n = i + j*uABC[0] + 0*uABC[0]*uABC[1];
              break;      
        }
        rgb[0] = 0; rgb[1] = 0; rgb[2] = 0;
        vec[0] = Sx[n]; vec[1] = Sy[n]; vec[2] = Sz[n];
        // vec[0] = 0; vec[1] = 0; vec[2] =1;
        HSVtoRGB(vec, rgb, InvertValue, InvertHue);///1,0
        rgb[0] *= 255;
        rgb[1] *= 255;
        rgb[2] *= 255;

        int ii = scale*i;
        for(int k1 = 0;k1<scale;k1++)
            for(int k2 = 0;k2<scale; k2++)
                image.set_pixel(ii+k1,   scale*((N2-1)-j)+k2,   (unsigned char)rgb[0], (unsigned char)rgb[1], (unsigned char)rgb[2]);
      }
    }
    image.save_image(bmp_filename);
    printf("Image has been saved to %s\n",bmp_filename);
  }
}

void Save_VTS_b4(double* Sx, double* Sy, double* Sz, float * Px, float * Py, float * Pz, float box[][3], char vts_filename[64]){
    // float temp0;
    // float temp1;
    // float temp2;
    // float temp3;
    float a_lattice = 1.0e-9; 
    //char ovf_filename[64] = "mode.ovf";
    FILE * pFile = fopen (vts_filename,"wb");
    if(pFile!=NULL) {   
        fputs ("<?xml version=\"1.0\"?>\n",pFile);
        fputs ("<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n",pFile);

        // temp0  = abc[0][0]*abc[0][0];
        // temp0 += abc[0][1]*abc[0][1];
        // temp0 += abc[0][2]*abc[0][2];
        // temp1 = sqrt(temp0);

        // temp0  = abc[1][0]*abc[1][0];
        // temp0 += abc[1][1]*abc[1][1];
        // temp0 += abc[1][2]*abc[1][2];
        // temp2 = sqrt(temp0); 

        // temp0  = abc[2][0]*abc[2][0];
        // temp0 += abc[2][1]*abc[2][1];
        // temp0 += abc[2][2]*abc[2][2];
        // temp3 = sqrt(temp0);  
        snprintf(shortBufer,80,"\t<StructuredGrid WholeExtent=\"0 %d 0 %d 0 %d\">\n",uABC[0]-1, uABC[1]-1, uABC[2]-1);
        fputs (shortBufer,pFile);
        snprintf(shortBufer,80,"\t\t<Piece Extent=\"0 %d 0 %d 0 %d \">\n",uABC[0]-1, uABC[1]-1, uABC[2]-1);
        fputs (shortBufer,pFile);
        fputs ("\t\t\t<PointData Vectors=\"m\">\n",pFile);
        fputs ("\t\t\t\t<DataArray type=\"Float32\" Name=\"m\" NumberOfComponents=\"3\" format=\"binary\">\n",pFile);
        fputs ("\t\t\t\t\t",pFile);
        float Tr[3] = { 
                        (box[0][0]+box[1][0]+box[2][0])/2.f,
                        (box[0][1]+box[1][1]+box[2][1])/2.f,
                        (box[0][2]+box[1][2]+box[2][2])/2.f
                      };
        for (int cn = 0; cn<uABC[2]; cn++){
            for (int bn = 0; bn<uABC[1]; bn++){
                for (int an = 0; an<uABC[0]; an++){
                    int n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
                    n = n*AtomsPerBlock;//index of the first spin in the block
                    for (int atom=0; atom<AtomsPerBlock; atom++){
                        int N = n + atom;
                        float Temp[]= {(Px[N]+Tr[0])*a_lattice, (Py[N]+Tr[1])*a_lattice, (Pz[N]+Tr[2])*a_lattice};  
                        // fwrite (Temp , sizeof(float), 3, pFile);
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
        for (int cn = 0; cn<uABC[2]; cn++){
            for (int bn = 0; bn<uABC[1]; bn++){
                for (int an = 0; an<uABC[0]; an++){
                    int n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
                    n = n*AtomsPerBlock;//index of the first spin in the block
                    for (int atom=0; atom<AtomsPerBlock; atom++){
                        int N = n + atom;
                        float Temp[]= {(float)Sx[N], (float)Sy[N], (float)Sz[N]}; 
                        // fwrite (Temp , sizeof(float), 3, pFile);
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

void Save_VTS_ascii(double* Sx, double* Sy, double* Sz, float * Px, float * Py, float * Pz, float box[][3], char vts_filename[64]){
    // float temp0 = 0;
    // float temp1 = 0;
    // float temp2 = 0;
    // float temp3 = 0;
    float a_lattice = 1.0;//.0e-9; 
    //char ovf_filename[64] = "mode.ovf";
    FILE * pFile = fopen (vts_filename,"wb");
    if(pFile!=NULL) {   
        fputs ("<?xml version=\"1.0\"?>\n",pFile);
        fputs ("<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n",pFile);

        // temp0  = abc[0][0]*abc[0][0];
        // temp0 += abc[0][1]*abc[0][1];
        // temp0 += abc[0][2]*abc[0][2];
        // temp1 = sqrt(temp0);

        // temp0  = abc[1][0]*abc[1][0];
        // temp0 += abc[1][1]*abc[1][1];
        // temp0 += abc[1][2]*abc[1][2];
        // temp2 = sqrt(temp0); 

        // temp0  = abc[2][0]*abc[2][0];
        // temp0 += abc[2][1]*abc[2][1];
        // temp0 += abc[2][2]*abc[2][2];
        // temp3 = sqrt(temp0);  
        snprintf(shortBufer,80,"\t<StructuredGrid WholeExtent=\"0 %d 0 %d 0 %d\">\n",uABC[0]-1, uABC[1]-1, uABC[2]-1);
        fputs (shortBufer,pFile);
        snprintf(shortBufer,80,"\t\t<Piece Extent=\"0 %d 0 %d 0 %d \">\n",uABC[0]-1, uABC[1]-1, uABC[2]-1);
        fputs (shortBufer,pFile);
        fputs ("\t\t\t<PointData Vectors=\"m\">\n",pFile);
        fputs ("\t\t\t\t<DataArray type=\"Float32\" Name=\"m\" NumberOfComponents=\"3\" format=\"ascii\">\n",pFile);
        fputs ("\t\t\t\t\t",pFile);
        float Tr[3] = { 
                        (box[0][0]+box[1][0]+box[2][0])/2.f,
                        (box[0][1]+box[1][1]+box[2][1])/2.f,
                        (box[0][2]+box[1][2]+box[2][2])/2.f
                      };
        for (int cn = 0; cn<uABC[2]; cn++){
            for (int bn = 0; bn<uABC[1]; bn++){
                for (int an = 0; an<uABC[0]; an++){
                    int n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
                    n = n*AtomsPerBlock;//index of the first spin in the block
                    for (int atom=0; atom<AtomsPerBlock; atom++){
                        int N = n + atom;
                        // float Temp[]= {(Px[N]+Tr[0])*a_lattice, (Py[N]+Tr[1])*a_lattice, (Pz[N]+Tr[2])*a_lattice};  
                        // fwrite (Temp , sizeof(float), 3, pFile);
                        snprintf(shortBufer,80,"%.6g %.6g %.6g ",(Px[N]+Tr[0])*a_lattice, (Py[N]+Tr[1])*a_lattice, (Pz[N]+Tr[2])*a_lattice);
                        fputs (shortBufer,pFile);
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
        for (int cn = 0; cn<uABC[2]; cn++){
            for (int bn = 0; bn<uABC[1]; bn++){
                for (int an = 0; an<uABC[0]; an++){
                    int n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
                    n = n*AtomsPerBlock;//index of the first spin in the block
                    for (int atom=0; atom<AtomsPerBlock; atom++){
                        int N = n + atom;
                        // float Temp[]= {Sx[N], Sy[N], Sz[N]}; 
                        // fwrite (Temp , sizeof(float), 3, pFile);
                        snprintf(shortBufer,80,"%.6g %.6g %.6g ",Sx[N], Sy[N], Sz[N]);
                        fputs (shortBufer,pFile);
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
float end_swap(const float val, const int size)
    {
    float ret = 0;
    for(int i = 0; i < size; ++i)
        ((char*)&ret)[size-1-i] = ((char*)&val)[i];
    return ret;
    };

double end_swap(const double val, const int size)
    {
    double ret = 0;
    for(int i = 0; i < size; ++i)
        ((char*)&ret)[size-1-i] = ((char*)&val)[i];
    return ret;
    };

int end_swap(const int val, const int size)
    {
    int ret = 0;
    for(int i = 0; i < size; ++i)
        ((char*)&ret)[size-1-i] = ((char*)&val)[i];
    return ret;
    };

void Save_VTK(double* Sx, double* Sy, double* Sz, const int mode, char vtk_filename[64])
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
        snprintf(shortBufer,80,"%s\n",(mode==0 ? "BINARY" : "ASCII" ));
        fputs (shortBufer,pFile);
        fputs ("\n",pFile);
        fputs ("DATASET STRUCTURED_POINTS\n",pFile);
        snprintf(shortBufer,80,"DIMENSIONS %d %d %d \n",uABC[0], uABC[1], uABC[2]);
        fputs (shortBufer,pFile);
        fputs ("ORIGIN 0 0 0\n",pFile);

        temp0  = abc[0][0]*abc[0][0];
        temp0 += abc[0][1]*abc[0][1];
        temp0 += abc[0][2]*abc[0][2];
        temp1 = sqrt(temp0);

        temp0  = abc[1][0]*abc[1][0];
        temp0 += abc[1][1]*abc[1][1];
        temp0 += abc[1][2]*abc[1][2];
        temp2 = sqrt(temp0); 

        temp0  = abc[2][0]*abc[2][0];
        temp0 += abc[2][1]*abc[2][1];
        temp0 += abc[2][2]*abc[2][2];
        temp3 = sqrt(temp0);  

        //metka
        // snprintf(shortBufer,80,"SPACING %.6g %.6g %.6g \n",temp1*a_lattice, temp2*a_lattice, temp3*a_lattice);
        snprintf(shortBufer,80,"SPACING %.6g %.6g %.6g \n",temp1, temp2, temp3);
        fputs (shortBufer,pFile);

        snprintf(shortBufer,80,"POINT_DATA %d \n",uABC[0]*uABC[1]*uABC[2]);
        fputs (shortBufer,pFile);
        fputs ("\n",pFile);
        // //metka popytka dobavit' color v vtk file
        // fputs ("VECTORS spins float\n",pFile);
        // for (int cn = 0; cn<uABC[2]; cn++){
        //     for (int bn = 0; bn<uABC[1]; bn++){
        //         for (int an = 0; an<uABC[0]; an++){
        //             int n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
        //             n = n*AtomsPerBlock;//index of the first spin in the block
        //             temp0 = 0.f;
        //             temp1 = 0.f;
        //             temp2 = 0.f;
        //             for (int atom=0; atom<AtomsPerBlock; atom++){
        //                 int N = n + atom;
        //                 temp0 += Sx[N];
        //                 temp1 += Sy[N];
        //                 temp2 += Sz[N]; 
        //             }
        //             temp0 = temp0/AtomsPerBlock;
        //             temp1 = temp1/AtomsPerBlock;
        //             temp2 = temp2/AtomsPerBlock;
        //             if (mode==0){
        //                 temp0 = end_swap(temp0, sizeof(temp0));
        //                 temp1 = end_swap(temp1, sizeof(temp1));
        //                 temp2 = end_swap(temp2, sizeof(temp2));
        //                 fwrite ((const char*)&temp0, sizeof(float), 1, pFile);
        //                 fwrite ((const char*)&temp1, sizeof(float), 1, pFile);
        //                 fwrite ((const char*)&temp2, sizeof(float), 1, pFile);
        //             }else{
        //                 snprintf(shortBufer,80,"%.6g %.6g %.6g ",temp0, temp1, temp2);
        //                 fputs (shortBufer,pFile);
        //             }
        //         }
        //     }
        // }

        // snprintf(shortBufer,200,"COLOR_SCALARS sccol 3\n");
        // fputs (shortBufer,pFile);
        // float RGB[3];
        // float temp[3];

        // for (int i = 0; i<NOS; i++){
        //     temp[0]=Sx[i];
        //     temp[1]=Sy[i];
        //     temp[2]=Sz[i];
        //     HSVtoRGB(temp, RGB, 1, 0 );
        //     snprintf(shortBufer,200,"%.6g %.6g %.6g \n", RGB[0], RGB[1], RGB[2]);
        //     fputs (shortBufer,pFile);   
        // }
        fputs ("SCALARS m float 3\n",pFile);
        fputs ("LOOKUP_TABLE default\n",pFile);
        for (int cn = 0; cn<uABC[2]; cn++){
            for (int bn = 0; bn<uABC[1]; bn++){
                for (int an = 0; an<uABC[0]; an++){
                    int n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
                    n = n*AtomsPerBlock;//index of the first spin in the block
                    temp0 = 0.f;
                    temp1 = 0.f;
                    temp2 = 0.f;
                    for (int atom=0; atom<AtomsPerBlock; atom++){
                        int N = n + atom;
                        temp0 += Sx[N];
                        temp1 += Sy[N];
                        temp2 += Sz[N]; 
                    }
                    temp0 = temp0/AtomsPerBlock;
                    temp1 = temp1/AtomsPerBlock;
                    temp2 = temp2/AtomsPerBlock;
                    if (mode==0){
                        temp0 = end_swap(temp0, sizeof(temp0));
                        temp1 = end_swap(temp1, sizeof(temp1));
                        temp2 = end_swap(temp2, sizeof(temp2));
                        fwrite ((const char*)&temp0, sizeof(float), 1, pFile);
                        fwrite ((const char*)&temp1, sizeof(float), 1, pFile);
                        fwrite ((const char*)&temp2, sizeof(float), 1, pFile);
                    }else{
                        snprintf(shortBufer,80,"%.6g %.6g %.6g ",temp0, temp1, temp2);
                        fputs (shortBufer,pFile);
                    }
                }
            }
        }
        fputs ("\n",pFile);
        fputs ("SCALARS F1 float 1\n",pFile);
        fputs ("LOOKUP_TABLE default\n",pFile);
        for (int cn = 0; cn<uABC[2]; cn++){
            for (int bn = 0; bn<uABC[1]; bn++){
                for (int an = 0; an<uABC[0]; an++){
                    int n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
                    n = n*AtomsPerBlock;//index of the first spin in the block
                    temp0 = 0.f;
                    temp1 = 0.f;
                    temp2 = 0.f;
                    for (int atom=0; atom<AtomsPerBlock; atom++){
                        int N = n + atom;
                        temp0 += Sx[N];
                        temp1 += Sy[N];
                        temp2 += Sz[N];
                    }
                        temp0 = temp0/AtomsPerBlock;
                        temp1 = temp1/AtomsPerBlock;
                        temp2 = temp2/AtomsPerBlock;
                        temp3 = sqrt(temp0*temp0+temp1*temp1+temp2*temp2);
                        temp0 = temp0/temp3;
                        temp1 = temp1/temp3;
                        temp2 = temp2/temp3;
                        temp0 = atan2(temp1,temp0) * 180 / PI;
                    if (mode==0){
                        temp0 = end_swap(temp0, sizeof(temp0));
                        fwrite ((const char*)&temp0, sizeof(float), 1, pFile);
                    }else{
                        snprintf(shortBufer,80,"%.6g ",temp0);
                        fputs (shortBufer,pFile);
                    }
                }
            }
        }
        fputs ("\n",pFile);
        fputs ("SCALARS F2 float 1\n",pFile);
        fputs ("LOOKUP_TABLE default\n",pFile);
        for (int cn = 0; cn<uABC[2]; cn++){
            for (int bn = 0; bn<uABC[1]; bn++){
                for (int an = 0; an<uABC[0]; an++){
                    int n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
                    n = n*AtomsPerBlock;//index of the first spin in the block
                    temp0 = 0.f;
                    temp1 = 0.f;
                    temp2 = 0.f;
                    for (int atom=0; atom<AtomsPerBlock; atom++){
                        int N = n + atom;
                        temp0 += Sx[N];
                        temp1 += Sy[N];
                        temp2 += Sz[N];
                    }
                        temp0 = temp0/AtomsPerBlock;
                        temp1 = temp1/AtomsPerBlock;
                        temp2 = temp2/AtomsPerBlock;
                        temp3 = sqrt(temp0*temp0+temp1*temp1+temp2*temp2);
                        temp0 = temp0/temp3;
                        temp1 = temp1/temp3;
                        temp2 = temp2/temp3;
                        temp0 = atan2(temp1,-temp0) * 180 / PI;
                    if (mode==0){
                        temp0 = end_swap(temp0, sizeof(temp0));
                        fwrite ((const char*)&temp0, sizeof(float), 1, pFile);
                    }else{
                        snprintf(shortBufer,80,"%.6g ",temp0);
                        fputs (shortBufer,pFile);
                    }
                }
            }
        }
        // fputs ("\n",pFile);
        // fputs ("SCALARS abs(m) float 1\n",pFile);
        // fputs ("LOOKUP_TABLE default\n",pFile);
        // for (int cn = 0; cn<uABC[2]; cn++){
        //     for (int bn = 0; bn<uABC[1]; bn++){
        //         for (int an = 0; an<uABC[0]; an++){
        //             int n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
        //             n = n*AtomsPerBlock;//index of the first spin in the block
        //             temp0 = 0.f;
        //             for (int atom=0; atom<AtomsPerBlock; atom++){
        //                 int N = n + atom;
        //                 temp0 += sqrt(Sx[N]*Sx[N]+Sy[N]*Sy[N]+Sz[N]*Sz[N]);
        //             }
        //             temp0 = temp0/AtomsPerBlock;
        //             if (mode==0){
        //                 temp0 = end_swap(temp0, sizeof(temp0));
        //                 fwrite ((const char*)&temp0, sizeof(float), 1, pFile);
        //             }else{
        //                 snprintf(shortBufer,80,"%.6g ",temp0);
        //                 fputs (shortBufer,pFile);
        //             }
        //         }
        //     }
        // }
        fclose (pFile);
        printf("Recording to the file %s is done!\n", vtk_filename);
    }
}

void Save_VTK_6(double* Sx, double* Sy, double* Sz,
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
        snprintf(shortBufer,80,"%s\n",(mode==0 ? "BINARY" : "ASCII" ));
        fputs (shortBufer,pFile);
        fputs ("\n",pFile);
        fputs ("DATASET STRUCTURED_POINTS\n",pFile);
        snprintf(shortBufer,80,"DIMENSIONS %d %d %d \n",uABC[0], uABC[1], uABC[2]);
        fputs (shortBufer,pFile);
        fputs ("ORIGIN 0 0 0\n",pFile);

        temp0  = abc[0][0]*abc[0][0];
        temp0 += abc[0][1]*abc[0][1];
        temp0 += abc[0][2]*abc[0][2];
        temp1 = sqrt(temp0);

        temp0  = abc[1][0]*abc[1][0];
        temp0 += abc[1][1]*abc[1][1];
        temp0 += abc[1][2]*abc[1][2];
        temp2 = sqrt(temp0); 

        temp0  = abc[2][0]*abc[2][0];
        temp0 += abc[2][1]*abc[2][1];
        temp0 += abc[2][2]*abc[2][2];
        temp3 = sqrt(temp0);  

        snprintf(shortBufer,80,"SPACING %.6g %.6g %.6g \n",temp1*a_lattice, temp2*a_lattice, temp3*a_lattice);
        fputs (shortBufer,pFile);

        snprintf(shortBufer,80,"POINT_DATA %d \n",uABC[0]*uABC[1]*uABC[2]);
        fputs (shortBufer,pFile);
        fputs ("\n",pFile);
        fputs ("SCALARS m float 3\n",pFile);
        fputs ("LOOKUP_TABLE default\n",pFile);
        for (int cn = 0; cn<uABC[2]; cn++){
            for (int bn = 0; bn<uABC[1]; bn++){
                for (int an = 0; an<uABC[0]; an++){
                    int n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
                    n = n*AtomsPerBlock;//index of the first spin in the block
                    temp0 = 0.f;
                    temp1 = 0.f;
                    temp2 = 0.f;
                    for (int atom=0; atom<AtomsPerBlock; atom++){
                        int N = n + atom;
                        temp0 += Sx[N];
                        temp1 += Sy[N];
                        temp2 += Sz[N]; 
                    }
                    temp0 = temp0/AtomsPerBlock;
                    temp1 = temp1/AtomsPerBlock;
                    temp2 = temp2/AtomsPerBlock;
                    if (mode==0){
                        temp0 = end_swap(temp0, sizeof(temp0));
                        temp1 = end_swap(temp1, sizeof(temp1));
                        temp2 = end_swap(temp2, sizeof(temp2));
                        fwrite ((const char*)&temp0, sizeof(float), 1, pFile);
                        fwrite ((const char*)&temp1, sizeof(float), 1, pFile);
                        fwrite ((const char*)&temp2, sizeof(float), 1, pFile);
                    }else{
                        snprintf(shortBufer,80,"%.6g %.6g %.6g ",temp0, temp1, temp2);
                        fputs (shortBufer,pFile);
                    }
                }
            }
        }

        fputs ("\n",pFile);
        fputs ("SCALARS dm float 3\n",pFile);
        fputs ("LOOKUP_TABLE default\n",pFile);
        for (int cn = 0; cn<uABC[2]; cn++){
            for (int bn = 0; bn<uABC[1]; bn++){
                for (int an = 0; an<uABC[0]; an++){
                    int n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
                    n = n*AtomsPerBlock;//index of the first spin in the block
                    temp0 = 0.f;
                    temp1 = 0.f;
                    temp2 = 0.f;
                    for (int atom=0; atom<AtomsPerBlock; atom++){
                        int N = n + atom;
                        temp0 += dSx[N];
                        temp1 += dSy[N];
                        temp2 += dSz[N]; 
                    }
                    temp0 = temp0/AtomsPerBlock;
                    temp1 = temp1/AtomsPerBlock;
                    temp2 = temp2/AtomsPerBlock;
                    if (mode==0){
                        temp0 = end_swap(temp0, sizeof(temp0));
                        temp1 = end_swap(temp1, sizeof(temp1));
                        temp2 = end_swap(temp2, sizeof(temp2));
                        fwrite ((const char*)&temp0, sizeof(float), 1, pFile);
                        fwrite ((const char*)&temp1, sizeof(float), 1, pFile);
                        fwrite ((const char*)&temp2, sizeof(float), 1, pFile);
                    }else{
                        snprintf(shortBufer,80,"%.6g %.6g %.6g ",temp0, temp1, temp2);
                        fputs (shortBufer,pFile);
                    }
                }
            }
        }

        fputs ("\n",pFile);
        fputs ("SCALARS F1 float 1\n",pFile);
        fputs ("LOOKUP_TABLE default\n",pFile);
        for (int cn = 0; cn<uABC[2]; cn++){
            for (int bn = 0; bn<uABC[1]; bn++){
                for (int an = 0; an<uABC[0]; an++){
                    int n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
                    n = n*AtomsPerBlock;//index of the first spin in the block
                    temp0 = 0.f;
                    temp1 = 0.f;
                    temp2 = 0.f;
                    for (int atom=0; atom<AtomsPerBlock; atom++){
                        int N = n + atom;
                        // temp0 += Sx[N];
                        // temp1 += Sy[N];
                        // temp2 += Sz[N];
                        temp0 += dSx[N];
                        temp1 += dSy[N];
                        temp2 += dSz[N];
                    }
                        temp0 = temp0/AtomsPerBlock;
                        temp1 = temp1/AtomsPerBlock;
                        temp2 = temp2/AtomsPerBlock;
                        temp3 = sqrt(temp0*temp0+temp1*temp1+temp2*temp2);
                        temp0 = temp0/temp3;
                        temp1 = temp1/temp3;
                        temp2 = temp2/temp3;
                        temp0 = atan2(temp1,temp0) * 180 / PI;
                    if (mode==0){
                        temp0 = end_swap(temp0, sizeof(temp0));
                        fwrite ((const char*)&temp0, sizeof(float), 1, pFile);
                    }else{
                        snprintf(shortBufer,80,"%.6g ",temp0);
                        fputs (shortBufer,pFile);
                    }
                }
            }
        }
        fputs ("\n",pFile);
        fputs ("SCALARS F2 float 1\n",pFile);
        fputs ("LOOKUP_TABLE default\n",pFile);
        for (int cn = 0; cn<uABC[2]; cn++){
            for (int bn = 0; bn<uABC[1]; bn++){
                for (int an = 0; an<uABC[0]; an++){
                    int n = an+bn*uABC[0]+cn*uABC[0]*uABC[1];// index of the block
                    n = n*AtomsPerBlock;//index of the first spin in the block
                    temp0 = 0.f;
                    temp1 = 0.f;
                    temp2 = 0.f;
                    for (int atom=0; atom<AtomsPerBlock; atom++){
                        int N = n + atom;
                        temp0 += Sx[N];
                        temp1 += Sy[N];
                        temp2 += Sz[N];
                    }
                        temp0 = temp0/AtomsPerBlock;
                        temp1 = temp1/AtomsPerBlock;
                        temp2 = temp2/AtomsPerBlock;
                        temp3 = sqrt(temp0*temp0+temp1*temp1+temp2*temp2);
                        temp0 = temp0/temp3;
                        temp1 = temp1/temp3;
                        temp2 = temp2/temp3;
                        temp0 = atan2(temp1,-temp0) * 180 / PI;
                    if (mode==0){
                        temp0 = end_swap(temp0, sizeof(temp0));
                        fwrite ((const char*)&temp0, sizeof(float), 1, pFile);
                    }else{
                        snprintf(shortBufer,80,"%.6g ",temp0);
                        fputs (shortBufer,pFile);
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
    
void Read_VTK(double* Sx, double* Sy, double* Sz, char vtk_filename[64]){
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

                printf("%s has wrong header or wrong file format! (not vtk 2.0 file format):\n", inputfilename);
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
                            if (k<uABC[2] && j<uABC[1] && i<uABC[0]){
                                n = i + j*xnodes + k*xnodes*ynodes; //index of the block!
                                //printf("n=%d\n", n);
                                if (binType==4){
                                    if(!fread(&temp4_x,binType,1,FilePointer)) break;
                                    if(!fread(&temp4_y,binType,1,FilePointer)) break;
                                    if(!fread(&temp4_z,binType,1,FilePointer)) break;
                                    for (int t=0; t<AtomsPerBlock; t++){
                                        int I=n*AtomsPerBlock+t;
                                        temp4_x = end_swap(temp4_x, sizeof(temp4_x));
                                        temp4_y = end_swap(temp4_y, sizeof(temp4_y));
                                        temp4_z = end_swap(temp4_z, sizeof(temp4_z));
                                        Sx[I]=bSx[I]=(double)temp4_x; 
                                        Sy[I]=bSy[I]=(double)temp4_y;
                                        Sz[I]=bSz[I]=(double)temp4_z;      
                                    }
                                }else if (binType==8){
                                    if(!fread(&temp8_x,binType,1,FilePointer)) break;
                                    if(!fread(&temp8_y,binType,1,FilePointer)) break;
                                    if(!fread(&temp8_z,binType,1,FilePointer)) break;
                                    for (int t=0; t<AtomsPerBlock; t++){
                                        int I=n*AtomsPerBlock+t;
                                        temp8_x = end_swap(temp8_x, sizeof(temp8_x));
                                        temp8_y = end_swap(temp8_y, sizeof(temp8_y));
                                        temp8_z = end_swap(temp8_z, sizeof(temp8_z));
                                        Sx[I]=bSx[I]=temp8_x; 
                                        Sy[I]=bSy[I]=temp8_y;
                                        Sz[I]=bSz[I]=temp8_z;       
                                    }
                                }
                            }   
                        }
                    }
                }
            }
        }else{
            printf("%s cannot read data in vtk file!\n", inputfilename);
        }         
        fclose(FilePointer);
        printf("Done!\n");// when everything is done
    }else{printf("Cannot open file: %s \n", inputfilename);}
}
