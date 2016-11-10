inline int sgn(float x) { return (x >= 0) - (x < 0);}

double 
Dot(double* V0, double* V1)
	{
	return (V0[0]*V1[0] + V0[1]*V1[1] + V0[2]*V1[2]);
	};

float  Dotf(float* V0, float* V1)
	{
	return (V0[0]*V1[0] + V0[1]*V1[1] + V0[2]*V1[2]);
	};

void 
Cross(double* V0, double* V1, double* R)
	{
	R[0] = (V0[1]*V1[2] - V0[2]*V1[1]);
	R[1] = (V0[2]*V1[0] - V0[0]*V1[2]);
	R[2] = (V0[0]*V1[1] - V0[1]*V1[0]);
	};

void 
Crossf(float* V0, float* V1, float* R)
	{
	R[0] = (V0[1]*V1[2] - V0[2]*V1[1]);
	R[1] = (V0[2]*V1[0] - V0[0]*V1[2]);
	R[2] = (V0[0]*V1[1] - V0[1]*V1[0]);
	};

double
Unit( double vin[3], double vout[3] )
{
	double dist = vin[0]*vin[0] + vin[1]*vin[1] + vin[2]*vin[2];

	if( dist > 0.0 )
	{
		dist = sqrt( dist );
		vout[0] = vin[0] / dist;
		vout[1] = vin[1] / dist;
		vout[2] = vin[2] / dist;
	}
	else
	{
		vout[0] = vin[0];
		vout[1] = vin[1];
		vout[2] = vin[2];
	}

	return dist;
}

float
Unitf( float vin[3], float vout[3] )
{
	float dist = vin[0]*vin[0] + vin[1]*vin[1] + vin[2]*vin[2];

	if( dist > 0.0 )
	{
		dist = sqrt( dist );
		vout[0] = vin[0] / dist;
		vout[1] = vin[1] / dist;
		vout[2] = vin[2] / dist;
	}
	else
	{
		vout[0] = vin[0];
		vout[1] = vin[1];
		vout[2] = vin[2];
	}

	return dist;
}

void
Enorm( float v1[3], float v2[3], float v3[3], float vout[3] )
{// normal vector to the element based on three vertexes v1, v2, v3	
	float tmp1[3];
	float tmp2[3];
	tmp1[0]=v2[0]-v1[0];
	tmp1[1]=v2[1]-v1[1];
	tmp1[2]=v2[2]-v1[2];

	tmp2[0]=v3[0]-v1[0];
	tmp2[1]=v3[1]-v1[1];
	tmp2[2]=v3[2]-v1[2];
	//Crossf( tmp1, tmp2, tmp3 );
	vout[0] = (tmp1[1]*tmp2[2] - tmp1[2]*tmp2[1]);
	vout[1] = (tmp1[2]*tmp2[0] - tmp1[0]*tmp2[2]);
	vout[2] = (tmp1[0]*tmp2[1] - tmp1[1]*tmp2[0]);
	(void)Unitf( vout, vout);
}

void
Enorm2( float * v1, float * v2, float * v3, float * vout)
{// normal vector to the element based on three vertexes v1, v2, v3	
	float tmp1[3];
	float tmp2[3];
	tmp1[0]=v2[0]-v1[0];
	tmp1[1]=v2[1]-v1[1];
	tmp1[2]=v2[2]-v1[2];

	tmp2[0]=v3[0]-v1[0];
	tmp2[1]=v3[1]-v1[1];
	tmp2[2]=v3[2]-v1[2];
	//Crossf( tmp1, tmp2, tmp3 );
	vout[0] = (tmp1[1]*tmp2[2] - tmp1[2]*tmp2[1]);
	vout[1] = (tmp1[2]*tmp2[0] - tmp1[0]*tmp2[2]);
	vout[2] = (tmp1[0]*tmp2[1] - tmp1[1]*tmp2[0]);
	(void)Unitf( vout, vout);
}

void
RotateVector(float v1x, float v1y, float v1z, float v2x, float v2y, float v2z, float T, double vout[3])
{
// Rotating the vector <v1x,v1y,v1z> about the axis <v2x,v2y,v2z> by the angle T in degree	
//                     <  x,  y,  z>                  <u,  v,  w>              θ
// | x'| | [u(ux+vy+wz)(1-cosθ) + U*x*cosθ + sqrt(U)*(-wy+vz)*sinθ]/U |
// | y'|=| [v(ux+vy+wz)(1-cosθ) + U*y*cosθ + sqrt(U)*( wx-uz)*sinθ]/U |
// | z'| | [w(ux+vy+wz)(1-cosθ) + U*z*cosθ + sqrt(U)*(-vx+uy)*sinθ]/U |
// where U=u^2+v^2+w^2
	float U = v2x*v2x + v2y*v2y + v2z*v2z;
	float iU= 1.0/U;
	float sU= sqrt(U);
	float cosT = cos(PI*T/180);
	float sinT = sin(PI*T/180);
	//
	float ux = v2x * v1x;
	float vy = v2y * v1y;
	float wz = v2z * v1z;
	//
	float wy = v2z * v1y;
	float uy = v2x * v1y;
	float vz = v2y * v1z;
	float uz = v2x * v1z;
	float vx = v2y * v1x;
	float wx = v2z * v1x;
	//output:
	vout[0] = (v2x* (ux+vy+wz)*(1-cosT) + U*v1x*cosT + sU*(-wy+vz) *sinT) * iU ;
	vout[1] = (v2y* (ux+vy+wz)*(1-cosT) + U*v1y*cosT + sU*( wx-uz) *sinT) * iU ;
	vout[2] = (v2z* (ux+vy+wz)*(1-cosT) + U*v1z*cosT + sU*(-vx+uy) *sinT) * iU ;		
}

void
NewBasisCartesian(float vin1[3], float vin2[3], float vout[3] )
{
// Rotating the point (x,y,z) about the axis ⟨u,v,w⟩ by the angle θ	
// | x'| | [u(ux+vy+wz)(1-cosθ) + U*x*cosθ + sqrt(U)*(-wy+vz)*sinθ]/U |
// | y'|=| [v(ux+vy+wz)(1-cosθ) + U*y*cosθ + sqrt(U)*( wx-uz)*sinθ]/U |
// | z'| | [w(ux+vy+wz)(1-cosθ) + U*z*cosθ + sqrt(U)*(-vx+uy)*sinθ]/U |
// where U=u^2+v^2+w^2
// Let define vector (X,Y,Z) towards which we rotate our point (x,y,z).
// Here we assume the rotation axis is orthogonal
// to the projection of vector (X,Y,Z) on xy-plane
// meaning u=-Y, v=+X, w=0 =>
// U = u^2+v^2 = Y*Y+X*X and => sinθ=sqrt(U), cosθ=Z
// sqrt(U)*sinθ=U 
//  (ux+vy+wz) = (-Y*x+X*y);
// -wy+vz = X*z;
//  wx-uz = Y*z;
// -vx+uy =-X*x-Y*y;
// | x'| |[-Y*(-Y*x+X*y)(1-Z) + U*x*Z + sqrt(U)*(  X*z  )*sqrt(U)]/U |   |[Y*(Y*x-X*y)(1-Z)]/U + x*Z + (  X*z  ) |
// | y'|=|[ X*(-Y*x+X*y)(1-Z) + U*y*Z + sqrt(U)*(  Y*z  )*sqrt(U)]/U | = |[Y*(Y*x-X*y)(1-Z)]/U + y*Z + (  Y*z  ) |
// | z'| |[                     U*z*Z - sqrt(U)*(X*x+Y*y)*sqrt(U)]/U |   |[Y*(Y*x-X*y)(1-Z)]/U + z*Z - (X*x+Y*y) |
// final equation
// | x'|  |-Y*A + x*Z + X*z|  
// | y'|= | X*A + y*Z + Y*z|, U=Y*Y+X*X, A=(-Y*x+X*y)(1-Z)]/U
// | z'|  |       z*Z -  U |
	float x = vin1[0]; float X = vin2[0];
	float y = vin1[1]; float Y = vin2[1];
	float z = vin1[2]; float Z = vin2[2];
	float U = Y*Y+X*X+(1e-37f); //<--minimal float is added to avoid division by zero
	float A = (-Y*x + X*y)*(1. - Z)/U; 
	vout[0] =-Y*A + x*Z + X*z;
	vout[1] = X*A + y*Z + Y*z;
	vout[2] = z*Z - (X*x+Y*y);		
}

float atan2int( float y, float x )
{
    static const uint32_t sign_mask = 0x80000000;

    // Extract the sign bits
    uint32_t ux_s  = sign_mask & (uint32_t &)x;
    uint32_t uy_s  = sign_mask & (uint32_t &)y;

    // Determine the quadrant offset
    float q = (float)( ( ~ux_s & uy_s ) >> 29 | ux_s >> 30 ); 

    // Calculate the arctangent in the first quadrant
    float bxy_a = fabs( 0.596227f * x * y );
    float num = bxy_a + y * y;
    float atan_1q =  num / ( x * x + bxy_a + num );

    // Translate it to the proper quadrant
    uint32_t uatan_2q = (ux_s ^ uy_s) | (uint32_t &)atan_1q;
    return (q + (float &)uatan_2q)*90.f; // Pi/2*180/Pi=90;
}