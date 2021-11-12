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
	float tmp1[3]={v2[0]-v1[0],v2[1]-v1[1],v2[2]-v1[2]};
	float tmp2[3]={v3[0]-v1[0],v3[1]-v1[1],v3[2]-v1[2]};
	// tmp1[0]=v2[0]-v1[0];
	// tmp1[1]=v2[1]-v1[1];
	// tmp1[2]=v2[2]-v1[2];

	// tmp2[0]=v3[0]-v1[0];
	// tmp2[1]=v3[1]-v1[1];
	// tmp2[2]=v3[2]-v1[2];
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


// Routine to set a quaternion from a rotation axis and angle
// ( input axis = float[3] angle = float  output: quat = float[4] )
void SetQuaternionFromAxisAngle(float *axis, float angle, float *quat)
{
    float sina2, norm;
    sina2 = (float)sin(0.5f * angle);
    norm = (float)sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
    quat[0] = sina2 * axis[0] / norm;
    quat[1] = sina2 * axis[1] / norm;
    quat[2] = sina2 * axis[2] / norm;
    quat[3] = (float)cos(0.5f * angle);
}

void SetQuaternionFromVector(float *vector, float *quat)
{
    quat[0] = vector[0];//i
    quat[1] = vector[1];//j
    quat[2] = vector[2];//k
    quat[3] = 0;//Re
}

void GetVectorFromQuaternion(float *quat, float *vector)
{
    vector[0] = quat[0];
    vector[1] = quat[1];
    vector[2] = quat[2];
}

void GetEulerFromQuaternion(float *quat, float *vector)
{
    vector[0] = R2D*atan2(2*quat[0]*quat[3]-2*quat[1]*quat[2] , 1 - 2*quat[0]*quat[0] - 2*quat[2]*quat[2]);
    // bank = atan2(2*qx*qw-2*qy*qz , 1 - 2*qx^2 - 2*qz^2)
    vector[1] = R2D*atan2(2*quat[1]*quat[3]-2*quat[0]*quat[2] , 1 - 2*quat[1]*quat[1] - 2*quat[2]*quat[2]);
    //heading = atan2(2*qy*qw                      -2*qx*qz , 1 - 2*qy^2 - 2*qz^2)
    vector[2] = R2D*asin(2*quat[0]*quat[1]+2*quat[2]*quat[3]);
}

void GetQuaternionFromEuler(float *q, float *Rot)
{
    const float D2R2=0.008726646255;//0.5*Pi/180 degrees to radians
    // Assuming the angles are in radians.
    float c1 = cos(D2R2 * Rot[1]);//heading
    float s1 = sin(D2R2 * Rot[1]);//heading
    float c2 = cos(D2R2 * Rot[2]);//attitude
    float s2 = sin(D2R2 * Rot[2]);//attitude
    float c3 = cos(D2R2 * Rot[0]);//bank
    float s3 = sin(D2R2 * Rot[0]);//bank
    float c1c2 = c1*c2;
    float s1s2 = s1*s2;
    q[3] = c1c2*c3 - s1s2*s3;
    q[0] = c1c2*s3 + s1s2*c3;
    q[1] = s1*c2*c3 + c1*s2*s3;
    q[2] = c1*s2*c3 - s1*c2*s3;
}
//metka
void GetNormQuaternion(float *q, float *qout)
{
float magnitude = sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
qout[0] = q[0] / magnitude;
qout[1] = q[1] / magnitude;
qout[2] = q[2] / magnitude;
qout[3] = q[3] / magnitude;
}

void GetQuaternionFromTwoVectors(float *q, float *V0, float *V1)
{
        float ex[3]={1,0,0};
        float ey[3]={0,1,0};
        float vt[3]={0,0,0};
        float dot = Dotf(V0, V1);
        if (dot < -0.999999) {
            Crossf(V0, ex, vt);
            if (vt[0]*vt[0]+vt[1]*vt[1]+vt[2]*vt[2] < 1e-12)
                Crossf(V0, ey, vt);
            Unitf( vt, vt );
            SetQuaternionFromAxisAngle(vt, PI, q);
        } else if (dot > 0.999999) {
            q[0] = 0;
            q[1] = 0;
            q[2] = 0;
            q[3] = 1;
        } else {
            Crossf(V0, V1, vt);
            q[0] = vt[0];
            q[1] = vt[1];
            q[2] = vt[2];
            q[3] = 1 + dot;
            GetNormQuaternion(q, q);
        }
}

// Routine to convert a quaternion to a 4x4 matrix
// ( input: quat = float[4]  output: mat = float[4*4] )
void ConvertQuaternionToMatrix(float *quat, float *mat)
{
    float yy2 = 2.0f * quat[1] * quat[1];
    float xy2 = 2.0f * quat[0] * quat[1];
    float xz2 = 2.0f * quat[0] * quat[2];
    float yz2 = 2.0f * quat[1] * quat[2];
    float zz2 = 2.0f * quat[2] * quat[2];
    float wz2 = 2.0f * quat[3] * quat[2];
    float wy2 = 2.0f * quat[3] * quat[1];
    float wx2 = 2.0f * quat[3] * quat[0];
    float xx2 = 2.0f * quat[0] * quat[0];
    mat[0*4+0] = - yy2 - zz2 + 1.0f;
    mat[0*4+1] = xy2 + wz2;
    mat[0*4+2] = xz2 - wy2;
    mat[0*4+3] = 0;
    mat[1*4+0] = xy2 - wz2;
    mat[1*4+1] = - xx2 - zz2 + 1.0f;
    mat[1*4+2] = yz2 + wx2;
    mat[1*4+3] = 0;
    mat[2*4+0] = xz2 + wy2;
    mat[2*4+1] = yz2 - wx2;
    mat[2*4+2] = - xx2 - yy2 + 1.0f;
    mat[2*4+3] = 0;
    mat[3*4+0] = mat[3*4+1] = mat[3*4+2] = 0;
    mat[3*4+3] = 1;
}

void RotateVectorByQuaternion(float *v, float *q)
{
    float x=v[0];
    float y=v[1];
    float z=v[2];
    //
    float qxqx=q[0]*q[0];
    float qyqy=q[1]*q[1];
    float qzqz=q[2]*q[2];
    float qwqw=q[3]*q[3];
    //
    float qxqy=q[0]*q[1];
    float qxqz=q[0]*q[2];
    float qyqz=q[1]*q[2];
    float qwqx=q[3]*q[0];
    float qwqy=q[3]*q[1];
    float qwqz=q[3]*q[2];

    v[0] = x*(   qxqx +   qwqw-qyqy- qzqz) + y*(2*qxqy- 2*qwqz) + z*(2*qxqz+ 2*qwqy);
    v[1] = x*( 2*qwqz + 2*qxqy) + y*(  qwqw - qxqx + qyqy - qzqz)+ z*(-2*qwqx+ 2*qyqz);
    v[2] = x*(-2*qwqy + 2*qxqz) + y*(2*qwqx + 2*qyqz) + z*(qwqw - qxqx- qyqy+ qzqz);
}
// Routine to multiply 2 quaternions (ie, compose rotations)
// ( input q1 = float[4] q2 = float[4]  output: qout = float[4] )
void MultiplyQuaternions(float *q1, float *q2, float *qout)
{
    float qr[4];
    qr[0] = q1[3]*q2[0] + q1[0]*q2[3] + q1[1]*q2[2] - q1[2]*q2[1];
    qr[1] = q1[3]*q2[1] + q1[1]*q2[3] + q1[2]*q2[0] - q1[0]*q2[2];
    qr[2] = q1[3]*q2[2] + q1[2]*q2[3] + q1[0]*q2[1] - q1[1]*q2[0];
    qr[3] = q1[3]*q2[3] -(q1[0]*q2[0] + q1[1]*q2[1] + q1[2]*q2[2]);
    qout[0] = qr[0]; qout[1] = qr[1]; qout[2] = qr[2]; qout[3] = qr[3];
}
