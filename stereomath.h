#ifndef STEREOMATH_H
#define STEREOMATH_H

#include <math.h>


static inline unsigned char getBit(unsigned char byte, int position) // position in range 0-7
{
    return (byte >> position) & 0x01;
}

static inline unsigned char setBit(unsigned char x, unsigned char k, unsigned char b) 
{// x: 8-bit value, bit position in range 0-7,  b: set bit either to 1 or 0
   return (b ?  (x | (0x01 << k))  :  (x & ~(0x01 << k)) );
              //   Set bit to 1           Set bit to 0
}

static inline void gamma_to_xyz(double g1, double g2, unsigned char T, double x, double y, double z)
{ 
	double etta = 2.0f/(1.0f + g1*g1 + g2*g2);
	x = g1*etta;
	y = g2*etta;
	//z = ( getBit(T, 0) ? (etta - 1.0f) : (1.0f - etta) );
	//                     if T = 1        if T = 0
	//                    South pole      North pole
	z = 2 * (getBit(T, 0)-0.5);
}

static inline void xyz_to_gammaS(double x, double y, double z, double g1, double g2, unsigned char T)
{
	double temp = 1.0f/(1.0f + z);
	g1 = x * temp;
	g2 = y * temp;
	T = setBit(T, 0, 1);
}

static inline void xyz_to_gammaN(double x, double y, double z, double g1, double g2, unsigned char T)
{
	double temp = 1.0f/(1.0f + z);
	g1 = x * temp;
	g2 = y * temp;
	T = setBit(T, 0, 0);
}

static inline void xyz_to_gamma(double x, double y, double z, double g1, double g2, unsigned char T, bool b)
{// b true (false) -> adapt (keep) pole projection
	
	if (getBit(T, 0)){//South pole
		if (b) {//adapt projection
			if (z<-0.6){//cahnge projection to N pole
				xyz_to_gammaN(x,y,z,g1,g2,T);
			}else{//keep projection from S pole
				xyz_to_gammaS(x,y,z,g1,g2,T);
			}
		}else{//keep projection from S pole
			xyz_to_gammaS(x,y,z,g1,g2,T);
		}
	}else{//North pole
		if (b) {//adapt projection
			if (z>0.6){//cahnge projection to S pole
				xyz_to_gammaS(x,y,z,g1,g2,T);
			}else{//keep projection from N pole
				xyz_to_gammaN(x,y,z,g1,g2,T);
			}
		}else{//keep projection from N pole
			xyz_to_gammaN(x,y,z,g1,g2,T);
		}
	}
} 

#endif
