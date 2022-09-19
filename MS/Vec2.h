//
//      Vec2<TPrecision>.h        Basic vector 2D class (TPrecision precision)
//
#pragma once

#ifndef _VEC2_H
#define _VEC2_H

#include <cmath>

namespace VEC2
{
	// CF_FAST_MATH implies imprecision in the methods:
	// -  "normalize()"
	// -  "length()"
	// - operator "=="
	// - operator "!="
	//#define CF_FAST_MATH

	// 32 bits = 6 decimal digits
	// 64 bits = 15 decimal digits
#ifdef CF_FAST_MATH
		//from irrMath.h:
		const float  ROUNDING_ERROR_32	= 0.00005f;
		const double ROUNDING_ERROR_64	= 0.000005f;
#else

		// Best possible precision
		//const float  ROUNDING_ERROR_32	= 0.0000001f;
		//const double ROUNDING_ERROR_64	= 0.0000000000000001;

		//used in irrMath.h (irrlitch): 
		const static float  ROUNDING_ERROR_32	= 0.000001f;
		const static double ROUNDING_ERROR_64	= 0.00000001f;

#endif

	//! Constant for PI.
	const float PI = 3.14159265359f;
	//! 32bit Constant for converting from degrees to radians
	const float DEGTORAD = PI / 180.0f;
	//! 32bit constant for converting from radians to degrees (formally known as GRAD_PI)
	const float RADTODEG   = 180.0f / PI;

	//! Constant for 64bit PI.
	const static double PI64 = 3.1415926535897932384626433832795028841971693993751;
	//! 64bit constant for converting from degrees to radians (formally known as GRAD_PI2)
	const static double DEGTORAD64 = PI64 / 180.0;
	//! 64bit constant for converting from radians to degrees
	const static double RADTODEG64 = 180.0 / PI64;
};

inline bool isNAN(const float value) { return value!=value; };

//-----------------------------------------------------------------------------
template <class TPrecision>
class Vec2
{
public:

	// Constructor that sets all components of the vector to zero
    Vec2() { _v[0]=0.0; _v[1]=0.0; }
	
	Vec2(const TPrecision& x, const TPrecision& y) { _v[0] = x; _v[1] = y; }

	inline Vec2(const Vec2& vec) { _v[0] = vec._v[0]; _v[1] = vec._v[1]; }

	~Vec2() {};

	// To be used with OpenGL functions
	inline TPrecision* ptr() { return _v; }
    inline const TPrecision* ptr() const { return _v; }

    inline void set(const TPrecision& x, const TPrecision& y) { _v[0]=x; _v[1]=y; }

    inline void set( const Vec2<TPrecision>& rhs) { _v[0] = rhs._v[0]; _v[1] = rhs._v[1]; }

    inline TPrecision& operator [] (int i) { return _v[i]; }
    inline TPrecision  operator [] (int i) const { return _v[i]; }

    inline TPrecision& x() { return _v[0]; }
    inline TPrecision& y() { return _v[1]; }

    inline TPrecision x() const { return _v[0]; }
    inline TPrecision y() const { return _v[1]; }

	// Dot product.
    inline TPrecision operator * (const Vec2<TPrecision>& rhs) const
    {
        return _v[0]*rhs._v[0] + _v[1]*rhs._v[1];
    }

    // Multiply by scalar. 
    inline const Vec2<TPrecision> operator * (const TPrecision rhs) const
    {
        return Vec2<TPrecision>(_v[0]*rhs, _v[1]*rhs);
    }

    // Unary multiply by scalar. 
    inline Vec2<TPrecision>& operator *= (const TPrecision rhs)
    {
        _v[0]*=rhs;
        _v[1]*=rhs;
        return *this;
    }

	// Divide by scalar. 
    inline const Vec2<TPrecision> operator / (const TPrecision rhs) const
    {
        return Vec2<TPrecision>(_v[0]/rhs, _v[1]/rhs);
    }

    // Unary divide by scalar.
    inline Vec2<TPrecision>& operator /= (const TPrecision rhs)
    {
        _v[0]/=rhs;
        _v[1]/=rhs;
        return *this;
    }

	 // Binary vector add. 
    inline const Vec2<TPrecision> operator + (const Vec2<TPrecision>& rhs) const
    {
        return Vec2<TPrecision>(_v[0]+rhs._v[0], _v[1]+rhs._v[1]);
    }

    // Unary vector add. Slightly more efficient because no temporary
    // intermediate object.
    inline Vec2<TPrecision>& operator += (const Vec2<TPrecision>& rhs)
    {
        _v[0] += rhs._v[0];
        _v[1] += rhs._v[1];
        return *this;
    }

    // Binary vector subtract. 
    inline const Vec2<TPrecision> operator - (const Vec2<TPrecision>& rhs) const
    {
        return Vec2<TPrecision>(_v[0]-rhs._v[0], _v[1]-rhs._v[1]);
    }

    // Unary vector subtract.
    inline Vec2<TPrecision>& operator -= (const Vec2<TPrecision>& rhs)
    {
        _v[0]-=rhs._v[0];
        _v[1]-=rhs._v[1];
        return *this;
    }

    // Negation operator. Returns the negative of the Vec2.
    inline const Vec2<TPrecision> operator - () const
    {
        return Vec2<TPrecision>(-_v[0], -_v[1]);
    }

	//inline bool operator == (const Vec2<TPrecision>& v) const { return (_v[0]==v._v[0]) && (_v[1]==v._v[1]) }
	//inline bool operator != (const Vec2<TPrecision>& v) const { return _v[0]!=v._v[0] || _v[1]!=v._v[1] }

	// The methods below are slightly faster than the built-in operator
	// 32 bits
	bool operator==(const Vec2<float>& other) const
	{
		return this->equals(other);
	}

	bool operator!=(const Vec2<float>& other) const
	{
		return !this->equals(other);
	}

	// 64 bits
	bool operator==(const Vec2<double>& other) const
	{
		return this->equals(other);
	}

	bool operator!=(const Vec2<double>& other) const
	{
		return !this->equals(other);
	}

	// It's just a definition
	inline bool operator <  (const Vec2<TPrecision>& v) const
    {
        if (_v[0]<v._v[0]) return true;
        else if (_v[0]>v._v[0]) return false;
        else return (_v[1]<v._v[1]);
    }

	// Length of the vector = sqrt( vec . vec ) 
	// Anyway, float has enough precision to length()
    inline TPrecision length() const
    {
        return sqrt(_v[0]*_v[0] + _v[1]*_v[1]);
    }

    // Length squared of the vector = vec . vec
    inline TPrecision length2() const
    {
        return _v[0]*_v[0] + _v[1]*_v[1];
    }

	//! Normalizes the vector. In case of the 0 vector the result
		//! is still 0, otherwise the length of the vector will be 1.
		//! Todo: 64 Bit template doesnt work.. need specialized template
	inline void normalize()
	{
		TPrecision norm2 = _v[0]*_v[0] + _v[1]*_v[1];
        if (norm2 > 0.0)
        {
			norm2 = inverse_sqrtf( static_cast<float>(norm2) );
			//norm2 = 1.0 / sqrt( norm2 );
            _v[0] *= norm2;
            _v[1] *= norm2;
        }                
	}

	//! Normalizes the vector. In case of the 0 vector the result
		//! is still 0, otherwise the length of the vector will be 1.
		//! Todo: 64 Bit template doesnt work.. need specialized template
	/*inline void normalize()
	{
		TPrecision norm2 = _v[0]*_v[0] + _v[1]*_v[1];
        if (norm2 > 0.0)
        {
			norm2 = 1.0/sqrt(norm2);
            _v[0] *= norm2;
            _v[1] *= norm2;
        }                
	}*/

	// multiply by vector components.
	inline void componentMultiply(const Vec2<TPrecision>& rhs)
	{
		_v[0]*=rhs[0];
		_v[1]*=rhs[1];
	}

	// divide rhs components by rhs vector components.
	inline void componentDivide(const Vec2<TPrecision>& rhs)
	{
		_v[0]/=rhs[0];
		_v[1]/=rhs[1];
	}	

	//! Returns distance from another point.
	/// Here, the vector is interpreted as point in 3 dimensional space.
	TPrecision distanceFrom(const Vec2<TPrecision>& other) const
	{
		return Vec2<TPrecision>(_v[0] - other._v[0], _v[1] - other._v[1]).length();
	}

	//! Returns squared distance from another point.
	// Here, the vector is interpreted as point in 3 dimensional space.
	TPrecision distance2From(const Vec2<TPrecision>& other) const
	{
		return Vec2<TPrecision>(_v[0] - other._v[0], _v[1] - other._v[1]).length2();
	}

	//! Returns if this vector interpreted as a point is on a line between two other points.
	// WARNING: It is assumed that the point is on the line.
	//! \param begin: Beginning vector to compare between.
	//! \param end: Ending vector to compare between.
	//! \return True if this vector is between begin and end.  False if not.
	bool isBetweenPoints(const Vec2<TPrecision>& begin, const Vec2<TPrecision>& end) const
	{
		 if (begin._v[0] != end._v[0])
		 {
				 return ((begin._v[0] <= _v[0] && _v[0] <= end._v[0]) ||
						 (begin._v[0] >= _v[0] && _v[0] >= end._v[0]));
		 }
		 else
		 {
				 return ((begin._v[1] <= _v[1] && _v[1] <= end._v[1]) ||
						 (begin._v[1] >= _v[1] && _v[1] >= end._v[1]));
		 }
	}

	Vec2<TPrecision>& rotateBy(TPrecision degrees, const Vec2<TPrecision>& center=Vec2<TPrecision>())
	{
		degrees *= static_cast<TPrecision>(VEC2::DEGTORAD);
		 const TPrecision cs = cos(degrees);
		 const TPrecision sn = sin(degrees);

		 _v[0] -= center._v[0];
		 _v[1] -= center._v[1];

		 set((_v[0]*cs - _v[1]*sn), (_v[0]*sn + _v[1]*cs));

		 _v[0] += center._v[0];
		 _v[1] += center._v[1];
		 return *this;
	}

	Vec2<TPrecision> getInterpolated(const Vec2<TPrecision>& other, const TPrecision d) const
	{
		 const TPrecision inv = 1.0f - d;
		 return Vec2<TPrecision>(other._v[0]*inv + _v[0]*d, other._v[1]*inv + _v[1]*d);
	}


	Vec2<TPrecision> getInterpolatedQuadratic(const Vec2<TPrecision>& v2, const Vec2<TPrecision>& v3, const TPrecision d) const
	{
		 // this*(1-d)*(1-d) + 2 * v2 * (1-d) + v3 * d * d;
		 const TPrecision inv = 1.0f - d;
		 const TPrecision mul0 = inv * inv;
		 const TPrecision mul1 = 2.0f * d * inv;
		 const TPrecision mul2 = d * d;

		 return Vec2<TPrecision> (_v[0] * mul0 + v2._v[0] * mul1 + v3._v[0] * mul2,
								 _v[1] * mul0 + v2._v[1] * mul1 + v3._v[1] * mul2);
	}


	Vec2<TPrecision>& interpolate(const Vec2<TPrecision>& a, const Vec2<TPrecision>& b, const TPrecision d)
	{
		 _v[0] = b._v[0] + ( ( a._v[0] - b._v[0] ) * d );
		 _v[1] = b._v[1] + ( ( a._v[1] - b._v[1] ) * d );
		 return *this;
	}

	//	Retorna o menor angulo entre este vetor e v. 
	//	O angulo retornado esta em radianos.
	inline TPrecision getMinAngleWithRAD(const Vec2& v) const
	{
		return acos( ((*this) * v) / (length() * v.length()) );
	}

	//	Retorna o menor angulo entre este vetor e v. 
	//	O angulo retornado esta em radianos.
	inline TPrecision getMinAngleWithDEG(const Vec2& v) const
	{
		return acos( ((*this) * v) / (length() * v.length()) )*static_cast<TPrecision>(VEC2::RADTODEG);
	}
	
	// Calculates the angle of this vector in degrees in the trigonometric sense.
	// 0 is to the right (3 o'clock), values increase counter-clockwise. This method has been suggested by Pr3t3nd3r.
	TPrecision getAngleTrigDEG() const
	{
		 if (_v[1] == 0.0f)
				 return _v[0] < 0.0f ? 180.0f : 0.0f;
		 else
		 if (_v[0] == 0.0f)
				 return _v[1] < 0.0f ? 270.0f : 90.0f;

		 if ( _v[1] > 0.0f)
				 if (_v[0] > 0.0f)
						 return atan(_v[1]/_v[0]) * static_cast<TPrecision>(VEC2::RADTODEG64);
				 else
						 return 180.0f - atan(_v[1]/-_v[0]) * static_cast<TPrecision>(VEC2::RADTODEG64);
		 else
				 if (_v[0] > 0.0f)
						 return 360.0f - atan(-_v[1]/_v[0]) * static_cast<TPrecision>(VEC2::RADTODEG64);
				 else
						 return 180.0f + atan(-_v[1]/-_v[0]) * static_cast<TPrecision>(VEC2::RADTODEG64);
	}

	TPrecision getAngleTrigRAD() const
	{
		 if (_v[1] == 0.0f)
				 return _v[0] < 0.0f ? static_cast<TPrecision>(VEC2::PI64) : 0.0f;
		 else
		 if (_v[0] == 0.0f)
				 return _v[1] < 0.0f ? 3.0f*static_cast<TPrecision>(VEC2::PI64)/2.0f : static_cast<TPrecision>(VEC2::PI64)/2.0f;

		 if ( _v[1] > 0.0f)
				 if (_v[0] > 0.0f)
						 return atan(_v[1]/_v[0]);
				 else
						 return static_cast<TPrecision>(VEC2::PI64) - atan(_v[1]/-_v[0]);
		 else
				 if (_v[0] > 0.0f)
						 return 2.0f*static_cast<TPrecision>(VEC2::PI64) - atan(-_v[1]/_v[0]);
				 else
						 return static_cast<TPrecision>(VEC2::PI64) + atan(-_v[1]/-_v[0]);
	}

	// Calculates the angle of this vector in degrees in the counter trigonometric sense.
	//0 is to the right (3 o'clock), values increase clockwise.
	inline TPrecision getAngleDEG() const
	{
		 if (_v[1] == 0.0f) // corrected thanks to a suggestion by Jox
				 return _v[0] < 0.0f ? 180.0f : 0.0f;
		 else if (_v[0] == 0.0f)
				 return _v[1] < 0.0f ? 90.0f : 270.0f;

		 TPrecision tmp = _v[1] / sqrt((_v[0]*_v[0] + _v[1]*_v[1]));
		 tmp = atan( sqrt(1 - tmp*tmp) / tmp) * static_cast<TPrecision>(VEC2::RADTODEG64);

		 if (_v[0]>0.0f && _v[1]>0.0f)
				 return tmp + 270.0f;
		 else
		 if (_v[0]>0.0f && _v[1]<0.0f)
				 return tmp + 90.0f;
		 else
		 if (_v[0]<0.0f && _v[1]<0.0f)
				 return 90.0f - tmp;
		 else
		 if (_v[0]<0.0f && _v[1]>0.0f)
				 return 270.0f - tmp;

		 return tmp;
	}


	// Calculates the angle of this vector in radians in the counter trigonometric sense.
	//0 is to the right (3 o'clock), values increase clockwise.
	inline TPrecision getAngleRAD() const
	{
		 if (_v[1] == 0.0f) // corrected thanks to a suggestion by Jox
				 return _v[0] < 0.0f ? static_cast<TPrecision>(VEC2::PI64) : 0.0f;
		 else if (_v[0] == 0.0f)
				 return _v[1] < 0.0f ? static_cast<TPrecision>(VEC2::PI64)/2.0f : 3.0f*static_cast<TPrecision>(VEC2::PI64)/2.0f;

		 TPrecision tmp = _v[1] / sqrt((_v[0]*_v[0] + _v[1]*_v[1]));
		 tmp = atan( sqrt(1 - tmp*tmp) / tmp);

		 if (_v[0]>0.0f && _v[1]>0.0f)
				 return tmp + 3.0f*static_cast<TPrecision>(VEC2::PI64)/2.0f;
		 else
		 if (_v[0]>0.0f && _v[1]<0.0f)
				 return tmp + static_cast<TPrecision>(VEC2::PI64)/2.0f;
		 else
		 if (_v[0]<0.0f && _v[1]<0.0f)
				 return static_cast<TPrecision>(VEC2::PI64)/2.0f - tmp;
		 else
		 if (_v[0]<0.0f && _v[1]>0.0f)
				 return 3.0f*static_cast<TPrecision>(VEC2::PI64)/2.0f - tmp;

		 return tmp;
	}

	//Calculates the angle between this vector and another one in degree. Returns a value between 0 and 90.
	// Use this instead of getMinAngleTo(.) if possible, it is faster.
	inline TPrecision getAngleWithDEG(const Vec2<TPrecision>& b) const
	{
		/*
		 TPrecision tmp = _v[0]*b._v[0] + _v[1]*b._v[1];

		 if (tmp == 0.0f)
				 return 90.0f;

		 tmp = tmp / sqrt((_v[0]*_v[0] + _v[1]*_v[1]) * (b._v[0]*b._v[0] + b._v[1]*b._v[1]));
		 if (tmp < 0.0f)
				 tmp = -tmp;

		 return atan(sqrt(1.0f - tmp*tmp) / tmp) * static_cast<TPrecision>(RADTODEG64);
		 */
		return fabs(getAngleTrigDEG() - b.getAngleTrigDEG());
	}


	// Calculates the angle between this vector and another one in radians. Returns a value between 0 and PI/2.
	// Use this instead of getMinAngleTo(.) if possible, it is faster.
	inline TPrecision getAngleWithRAD(const Vec2<TPrecision>& b) const
	{
		/*
		 TPrecision tmp = _v[0]*b._v[0] + _v[1]*b._v[1];

		 if (tmp == 0.0f)
				 return static_cast<TPrecision>(PI64)/2.0f;

		 tmp = tmp / sqrt((_v[0]*_v[0] + _v[1]*_v[1]) * (b._v[0]*b._v[0] + b._v[1]*b._v[1]));
		 if (tmp < 0.0f)
				 tmp = -tmp;

		 return atan(sqrt(1.0f - tmp*tmp) / tmp);
		 */
		return fabs(getAngleTrigRAD() - b.getAngleTrigRAD());
	}

protected:

	// Functions
	//====================  32 bits specific methods =========================

	//! returns if a equals b, taking possible rounding errors into account
	inline bool equalsf(const float& a, const float& b, const float tolerance = VEC2::ROUNDING_ERROR_32) const { return (a + tolerance >= b) && (a - tolerance <= b); }

	inline float inverse_sqrtf(const float x)
	{
	#ifdef CF_FAST_MATH
		// comes from Nvidia
	#if 1
		//!<	integer representation of 1.0
		#define IEEE_1_0			0x3f800000

		unsigned int tmp = (unsigned int(IEEE_1_0 << 1) + IEEE_1_0 - *(unsigned int*)&x) >> 1;   
		float y = *(float*)&tmp;                                             
		return y * (1.47f - 0.47f * x * y * y);
	#elif defined(_MSC_VER)
		// an sse2 version
		__asm
		{
			movss	xmm0, x
			rsqrtss	xmm0, xmm0
			movss	x, xmm0
		}
		return x;
	#endif
	#else // no fast math
		return 1.0f / sqrtf( x );
	#endif
	}

	//====================  64 bits specific methods =========================

	//! returns if a equals b, taking possible rounding errors into account
	inline bool equalsd(const double& a, const double& b, const double tolerance = VEC2::ROUNDING_ERROR_64) const { return (a + tolerance >= b) && (a - tolerance <= b); }


	// 32 bits
	//! returns if this vector equals the other one, taking floating point rounding errors into account
	bool equals(const Vec2<float>& other, const TPrecision tolerance = VEC2::ROUNDING_ERROR_32 ) const
	{
		return equalsf(_v[0], other._v[0], tolerance)
				&& equalsf(_v[1], other._v[1], tolerance);
	}

	// 64 btis
	//! returns if this vector equals the other one, taking floating point rounding errors into account
	bool equals(const Vec2<double>& other, const TPrecision tolerance = VEC2::ROUNDING_ERROR_64 ) const
	{
		return equalsd(_v[0], other._v[0], tolerance)
				&& equalsd(_v[1], other._v[1], tolerance);
	}
	
	// Members

	TPrecision _v[2];
};


#endif //VEC2_H