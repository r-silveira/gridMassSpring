//
//      Vec3<TPrecision>.h        Basic vector class (TPrecision precision)
//

#ifndef _VEC3_H
#define _VEC3_H

#include <cmath>

namespace VEC3 
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
		const float  ROUNDING_ERROR_32	= 0.000001f;
		const double ROUNDING_ERROR_64	= 0.00000001f;

#endif

	//! Constant for PI.
	//const float PI = 3.14159265359f;
	//! 32bit Constant for converting from degrees to radians
	//const float DEGTORAD = PI / 180.0f;
	//! 32bit constant for converting from radians to degrees (formally known as GRAD_PI)
	//const float RADTODEG   = 180.0f / PI;

	//! Constant for 64bit PI.
	const double PI64 = 3.1415926535897932384626433832795028841971693993751;
	//! 64bit constant for converting from degrees to radians (formally known as GRAD_PI2)
	const double DEGTORAD64 = PI64 / 180.0;
	//! 64bit constant for converting from radians to degrees
	const double RADTODEG64 = 180.0 / PI64;
};


//-----------------------------------------------------------------------------
template <class TPrecision>
class Vec3
{
public:

	// Constructor that sets all components of the vector to zero
    Vec3() { _v[0]=0.0; _v[1]=0.0; _v[2]=0.0; }
	Vec3(const TPrecision& x, const TPrecision& y, const TPrecision& z) { _v[0] = x; _v[1] = y; _v[2] = z; }

	inline Vec3(const Vec3& vec) { _v[0] = vec._v[0]; _v[1] = vec._v[1]; _v[2] = vec._v[2]; }

	~Vec3() {};

	// To be used with OpenGL functions
	inline TPrecision* ptr() { return _v; }
    inline const TPrecision* ptr() const { return _v; }

    inline void set(const TPrecision& x, const TPrecision& y, const TPrecision& z) { _v[0]=x; _v[1]=y; _v[2]=z; }
    inline void set( const Vec3<TPrecision>& rhs) { _v[0] = rhs._v[0]; _v[1] = rhs._v[1]; _v[2] = rhs._v[2]; }

    inline TPrecision& operator [] (int i) { return _v[i]; }
    inline TPrecision  operator [] (int i) const { return _v[i]; }

    inline TPrecision& x() { return _v[0]; }
    inline TPrecision& y() { return _v[1]; }
    inline TPrecision& z() { return _v[2]; }

    inline TPrecision x() const { return _v[0]; }
    inline TPrecision y() const { return _v[1]; }
    inline TPrecision z() const { return _v[2]; }

	// Dot product.
    inline TPrecision operator * (const Vec3<TPrecision>& rhs) const
    {
        return _v[0]*rhs._v[0] + _v[1]*rhs._v[1] + _v[2]*rhs._v[2];
    }


    // Cross product.
    inline const Vec3<TPrecision> operator ^ (const Vec3<TPrecision>& rhs) const
    {
        return Vec3<TPrecision>(_v[1]*rhs._v[2]-_v[2]*rhs._v[1],
								_v[2]*rhs._v[0]-_v[0]*rhs._v[2] ,
								_v[0]*rhs._v[1]-_v[1]*rhs._v[0]);
    }

    // Multiply by scalar. 
    inline const Vec3<TPrecision> operator * (const TPrecision rhs) const
    {
        return Vec3<TPrecision>(_v[0]*rhs, _v[1]*rhs, _v[2]*rhs);
    }

    // Unary multiply by scalar. 
    inline Vec3<TPrecision>& operator *= (const TPrecision rhs)
    {
        _v[0]*=rhs;
        _v[1]*=rhs;
        _v[2]*=rhs;
        return *this;
    }

	// Divide by scalar. 
    inline const Vec3<TPrecision> operator / (const TPrecision rhs) const
    {
        return Vec3<TPrecision>(_v[0]/rhs, _v[1]/rhs, _v[2]/rhs);
    }

    // Unary divide by scalar.
    inline Vec3<TPrecision>& operator /= (const TPrecision rhs)
    {
        _v[0]/=rhs;
        _v[1]/=rhs;
        _v[2]/=rhs;
        return *this;
    }

	 // Binary vector add. 
    inline const Vec3<TPrecision> operator + (const Vec3<TPrecision>& rhs) const
    {
        return Vec3<TPrecision>(_v[0]+rhs._v[0], _v[1]+rhs._v[1], _v[2]+rhs._v[2]);
    }

    // Unary vector add. Slightly more efficient because no temporary
    // intermediate object.
    inline Vec3<TPrecision>& operator += (const Vec3<TPrecision>& rhs)
    {
        _v[0] += rhs._v[0];
        _v[1] += rhs._v[1];
        _v[2] += rhs._v[2];
        return *this;
    }

    // Binary vector subtract. 
    inline const Vec3<TPrecision> operator - (const Vec3<TPrecision>& rhs) const
    {
        return Vec3<TPrecision>(_v[0]-rhs._v[0], _v[1]-rhs._v[1], _v[2]-rhs._v[2]);
    }

    // Unary vector subtract.
    inline Vec3<TPrecision>& operator -= (const Vec3<TPrecision>& rhs)
    {
        _v[0]-=rhs._v[0];
        _v[1]-=rhs._v[1];
        _v[2]-=rhs._v[2];
        return *this;
    }


    // Negation operator. Returns the negative of the Vec3d.
    inline const Vec3<TPrecision> operator - () const
    {
        return Vec3<TPrecision>(-_v[0], -_v[1], -_v[2]);
    }

	// I tested in Windows and it works as expected
	//inline bool operator == (const Vec3<TPrecision>& v) const { return (_v[0]==v._v[0]) && (_v[1]==v._v[1]) && (_v[2]==v._v[2]); }
	//inline bool operator != (const Vec3<TPrecision>& v) const { return _v[0]!=v._v[0] || _v[1]!=v._v[1] || _v[2]!=v._v[2]; }

	// The methods below are slightly faster than the built-in operator
	// 32 bits
	bool operator==(const Vec3<float>& other) const
	{
		return this->equals(other);
	}

	bool operator!=(const Vec3<float>& other) const
	{
		return !this->equals(other);
	}

	// 64 bits
	bool operator==(const Vec3<double>& other) const
	{
		return this->equals(other);
	}

	bool operator!=(const Vec3<double>& other) const
	{
		return !this->equals(other);
	}

	// It's just a definition
	inline bool operator <  (const Vec3<TPrecision>& v) const
    {
        if (_v[0]<v._v[0]) return true;
        else if (_v[0]>v._v[0]) return false;
        else if (_v[1]<v._v[1]) return true;
        else if (_v[1]>v._v[1]) return false;
        else return (_v[2]<v._v[2]);
    }

	// Length of the vector = sqrt( vec . vec ) 
	// using float to works with TPrecision == intenger
	// Anyway, float has enough precision to length()
    inline TPrecision length() const
    {
        return sqrt(_v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2] );
    }

    // Length squared of the vector = vec . vec
    inline TPrecision length2() const
    {
        return _v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2];
    }

	//! Normalizes the vector. In case of the 0 vector the result
		//! is still 0, otherwise the length of the vector will be 1.
		//! Todo: 64 Bit template doesnt work.. need specialized template
	inline void normalize()
	{
		TPrecision norm2 = _v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2];
        if (norm2 > 0.0)
        {
			norm2 = static_cast<TPrecision>(inverse_sqrtf( static_cast<float>(norm2) ));
            _v[0] *= norm2;
            _v[1] *= norm2;
            _v[2] *= norm2;
        }                
	}

	// multiply by vector components.
	inline void componentMultiply(const Vec3<TPrecision>& rhs)
	{
		_v[0]*=rhs[0];
		_v[1]*=rhs[1];
		_v[2]*=rhs[2];
	}

	// divide rhs components by rhs vector components.
	inline void componentDivide(const Vec3<TPrecision>& rhs)
	{
		_v[0]/=rhs[0];
		_v[1]/=rhs[1];
		_v[2]/=rhs[2];
	}	

	//! Returns distance from another point.
	/// Here, the vector is interpreted as point in 3 dimensional space.
	TPrecision distanceFrom(const Vec3<TPrecision>& other) const
	{
		return Vec3<TPrecision>(_v[0] - other._v[0], _v[1] - other._v[1], _v[2] - other._v[2]).length();
	}

	//! Returns squared distance from another point.
	// Here, the vector is interpreted as point in 3 dimensional space.
	TPrecision distance2From(const Vec3<TPrecision>& other) const
	{
		return Vec3<TPrecision>(_v[0] - other._v[0], _v[1] - other._v[1], _v[2] - other._v[2]).length2();
	}

	//! Returns if this vector interpreted as a point is on a line between two other points.
	// WARNING: It is assumed that the point is on the line.
	//! \param begin: Beginning vector to compare between.
	//! \param end: Ending vector to compare between.
	//! \return True if this vector is between begin and end.  False if not.
	bool isBetweenPoints(const Vec3<TPrecision>& begin, const Vec3<TPrecision>& end) const
	{
		const TPrecision f = (end - begin).length2();
		return distance2From(begin) < f && distance2From(end) < f;
	}

	//! Rotates the vector by a specified number of degrees around the Y
	//! axis and the specified center.
	//! \param degrees: Number of degrees to rotate around the Y axis.
	//! \param center: The center of the rotation.
	void rotateY(TPrecision degrees, const Vec3<TPrecision>& center)
	{
		degrees *= static_cast<TPrecision>(VEC3::DEGTORAD64);
		const TPrecision cs = cos(degrees);
		const TPrecision sn = sin(degrees);
		_v[0] -= center._v[0];
		_v[2] -= center._v[2];
		set(_v[0]*cs - _v[2]*sn, _v[1], _v[0]*sn + _v[2]*cs);
		_v[0] += center._v[0];
		_v[2] += center._v[2];
	}

	//! Rotates the vector by a specified number of degrees around the Z
	//! axis and the specified center.
	//! \param degrees: Number of degrees to rotate around the Z axis.
	//! \param center: The center of the rotation.
	void rotateZ(TPrecision degrees, const Vec3<TPrecision>& center)
	{
		degrees *= static_cast<TPrecision>(VEC3::DEGTORAD64);
		const TPrecision cs = cos(degrees);
		const TPrecision sn = sin(degrees);
		_v[0] -= center._v[0];
		_v[1] -= center._v[1];
		set(_v[0]*cs - _v[1]*sn, _v[0]*sn + _v[1]*cs, _v[2]);
		_v[0] += center._v[0];
		_v[1] += center._v[1];
	}

	//! Rotates the vector by a specified number of degrees around the X
	//! axis and the specified center.
	//! \param degrees: Number of degrees to rotate around the X axis.
	//! \param center: The center of the rotation.
	void rotateX(TPrecision degrees, const Vec3<TPrecision>& center)
	{
		degrees *= static_cast<TPrecision>(VEC3::DEGTORAD64);
		const TPrecision cs = cos(degrees);
		const TPrecision sn = sin(degrees);
		_v[2] -= center._v[2];
		_v[1] -= center._v[1];
		set(_v[0], _v[1]*cs - _v[2]*sn, _v[1]*sn + _v[2]*cs);
		_v[2] += center._v[2];
		_v[1] += center._v[1];
	}

	//! Returns interpolated vector.
	// \param other: other vector to interpolate between
	//\param d: value between 0.0f and 1.0f. 
	Vec3<TPrecision> getInterpolated(const Vec3<TPrecision>& other, const TPrecision d) const
	{
		const TPrecision inv = static_cast<TPrecision>(1.0 - d);
		return Vec3<TPrecision>(other._v[0]*inv + _v[0]*d, other._v[1]*inv + _v[1]*d, other._v[2]*inv + _v[2]*d);
	}
	
	//! Returns interpolated vector. ( quadratic )
	// \param v2: second vector to interpolate with
	//\param v3: third vector to interpolate with
	//\param d: value between 0.0f and 1.0f. 
	Vec3<TPrecision> getInterpolatedQuadratic(const Vec3<TPrecision>& v2, const Vec3<TPrecision>& v3, const TPrecision d) const
	{
		// this*(1-d)*(1-d) + 2 * v2 * (1-d) + v3 * d * d;
		const TPrecision inv = static_cast<TPrecision>(1.0 - d);
		const TPrecision mul0 = inv * inv;
		const TPrecision mul1 = static_cast<TPrecision>(2.0 * d * inv);
		const TPrecision mul2 = d * d;

		return Vec3<TPrecision> ( _v[0] * mul0 + v2._v[0] * mul1 + v3._v[0] * mul2,
								_v[1] * mul0 + v2._v[1] * mul1 + v3._v[1] * mul2,
								_v[2] * mul0 + v2._v[2] * mul1 + v3._v[2] * mul2);
	}


	inline void zeros() { _v[0]=0.0; _v[1]=0.0; _v[2]=0.0; };
	/*
	//! Gets the Y and Z rotations of a vector.
	// Thanks to Arras on the Irrlicht forums to add this method.
	//return A vector representing the rotation in degrees of
	//this vector. The Z component of the vector will always be 0.
	//              
	//             y 
	//            |
	//            |___ x
	//         z /  
	//
	Vec3<TPrecision> getHorizontalAngle()
	{
		Vec3<TPrecision> angle;

		angle._v[1] = atan2(_v[0], _v[2]);
		angle._v[1] *= static_cast<TPrecision>(RADTODEG64);

		if (angle._v[1] < 0.0f) angle._v[1] += 360.0f;
		if (angle._v[1] >= 360.0f) angle._v[1] -= 360.0f;

		const TPrecision z1 = sqrt(_v[0]*_v[0] + _v[2]*_v[2]);

		angle._v[0] = atan2(z1, _v[1]);
		angle._v[0] *= static_cast<TPrecision>(RADTODEG64);
		angle._v[0] -= 90.0f;

		if (angle._v[0] < 0.0f) angle._v[0] += 360.0f;
		if (angle._v[0] >= 360.0f) angle._v[0] -= 360.0f;

		return angle;
	}*/

protected:

	// Functions

	//====================  32 bits specific methods =========================

	//! returns if a equals b, taking possible rounding errors into account
	inline bool equalsf(const float& a, const float& b, const float tolerance = VEC3::ROUNDING_ERROR_32) const { return (a + tolerance >= b) && (a - tolerance <= b); }

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
	inline bool equalsd(const double& a, const double& b, const double tolerance = VEC3::ROUNDING_ERROR_64) const { return (a + tolerance >= b) && (a - tolerance <= b); }

	// 32 bits
	//! returns if this vector equals the other one, taking floating point rounding errors into account
	bool equals(const Vec3<float>& other, const TPrecision tolerance = VEC3::ROUNDING_ERROR_32 ) const
	{
		return equalsf(_v[0], other._v[0], tolerance)
				&& equalsf(_v[1], other._v[1], tolerance)
				&& equalsf(_v[2], other._v[2], tolerance);
	}

	// 64 btis
	//! returns if this vector equals the other one, taking floating point rounding errors into account
	bool equals(const Vec3<double>& other, const TPrecision tolerance = VEC3::ROUNDING_ERROR_64 ) const
	{
		return equalsd(_v[0], other._v[0], tolerance)
				&& equalsd(_v[1], other._v[1], tolerance)
				&& equalsd(_v[2], other._v[2], tolerance);
	}
	

	// Members

	TPrecision _v[3];
};


#endif //VEC3_H

