#pragma once
#ifndef _UTILS_H
#define _UTILS_H

#include "Vec3.h"

#include <list>
#include <vector>
#include <cmath>

// 2D --> x, z
// works on plane x, z by default

//-----------------------------------------------------------------------------
#define printMe(e) (std::cout << #e" = " << (e) << std::endl << std::flush)

namespace Utils
{

//-----------------------------------------------------------------------------
// Generates a pseudorandom number in [min,max] range.
template<typename real_type>
inline static real_type rand_range(const real_type &min, const real_type &max)
{
	//assert( min <= max );
	return ((static_cast<real_type>( rand() ) / static_cast<real_type>( RAND_MAX )) * (max - min)) + min;
}


//-----------------------------------------------------------------------------
// line parametric equation 
static Vec3<float> parametricLine(const float t, const Vec3<float>& p1, const Vec3<float>& p2)
{
	return p1*(1.0f-t) + p2*t;
}


//-----------------------------------------------------------------------------
// distance in 3D, from point P to line P0_P1
static float distancePointFromLine3D(const Vec3<float>& p, const Vec3<float>& p0, const Vec3<float>& p1)
{
	const Vec3<float> v = p1 - p0;
	const Vec3<float> w = p - p0;

	const float c1 = w*v;
	const float c2 = v*v;
	const float b = c1 / c2;

	const Vec3<float> pb = p0 + v*b;

	return (p-pb).length();

	//or
	//return ((p-p0)^(p-p1)).length()/(p1-p0).length();
}


//-----------------------------------------------------------------------------
// Based on implementation found at:
// http://softsurfer.com/Archive/algorithm_0102/algorithm_0102.htm#dist_Point_to_Line
//-----------------------------------------------------------------------------
static float distancePointFromSegment3D(const Vec3<float>& p, const Vec3<float>& p0, const Vec3<float>& p1)
{
	Vec3<float> v = p1 - p0;
    Vec3<float> w = p - p0;

    const float c1 = w*v;
    if ( c1 <= 0 )
        return (p - p0).length();

    const float c2 = v*v;
    if ( c2 <= c1 )
        return (p - p1).length();

    const float b = c1 / c2;
    Vec3<float> pb = p0 + v*b;
    
	return (p - pb).length();
}

//-----------------------------------------------------------------------------
// Catmull-Rom Spline
static void catmullRom(const int quality, const std::vector< Vec3<float> >& points, std::list< Vec3<float> >& total_points)
{
	//kontrola, ci je dostatocny pocet bodov
	if (points.size() < 4) return;

	float step = 1.0f/quality;	//krok vypoctu
	total_points.clear();		//vycistime pole
	Vec3<float> tmp;			//pomocny bod
	float F[4];					//pomocny vektor
	float t3,t2,t;				//pomocne premenne, t^3 a t^2

	//vypocet poctu segmentov krivky
	int segments = points.size() - 3;

	//vonkajsi cyklus, prebehne tolkokrat, kolko je segmentov
	for (int i=0; i < segments-1; ++i)
	{
		t = 0.0f;	//parameter t
		//kvoli efektivite sa vypocita pocet krokov a pouzije celociselny cyklus
		for (int j = 0; j <= quality; ++j)
		{
			//vypocet jednotlivych parametrov t
			t2 = t*t;
			t3 = t*t*t;
			//vektorovy sucin s maticou
			F[0] = -t3 + 2*t2 - t;
			F[1] = 3*t3 - 5*t2 + 2;
			F[2] = -3*t3 + 4*t2 + t;
			F[3] = t3 - t2;
			
			//sucin s riadiacimi bodmi, zapis do zoznamu bodov
			tmp[0] = (F[0]*points[i].x() +
							F[1]*points[i+1].x() +
							F[2]*points[i+2].x() +
							F[3]*points[i+3].x())/2.0f;

			tmp[2] = (F[0]*points[i].z() +
							F[1]*points[i+1].z() +
							F[2]*points[i+2].z() +
							F[3]*points[i+3].z())/2.0f;

			total_points.push_back(tmp);	//vlozenie do zoznamu bodov
			t+=step;						//inkrementacia kroku
		}
	}
}


//-----------------------------------------------------------------------------
enum IntersectResult { PARALLEL, COINCIDENT, NOT_INTERESECTING, INTERESECTING };
//-----------------------------------------------------------------------------
// Line x Line intersection
template<typename Real_Type>
IntersectResult intersect(	const Vec3<Real_Type>& line1Begin,
							const Vec3<Real_Type>& line1End,
							const Vec3<Real_Type>& line2Begin,
							const Vec3<Real_Type>& line2End, 
							Vec3<Real_Type>& intersection )
{                   
	const float denom = ((line2End.y() - line2Begin.y())*(line1End.x() - line1Begin.x())) 
						- ((line2End.x() - line2Begin.x())*(line1End.y() - line1Begin.y()));

    const float nume_a = ((line2End.x() - line2Begin.x())*(line1Begin.y() - line2Begin.y())) 
						- ((line2End.y() - line2Begin.y())*(line1Begin.x() - line2Begin.x()));

    const float nume_b = ((line1End.x() - line1Begin.x())*(line1Begin.y() - line2Begin.y())) 
						- ((line1End.y() - line1Begin.y())*(line1Begin.x() - line2Begin.x()));

    if (denom == 0.0f)
    {
        if (nume_a == 0.0f && nume_b == 0.0f)
        {
            return COINCIDENT;
        }
        return PARALLEL;
    }

    const float ua = nume_a / denom;
    const float ub = nume_b / denom;

    if(ua >= 0.0f && ua <= 1.0f && ub >= 0.0f && ub <= 1.0f)
    {
        // Get the intersection point.
        intersection[0] = line1Begin.x() + ua*(line1End.x() - line1Begin.x());
        intersection[1] = line1Begin.y() + ua*(line1End.y() - line1Begin.y());

        return INTERESECTING;
    }

    return NOT_INTERESECTING;
}


//-----------------------------------------------------------------------------
static bool inline getIntersection( float fDst1, float fDst2, const Vec3<float>& p1, const Vec3<float>& p2, Vec3<float>& hit) 
{
	if ( (fDst1 * fDst2) >= 0.0f) return false;
	if ( fDst1 == fDst2) return false; 
	
	hit = p1 + (p2-p1) * ( -fDst1/(fDst2-fDst1) );

	return true;
}


//-----------------------------------------------------------------------------
static bool inline inBox( const Vec3<float>& hit, const Vec3<float>& b1, const Vec3<float>& b2, const int axis) 
{
	if ( axis==1 && hit.z() > b1.z() && hit.z() < b2.z() && hit.y() > b1.y() && hit.y() < b2.y()) return true;
	if ( axis==2 && hit.z() > b1.z() && hit.z() < b2.z() && hit.x() > b1.x() && hit.x() < b2.x()) return true;
	if ( axis==3 && hit.x() > b1.x() && hit.x() < b2.x() && hit.y() > b1.y() && hit.y() < b2.y()) return true;

	return false;
}


// Returns true if line (L1, L2) intersects with the box (B1, B2)
// Returns intersection point in Hit
//-----------------------------------------------------------------------------
static bool checkLineBox( const Vec3<float>& b1, const Vec3<float>& b2, const Vec3<float>& l1, const Vec3<float>& l2, Vec3<float>& hit)
{
	if (l2.x() < b1.x() && l1.x() < b1.x()) return false;
	if (l2.x() > b2.x() && l1.x() > b2.x()) return false;
	if (l2.y() < b1.y() && l1.y() < b1.y()) return false;
	if (l2.y() > b2.y() && l1.y() > b2.y()) return false;
	if (l2.z() < b1.z() && l1.z() < b1.z()) return false;
	if (l2.z() > b2.z() && l1.z() > b2.z()) return false;
	
	if (l1.x() > b1.x() && l1.x() < b2.x() &&
		l1.y() > b1.y() && l1.y() < b2.y() &&
		l1.z() > b1.z() && l1.z() < b2.z()) 
	{
			hit = l1; 
			return true;
	}

	if ( (getIntersection( l1.x()-b1.x(), l2.x()-b1.x(), l1, l2, hit) && inBox( hit, b1, b2, 1 ))
	  || (getIntersection( l1.y()-b1.y(), l2.y()-b1.y(), l1, l2, hit) && inBox( hit, b1, b2, 2 )) 
	  || (getIntersection( l1.z()-b1.z(), l2.z()-b1.z(), l1, l2, hit) && inBox( hit, b1, b2, 3 )) 
	  || (getIntersection( l1.x()-b2.x(), l2.x()-b2.x(), l1, l2, hit) && inBox( hit, b1, b2, 1 )) 
	  || (getIntersection( l1.y()-b2.y(), l2.y()-b2.y(), l1, l2, hit) && inBox( hit, b1, b2, 2 )) 
	  || (getIntersection( l1.z()-b2.z(), l2.z()-b2.z(), l1, l2, hit) && inBox( hit, b1, b2, 3 )))
		return true;

	return false;
}


/**
*/
static bool isFloatValid(float n)
{
	switch (_fpclass(n))
	{
	// Negative normalized non-zero.
	case _FPCLASS_NN:
	// Negative denormalized.
	case _FPCLASS_ND:
	// Negative zero (-0).
	case _FPCLASS_NZ:
	// Positive 0 (+0).
	case _FPCLASS_PZ:
	// Positive denormalized.
	case _FPCLASS_PD:
	// Positive normalized non-zero.
	case _FPCLASS_PN:
		return true;
	// Signaling NaN.
	case _FPCLASS_SNAN:
	// Quiet NaN.
	case _FPCLASS_QNAN:
	// Negative infinity (-INF).
	case _FPCLASS_NINF:
	// Positive infinity (+INF).
	case _FPCLASS_PINF:
		return false;
	}

	return false;
}

};

#endif // _UTILS_H