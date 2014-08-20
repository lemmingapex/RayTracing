#ifndef POINT_H
#define	POINT_H

#include <math.h>

class Point
{
	private:
		double x,y,z;
	public:	
		double X() { return x; }
		double Y() { return y; }
		double Z() { return z; }
	
		Point()
		{ }

		Point(double A[3])
		{
			x=A[0];
			y=A[1];
			z=A[2];
		}

		Point(const double& A, const double& B, const double& C)
		{
			x=A;
			y=B;
			z=C;
		}

		Point operator+ (const Point& a) const
		{
			Point result;
			result.x = x+a.x;
			result.y = y+a.y;
			result.z = z+a.z;
			return result;
		}

		Point operator+ (const double& a) const
		{
			Point result;
			result.x = x+a;
			result.y = y+a;
			result.z = z+a;
			return result;
		}

		Point operator- (const Point& a) const
		{
			Point result;
			result.x = x-a.x;
			result.y = y-a.y;
			result.z = z-a.z;
			return result;
		}

		Point operator- (const double& a) const
		{
			Point result;
			result.x = x-a;
			result.y = y-a;
			result.z = z-a;
			return result;
		}

		Point operator* (const Point& a) const
		{
			Point result;
			result.x = x*a.x;
			result.y = y*a.y;
			result.z = z*a.z;
			return result;
		}

		Point operator* (const double& a) const
		{
			Point result;
			result.x = x*a;
			result.y = y*a;
			result.z = z*a;
			return result;
		}

		Point operator/ (const Point& a) const
		{
			Point result;
			result.x = x/a.x;
			result.y = y/a.y;
			result.z = z/a.z;
			return result;
		}

		Point operator/ (const double& a) const
		{
			Point result;
			result.x = x/a;
			result.y = y/a;
			result.z = z/a;
			return result;
		}

		Point normalize()
		{
			Point result = Point(x,y,z);
			result = result/(result.magnitude());
			return result;
		}

		double magnitude()
		{
			return sqrt((x*x)+(y*y)+(z*z));
		}

		double dot(const Point& P)
		{
			return x*P.x+y*P.y+z*P.z;
		}

		Point cross(const Point& P)
		{
			return Point((y*P.z-z*P.y),(z*P.x-x*P.z),(x*P.y-y*P.x));
		}
};

#endif	/* POINT_H */