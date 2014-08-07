#include <math.h>

class point
{
	private:
		double x,y,z;
	public:	
		double X() { return x; }
		double Y() { return y; }
		double Z() { return z; }
	
		point()
		{ }

		point(double A[3])
		{
			x=A[0];
			y=A[1];
			z=A[2];
		}

		point(const double& A, const double& B, const double& C)
		{
			x=A;
			y=B;
			z=C;
		}

		point operator+ (const point& a) const
		{
			point result;
			result.x = x+a.x;
			result.y = y+a.y;
			result.z = z+a.z;
			return result;
		}

		point operator+ (const double& a) const
		{
			point result;
			result.x = x+a;
			result.y = y+a;
			result.z = z+a;
			return result;
		}

		point operator- (const point& a) const
		{
			point result;
			result.x = x-a.x;
			result.y = y-a.y;
			result.z = z-a.z;
			return result;
		}

		point operator- (const double& a) const
		{
			point result;
			result.x = x-a;
			result.y = y-a;
			result.z = z-a;
			return result;
		}

		point operator* (const point& a) const
		{
			point result;
			result.x = x*a.x;
			result.y = y*a.y;
			result.z = z*a.z;
			return result;
		}

		point operator* (const double& a) const
		{
			point result;
			result.x = x*a;
			result.y = y*a;
			result.z = z*a;
			return result;
		}

		point operator/ (const point& a) const
		{
			point result;
			result.x = x/a.x;
			result.y = y/a.y;
			result.z = z/a.z;
			return result;
		}

		point operator/ (const double& a) const
		{
			point result;
			result.x = x/a;
			result.y = y/a;
			result.z = z/a;
			return result;
		}

		point normalize()
		{
			point result = point(x,y,z);
			result = result/(result.magnitude());
			return result;
		}

		double magnitude()
		{
			return sqrt((x*x)+(y*y)+(z*z));
		}

		double dot(const point& P)
		{
			return x*P.x+y*P.y+z*P.z;
		}

		point cross(const point& P)
		{
			return point((y*P.z-z*P.y),(z*P.x-x*P.z),(x*P.y-y*P.x));
		}
};