#pragma once
#ifndef ____PBE_OBJECTS_H____
#define ____PBE_OBJECTS_H____

#include "Volume.h"
#include <cmath>

float PI = 3.1415926;

namespace lux
{
	class ConstantField : public Volume<float>
	{

	public:

		ConstantField(float f) : F(f) {}
		~ConstantField() {}//;

		const float eval(const Vector& P) const { float base; return base; }
		const Vector grad(const Vector& P) const { return Vector(0, 0, 0); }


	private:

		float F;
	};

	class ObjectSphere : public Volume<float>
	{

	public:

		ObjectSphere(float r) : R(r) {}
		~ObjectSphere() {};

		float eval(Vector& X) { return R - X.magnitude(); }
		//Vector grad(const Vector& P) const { return Vector(0, 0, 0); }

	private:

		float R;
	};

	class ObjectEllipse : public Volume<float>
	{

	public:

		ObjectEllipse(Vector n, float r1, float r2) : N(n), R1(r1), R2(r2){}
		~ObjectEllipse() {};

		float eval(Vector& X) 
		{
			float Z = X * N;
			float x_per = (X - Z * N).magnitude();
			x_per *= x_per;
			Z *= Z;
			float final = 1.0 - Z / (R1 * R1) - x_per / (R2 * R2);
			return final;
		}

	private:
		Vector N = Vector(1, 0, 0);
		float R1, R2;//R1:major; R2:minor
	};

	class ObjectTorus : public Volume<float>
	{

	public:

		ObjectTorus(Vector n, float r1, float r2) : N(n), R1(r1), R2(r2) {}
		~ObjectTorus() {};

		float eval(Vector& X)
		{
			float x_per = (X - (X * N)*N).magnitude();
			x_per *= x_per;
			float final = 4.0*R1*R1*x_per-(     X.magnitude()* X.magnitude() + (R1 * R1) - (R2 * R2)        ) * (     X.magnitude() * X.magnitude() + (R1 * R1) - (R2 * R2)    );
			return final;
		}

	private:
		Vector N = Vector(1, 0, 0);
		float R1, R2;//R1:major; R2:minor
	};

	class ObjectBox : public Volume<float>
	{

	public:

		ObjectBox(float r, float q) : R(r), Q(q) {}
		~ObjectBox() {};

		float eval(Vector& X)
		{

			float final = pow(R, Q*2) - pow(X.X(), Q*2) - pow(X.Y(), Q * 2) - pow(X.Z(), Q * 2);
			return final;
		}

	private:
		float Q;
		float R;
	};

	class ObjectPlane : public Volume<float>
	{

	public:

		ObjectPlane(Vector x0, Vector n) : X0(x0), N(n) {}
		~ObjectPlane() {};

		float eval(Vector& X)
		{

			float final = -1.0 * (X - X0) * N;
			return final;
		}

	private:
		Vector X0;
		Vector N;
	};

	class ObjectCone : public Volume<float>
	{

	public:

		ObjectCone(float h, float thi_max, Vector n) : H(h), thi(thi_max), N(n) {}
		~ObjectCone() {};

		float eval(Vector& X)
		{
			float temp = X * N;
			if (temp < 0)
			{
				return temp;
			}
			else if (temp > H)
			{
				return H - temp;
			}
			else if (0 < temp < H)
			{
				return X * N - X.magnitude() * cos(thi);
			}
			//return 0;
		}

	private:
		float thi, H;
		Vector N;
	};

	class ObjectCylinder : public Volume<float>
	{

	public:

		ObjectCylinder(float R, Vector N) : r(R), n(N) {}
		~ObjectCylinder() {};

		float eval(Vector& x)
		{

			float final = r - (x - (x * n) * n).magnitude();
			return final;
		}

	private:
		float r;
		Vector n;
	};

	class ObjectIcosahedron : public Volume<float>
	{

	public:

		ObjectIcosahedron(float T) : t(T) {}
		~ObjectIcosahedron() {};

		float eval(Vector& x)
		{
			float temp = x.magnitude();
			if (temp <= 1.8 * PI)
			{
				float final = cos(x.X() + t * x.Y()) + cos(x.X() - t * x.Y())
							+ cos(x.Y() + t * x.Z()) + cos(x.Y() - t * x.Z())
							+ cos(x.Z() - t * x.X()) + cos(x.Z() + t * x.X()) -2;
			}

			else if (temp > 1.8 * PI)
			{
				return -1.8*PI;
			}
		}

	private:
		float t;
	};

	class ObjectSteinerPatch : public Volume<float>
	{

	public:

		ObjectSteinerPatch()  {}
		~ObjectSteinerPatch() {};

		float eval(Vector& x)
		{

			float final = -1.0 * (x.X() * x.X() * x.Y() * x.Y() + x.X() * x.X() * x.Z() * x.Z()
						  +x.Y() * x.Y() * x.Z() * x.Z() - x.X() * x.Y() * x.Z());
			return final;
		}

	private:
		
	};

//==================================================================================================
	class ObjectSFTransform : public Volume<float>
	{

	public:

		ObjectSFTransform(Volume<float> *f, Vector dx) : elem(f), dx(dx) {}
		~ObjectSFTransform() {};

		float eval(Vector& x)
		{
			return elem->eval(x - dx);
		}

	private:
		Volume<float> *elem;
		Vector dx;
	};

	class ObjectSFRotation : public Volume<float>
	{

	public:

		ObjectSFRotation(Volume<float>* f, Vector axis, float angle) : elem(f), ax(axis.unitvector()), a(angle) {}//normalization or unitvector
		~ObjectSFRotation() {};

		float eval(Vector& p)
		{
			Vector vec_rotated = p * cos(a) + ax * (ax * p) * (1.0 - cos(a)) + p ^ ax * sin(a);
			float final = elem->eval(vec_rotated);
			return final;
		}

	private:
		Volume<float>* elem;
		Vector ax;
		float a;
	};

	class ObjectSFScale : public Volume<float>
	{

	public:

		ObjectSFScale(Volume<float>* f, Vector pivot, float scale) : elem(f), p(pivot), s(scale) {}
		~ObjectSFScale() {};

		float eval(Vector& x)
		{
			return elem->eval((x - p)/s+p);
		}

	private:
		Volume<float>* elem;
		Vector p;
		float s;
	};
	//f(  f(   f(x)    )  )

	//=======================================================================================
	//addition
	//cutout
	//intersection
	//blend
	//shell
	//Dilation
	//mask
	//clamp
	//========================================================================================

	// constant color
	class ConstantColor : public Volume<Color>
	{
	public:
		ConstantColor(Color color) : constantC(color) {}
		~ConstantColor() {}

		const Color eval(const Vector& x) const { return constantC; }

	private:
		Color constantC;//(0,0,0,0)
	};
}


#endif