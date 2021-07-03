#include "Objects.h"
#include "Camera.h"
#include "Color.h"
#include "Image.h"

int main(int argc, char* argv[])
{
	using namespace std;
	using namespace lux;
	
	Volume<float> densityVolume;
	Volume<Color> colorVolume;
	
	Camera cam;
	Image img;
	Color Lp;
	for (int u = 0;;)
	{
		for (int v = 0;;)
		{
			float T = 1.0;
			float K = 1.0;
			float del_s = 0.05;
			float s = 0.05;
			float alpha;
			float del_t;
			cam.setEyeViewUp(Vector(0, -0, -2.0), Vector(0, 0, 2.0), Vector(0, -0, -2));
			//cam.setRight(); //undefined?
			Vector Xp = cam.eye() + u * cam.right() + v * cam.up();
			Vector Np = (Xp - cam.eye()) / (cam.view() * (Xp - cam.eye())) - cam.view();
			Vector X = cam.eye() + Np * cam.nearPlane();
			while (s < cam.farPlane())
			{
				X = X + del_s * Np;
				del_t = exp(-1.0 * K * del_s * densityVolume.eval(X));
				Lp = Lp + colorVolume.eval(X) * T * (1 - del_t) / K;
				alpha = 1.0 - T;
				T *= del_t;
			}
		}
	}
	

	return 0;
}