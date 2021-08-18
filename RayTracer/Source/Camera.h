#pragma once
#include "Math/Vec3.h"
#include "Math/Ray.h"

class Camera
{
public:
	Camera() {
		lowerLeftCorner = Vec3(-2.0f, -1.5f, -1.0f);
		horizontal = Vec3(4.0f, 0.0f, 0.0f);
		vertical = Vec3(0.0f, 3.0f, 0.0f);
		origin = Vec3(0.0f, 0.0f, 0.0f);
	}
	Ray GetRay(float u, float v) { return Ray(origin, lowerLeftCorner + u*horizontal + v*vertical - origin); }

	Vec3 origin;
	Vec3 lowerLeftCorner;
	Vec3 horizontal;
	Vec3 vertical;
	
};

