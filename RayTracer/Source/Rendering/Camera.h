#pragma once
#include "Math/Vec3.h"
#include "Math/Ray.h"
#include "Math/MathUtility.h"
#include <cmath>

class Camera
{
public:
	Camera(Vec3 lookFrom, Vec3 lookAt, Vec3 vUp, float vfov, float aspect, float aperture, float focusDist) {
		lensRadius = aperture * 0.5f;
		float theta = vfov * F_PI / 180;
		float halfHeight = tan(theta*0.5f) * 2.f;
		float halfWidth = aspect * halfHeight;
		origin = lookFrom;
		w = (lookFrom - lookAt).Normalized();
		u = Vec3::Cross(vUp,w).Normalized();
		v = Vec3::Cross(w, u);
		horizontal = halfWidth * focusDist * u;
		vertical = halfHeight * focusDist * v;
		lowerLeftCorner = origin - horizontal * 0.5f - vertical * 0.5f - focusDist * w;
	}

	Ray GetRay(float s, float t) {
		Vec3 rd = lensRadius * RandomInUnitDisk();
		Vec3 offset = u * rd.x() + v * rd.y();
		return Ray(origin + offset, lowerLeftCorner + s*horizontal + t*vertical - origin - offset); 
	}

	Vec3 origin;
	Vec3 lowerLeftCorner;
	Vec3 horizontal;
	Vec3 vertical;
	Vec3 u, v, w;
	float lensRadius;
};

