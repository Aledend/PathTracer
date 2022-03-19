#pragma once
#include "Math/Vec3.h"
#define _USE_MATH_DEFINES
#include <math.h>

class Camera
{
public:
	Camera() = default;

	void Set(const Vec3& look_from, const Vec3& look_at, const Vec3& v_up, float v_fov, float aspect, float aperture, float focus_dist) 
	{
		lens_radius = aperture * 0.5f;
		const float half_theta = v_fov * static_cast<float>(M_PI / 180) * 0.5f;
		const float half_height = static_cast<float>(tan(half_theta)) * 2.f;
		const float half_width = aspect * half_height;
		origin = look_from;
		w = (look_from - look_at).Normalized();
		u = Vec3::Cross(v_up, w).Normalized();
		v = Vec3::Cross(w, u);
		horizontal = half_width * focus_dist * u;
		vertical = half_height * focus_dist * v;
		lowerLeftCorner = origin - horizontal * 0.5f - vertical * 0.5f - focus_dist * w;
	}

	Vec3 origin;
	Vec3 lowerLeftCorner;
	Vec3 horizontal;
	Vec3 vertical;
	Vec3 u, v, w;
	float lens_radius;
};

