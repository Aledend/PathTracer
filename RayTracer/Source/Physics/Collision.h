#pragma once
#include <math.h>
#include "Math/Vec3.h"
#include "Math/Ray.h"

class Collision
{
public:
	inline static float HitSphere(const Vec3& center, float radius, const Ray& r) {
		Vec3 oc = r.Origin() - center;
		float a = r.Direction().SqrMagnitude();
		float b = 2.0f * Vec3::Dot(oc, r.Direction());
		float c = Vec3::Dot(oc, oc) - radius * radius;
		float discriminant = b*b - 4*a*c;

		return discriminant < 0 ? -1.0f : (-b - sqrt(discriminant)) / (2.0f*a);
	}
};

