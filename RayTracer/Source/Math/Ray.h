#pragma once
#ifndef RAY_H
#define RAY_H

#include "Math/Vec3.h"

class Ray
{
public:
	Ray() = default;
	Ray(const Vec3& a, const Vec3& b) { A = a; B = b; }
	Vec3 Origin() const { return A; }
	Vec3 Direction() const { return B; }
	Vec3 PointAtParameter(float t) const { return A + B * t; }

	Vec3 A;
	Vec3 B;
};

#endif // !RAY_H