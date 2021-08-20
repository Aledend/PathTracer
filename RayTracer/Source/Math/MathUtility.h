#pragma once
#include <stdlib.h>
#include <math.h>
#include "Math/Vec3.h"

//////////////////// Constants ////////////////////
constexpr float F_PI = 3.1415926f;


//////////////////// Utility Functions ////////////////////
inline float RandNext() {
	return static_cast<float>(rand()) / (RAND_MAX + 1);
}

inline Vec3 RandomInUnitSphere() {

	float z = (RandNext() - 0.5f) * 1.99f;
	float rxy = sqrt(1 - z * z);
	float phi = RandNext() * 2 * F_PI;
	float x = rxy * cos(phi);
	float y = rxy * sin(phi);

	return Vec3(x, y, z).Normalized();
}

inline Vec3 RandomInUnitDisk() {
	Vec3 p = RandomInUnitSphere();
	p.z = 0;
	return p.Normalized();
}

inline Vec3 Reflect(const Vec3& v, const Vec3& n) {
	return v - 2 * Vec3::Dot(v, n) * n;
}

inline bool Refract(const Vec3& v, const Vec3& n, float ni_over_nt, Vec3& refracted) {
	Vec3 uv = v.Normalized();
	float dt = Vec3::Dot(uv, n);
	float discriminant = 1.0f - ni_over_nt * ni_over_nt * (1 - dt * dt);
	if (discriminant > 0)
	{
		refracted = ni_over_nt * (uv - n*dt) - n * sqrt(discriminant);
		return true;
	}
	return false;
}

inline float Schlick(float cosine, float refIdx) {
	float r0 = (1.f - refIdx) / (1.f + refIdx);
	r0 *= r0;
	return r0 + (1.f - r0) * static_cast<float>(pow((1.f - cosine), 5));
}
