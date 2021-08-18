#pragma once
#include <stdlib.h>
#include <math.h>
#include "Math/Vec3.h"

// Constants
constexpr float F_PI = 3.1415926f;



inline float RandNext() {
	return static_cast<float>(rand()) / (RAND_MAX + 1);
}

// Replace with gaussian distribution
inline Vec3 RandomInUnitSphere() {
	Vec3 p;
	do {
		float r1 = RandNext();
		float r2 = RandNext();
		float r3 = RandNext();
		p = 2.0 * Vec3(r1, r2, r3) - Vec3(1, 1, 1);
	} while (p.SqrMagnitude() >= 1.0f);
	return p;
}

inline Vec3 RandomInUnitDisk() {
	Vec3 p;
	do {
		float r1 = RandNext() * 2 - 1;
		float r2 = RandNext() * 2 - 1;
		p = Vec3(r1, r2, 0.f);
	} while (p.SqrMagnitude() >= 1.0f);
	return p;
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
