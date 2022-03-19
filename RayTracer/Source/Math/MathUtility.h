#pragma once
#include <stdlib.h>
#include "Vec3.h"
#define _USE_MATH_DEFINES
#include <math.h>

//////////////////// Utility Functions ////////////////////
inline float RandNext() {
	return static_cast<float>(rand()) / (RAND_MAX + 1);
}

inline Vec3 RandomDirection() {
	float z = (RandNext() - 0.5f) * 1.99f;
	float inv_sqr = 1 - z * z;
	float rxy = static_cast<float>(sqrt(inv_sqr));
	float phi = RandNext() * 2 * static_cast<float>(M_PI);
	float x = rxy * static_cast<float>(cos(phi));
	float y = rxy * static_cast<float>(sin(phi));

	return Vec3(x, y, z).Normalized();
}