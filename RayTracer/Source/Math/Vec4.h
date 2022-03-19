#pragma once
#include "Math/Vec3.h"

class Vec4
{
public:
	Vec4() = default;
	Vec4(float e0, float e1, float e2, float e3) { e[0] = e0; e[1] = e1; e[2] = e2; e[3] = e3; }
	Vec4(Vec3& v, float e3) { e[0] = v.x; e[1] = v.y; e[2] = v.z; e[3] = e3; }
	Vec4(Vec3& v) { e[0] = v.x; e[1] = v.y; e[2] = v.z; e[3] = 0.0f; }
	Vec4(const Vec3& v) { e[0] = v.x; e[1] = v.y; e[2] = v.z; e[3] = 0.0f; }

	inline Vec4& operator=(const Vec3& v) {
		x = v.x;
		y = v.y;
		z = v.z;
		return *this;
	}

	union {
		struct{ float e[4]; };
		struct{ float x, y, z, w; };
		struct{ float r, g, b, a; };
	};
};
