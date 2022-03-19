#pragma once
#include <math.h>

class Vec3
{
public:
	Vec3() = default;
	Vec3(float e0, float e1, float e2) { e[0] = e0; e[1] = e1; e[2] = e2; }

	inline const Vec3& operator+() const { return *this; }
	inline Vec3 operator-() const { return Vec3(-e[0], -e[1], -e[2]); }
	inline float operator[](int i) const { return e[i]; }
	inline float& operator[](int i) { return e[i]; }

	inline Vec3& operator+=(const Vec3& v) {
		e[0] += v.x;
		e[1] += v.y;
		e[2] += v.z;
		return *this;
	}
	inline Vec3& operator-=(const Vec3& v) {
		e[0] -= v.x;
		e[1] -= v.y;
		e[2] -= v.z;
		return *this;
	}
	inline Vec3& operator*=(const float t) {
		e[0] *= t;
		e[1] *= t;
		e[2] *= t;
		return *this;
	}
	inline Vec3& operator/=(const float t) {
		e[0] /= t;
		e[1] /= t;
		e[2] /= t;
		return *this;
	}

	inline Vec3 operator*(const float& t) const {
		return Vec3(x * t, y * t, z * t);
	}
	inline Vec3 operator*(const Vec3& v) const {
		return Vec3(x * v.x, y * v.y, z * v.z);
	}
	inline friend Vec3 operator*(float t, Vec3 v) {
		return v * t;
	}
	inline Vec3 operator/(const float& t) const {
		return Vec3(x / t, y / t, z / t);
	}
	inline Vec3 operator+(const Vec3& v) const {
		return Vec3(x + v.x, y + v.y, z + v.z);
	}
	inline Vec3 operator-(const Vec3& v) const {
		return Vec3(x - v.x, y - v.y, z - v.z);
	}

	
	inline static float Dot(const Vec3& v1, const Vec3& v2) {
		return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
	}
	inline static Vec3 Cross(const Vec3& v1, const Vec3& v2) {
		return Vec3(v1.y * v2.z - v1.z * v2.y,
			-(v1.x * v2.z - v1.z * v2.x),
			v1.x * v2.y - v1.y * v2.x);
	}
	inline float Magnitude() const {
		const float sqr_magnitude = x * x + y * y + z * z;
		return static_cast<float>(sqrt(sqr_magnitude));
	}
	inline Vec3 Normalized() const {
		return *this / this->Magnitude();
	}


	union {
		struct{ float x, y, z; };
		struct{ float r, g, b; };
		struct{ float e[3]; };
	};

};