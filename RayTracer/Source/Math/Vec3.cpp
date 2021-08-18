#include "Vec3.h"

inline std::istream& operator>>(std::istream& is, Vec3& v)
{
	is >> v.e[0] >> v.e[1] >> v.e[2];
	return is;
}

inline std::ostream& operator<<(std::ostream& os, const Vec3& v)
{
	os << v.e[0] << " " << v.e[1] << " " << v.e[2];
	return os;
}

inline Vec3 operator+(const Vec3& v1, const Vec3& v2)
{
	return Vec3(v1.x() + v2.x(), v1.y() + v2.y(), v1.z() + v2.z());
}

inline Vec3 operator-(const Vec3& v1, const Vec3& v2)
{
	return Vec3(v1.x() - v2.x(), v1.y() - v2.y(), v1.z() - v2.z());
}

inline Vec3 operator*(const Vec3& v1, const Vec3& v2)
{
	return Vec3(v1.x() * v2.x(), v1.y() * v2.y(), v1.z() * v2.z());
}

inline Vec3 operator/(const Vec3& v1, const Vec3& v2)
{
	return Vec3(v1.x() / v2.x(), v1.y() / v2.y(), v1.z() / v2.z());
}

inline Vec3 operator*(float t, const Vec3& v)
{
	return Vec3(v.x() * t, v.y() * t, v.z() * t);
}

inline Vec3 operator/(Vec3& v, float t)
{
	return Vec3(v.x() / t, v.y() / t, v.z() / t);
}

inline Vec3 operator*(const Vec3& v, float t)
{
	return Vec3(v.x() * t, v.y() * t, v.z() * t);
}



inline Vec3& Vec3::operator-=(const Vec3& v)
{
	e[0] -= v.x();
	e[1] -= v.y();
	e[2] -= v.z();
	return *this;
}

inline Vec3& Vec3::operator*=(const Vec3& v)
{
	e[0] *= v.x();
	e[1] *= v.y();
	e[2] *= v.z();
	return *this;
}

inline Vec3& Vec3::operator/=(const Vec3& v)
{
	e[0] /= v.x();
	e[1] /= v.y();
	e[2] /= v.z();
	return *this;
}

inline Vec3& Vec3::operator*=(const float t)
{
	e[0] *= t;
	e[1] *= t;
	e[2] *= t;
	return *this;
}




