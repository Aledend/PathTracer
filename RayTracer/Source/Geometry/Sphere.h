#pragma once
#include "Geometry/Hittable.h"
#include "Math/Vec3.h"

class Sphere : public Hittable
{
public:
	Sphere() {}
	Sphere(Vec3 cen, float r) : center(cen), radius(r) {};
	virtual bool Hit(const Ray& r, float tmin, float tMax, HitRecord& rec) const;
	Vec3 center;
	float radius;
};

