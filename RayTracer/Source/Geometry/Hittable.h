#pragma once
#include "Math/Ray.h"
#include "Math/Vec3.h"

struct HitRecord {
	float t;
	Vec3 p;
	Vec3 normal;
};

class Hittable
{
public:
	virtual bool Hit(const Ray& r, float tMin, float tMax, HitRecord& rec) const = 0;
};
