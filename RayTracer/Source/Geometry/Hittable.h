#pragma once
#include "Math/Ray.h"
#include "Math/Vec3.h"

class Material;

struct HitRecord {

	HitRecord() = default;
	
	float t;
	Vec3 p;
	Vec3 normal;
	Material* mat_ptr;
};

class Hittable
{
public:
	virtual bool Hit(const Ray& r, float tMin, float tMax, HitRecord& rec) const = 0;
};