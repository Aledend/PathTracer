#pragma once
#include "Math/Ray.h"
#include "Math/Vec3.h"

class Material;

struct HitRecord {

	HitRecord() : t(0.f), mat_ptr(0) {};
	
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