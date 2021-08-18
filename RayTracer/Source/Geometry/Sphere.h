#pragma once
#include "Rendering/Material.h"
#include "Geometry/Hittable.h"
#include "Math/Vec3.h"

class Sphere : public Hittable
{
public:
	Sphere() : radius(0), material(0) {}
	Sphere(Vec3 cen, float r, Material* mat) : center(cen), radius(r), material(mat) {};
	~Sphere() { delete material; material = nullptr; }

	virtual bool Hit(const Ray& r, float tmin, float tMax, HitRecord& rec) const;
	Material* material; 
	Vec3 center;
	float radius;
};

