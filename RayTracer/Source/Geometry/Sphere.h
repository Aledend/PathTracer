#pragma once
#include "Rendering/Material.h"
#include "Math/Vec3.h"

class Sphere
{
public:
	Sphere() = default;
	Sphere(const Vec3& cen, float r, const Material& mat) : center(cen), radius(r), material(mat) {};

	Vec3 center;
	float radius;
	Material material;
};

