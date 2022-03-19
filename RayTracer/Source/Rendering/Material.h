#pragma once
#include "Math/Vec3.h"

enum class MaterialType : int
{
	Lambertian = 0,
	Metal = 1,
	Dielectric = 2,
};

class Lambertian
{
public:
	Lambertian() = default;
	Lambertian(const Vec3& a) : albedo(a) {}
	Vec3 albedo;
};

class Metal
{
public:
	Metal() = default;
	Metal(const Vec3& a, const float f) : albedo(a)
	{ 
		fuzz = f > 1 ? 1 : f;
	}
	Vec3 albedo;
	float fuzz;
};

class Dielectric {
public:
	Dielectric() = default;
	Dielectric(const float ri) : reflection_index(ri) { }
	float reflection_index; 
};

class Material
{

public:
	Material() = default;
	Material(const Lambertian& mat) : lambertian(mat), matType(MaterialType::Lambertian) {}
	Material(const Metal& mat) : metal(mat), matType(MaterialType::Metal) {}
	Material(const Dielectric& mat) : dielectric(mat), matType(MaterialType::Dielectric) {}

	MaterialType matType;

	union
	{
		Lambertian lambertian;
		Metal metal;
		Dielectric dielectric;
	};
};