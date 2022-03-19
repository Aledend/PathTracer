#pragma once
#include "Rendering/Material.h"
#include "Math/Vec3.h"
#include <algorithm>

class Sphere
{
public:
	Sphere() = default;
	Sphere(const Vec3& cen, float r, const Material& mat, const Vec3& velocity) : center(cen), radius(r), material(mat), velocity(velocity) {};

	Vec3 center;
	float radius;
	Material material;

	void MoveWithinBoundingBox(const Vec3& bb_pos, const Vec3& bb_half_size, const float delta_time)
	{
		const Vec3 max_pos = bb_pos + bb_half_size - Vec3(radius, radius, radius);
		const Vec3 min_pos = bb_pos - bb_half_size + Vec3(radius, radius, radius);
		 
		Vec3 new_pos = center + velocity * delta_time;

		if (new_pos.x > max_pos.x)
		{
			new_pos.x = max_pos.x - (new_pos.x - max_pos.x);
			velocity.x *= -1;
		}
		else if (new_pos.x < min_pos.x)
		{
			new_pos.x = min_pos.x + (min_pos.x - new_pos.x);
			velocity.x *= -1;
		}

		if (new_pos.y > max_pos.y)
		{
			new_pos.y = max_pos.y - (new_pos.y - max_pos.y);
			velocity.y *= -1;
		}
		else if (new_pos.y < min_pos.y)
		{
			new_pos.y = min_pos.y + (min_pos.y - new_pos.y);
			velocity.y *= -1;
		}

		if (new_pos.z > max_pos.z)
		{
			new_pos.z = max_pos.z - (new_pos.z - max_pos.z);
			velocity.z *= -1;
		}
		else if (new_pos.z < min_pos.z)
		{
			new_pos.z = min_pos.z + (min_pos.z - new_pos.z);
			velocity.z *= -1;
		}

		center = new_pos;
	}

private:
	Vec3 velocity;
};

