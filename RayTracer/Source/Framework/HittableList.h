#pragma once
#include "Geometry/Hittable.h"

class HittableList : public Hittable
{
public:
	HittableList() {}
	HittableList(Hittable** l, int n) {list = l; listSize = n;}

	virtual bool Hit(const Ray& r, float tMin, float tMax, HitRecord& rec) const;

	Hittable** list;
	int listSize;
};

