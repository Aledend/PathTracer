#pragma once
#include "Geometry/Hittable.h"

class HittableList : public Hittable
{
public:
	HittableList() : list(0), listSize(0) {}
	HittableList(Hittable** l, int n) {list = l; listSize = n;}

	~HittableList() { for(int i = 0; i < listSize; i++) delete list[i]; }

	virtual bool Hit(const Ray& r, float tMin, float tMax, HitRecord& rec) const;

	Hittable** list;
	int listSize;
};

