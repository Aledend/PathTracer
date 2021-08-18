#include "HittableList.h"

bool HittableList::Hit(const Ray& r, float tMin, float tMax, HitRecord& rec) const {
	HitRecord tempRec;
	bool hitAnything = false;
	double closestSoFar = tMax;
	for (int i = 0; i < listSize; i++)
	{
		if (list[i]->Hit(r, tMin, static_cast<float>(closestSoFar), tempRec)) {
			hitAnything = true;
			closestSoFar = tempRec.t;
			rec = tempRec;
		}
	}
	return hitAnything;
}