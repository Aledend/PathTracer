#include <math.h>
#include "Sphere.h"

bool Sphere::Hit(const Ray& r, float tMin, float tMax, HitRecord& rec) const {
	Vec3 oc = r.Origin() - center;
	float a = r.Direction().SqrMagnitude();
	float b = Vec3::Dot(oc, r.Direction());
	float c = oc.SqrMagnitude() - radius * radius;
	float discriminant = b*b - a*c;
	if (discriminant > 0)
	{
		float temp = (-b - sqrt(b*b-a*c)) / a;
		if (temp < tMax && temp > tMin)
		{
			rec.t = temp;
			rec.p = r.PointAtParameter(rec.t);
			rec.normal = (rec.p - center) / radius;
			return true;
		}
		temp = (-b + sqrt(b*b-a*c)) / a;
		if (temp < tMax && temp > tMin) {
			rec.t = temp;
			rec.p = r.PointAtParameter(rec.t);
			rec.normal = (rec.p - center) / radius;
			return true;
		}
	}
	return false;
}