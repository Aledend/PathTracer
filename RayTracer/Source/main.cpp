#include <stdio.h>
#include <stdlib.h>
#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <GL/GL.h>
#include <limits>
#include "Math/MathUtility.h"
#include "Math/Vec3.h"
#include "Math/Vec4.h"
#include "Math/Ray.h"
#include "Physics/Collision.h"
#include "Geometry/Hittable.h"
#include "Geometry/Sphere.h"
#include "Framework/HittableList.h"
#include "Framework/Perlin.h"
#include "Rendering/Camera.h"




const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 600;
const int PIXEL_COUNT = WINDOW_HEIGHT * WINDOW_WIDTH;
const int WORK_GROUP_SIZE = 300;
const int antiAliasingIterations = 60;


GLuint camSSbo;
Camera cam;

struct computestruct {
	float u, v, parameter, t;
};
struct colorstruct {
	float r, g, b, a;
};

struct spherestruct {
	Vec4 center;
	Vec4 albedo;
	float fuzz;
	float reflectionIndex;
	int matType;
	float padd;
};

struct camerastruct {
	Vec4 origin;
	Vec4 lowerLeftCorner;
	Vec4 horizontal;
	Vec4 vertical;
	Vec4 u, v, w;
	float lensRadius;
};

#pragma region Compute Shader SRC
const char* CALC_SRC = R"(#version 430 compatibility
//#extension GL_ARB_compute_shader: enable
//#extension GL_ARB_shader_storage_buffer_object: enable

uniform sampler2D permTexture;

struct Camera {
	vec4 origin;
	vec4 lowerLeftCorner;
	vec4 horizontal;
	vec4 vertical;
	vec4 u;
	vec4 v;
	vec4 w;
	float lensRadius;
};


struct SphereData {
	vec4 center; // xyz = center, w = radius
	vec4 albedo;
	float fuzz;
	float reflectionIndex;
	int matType;
	float padd;
};

struct Ray {
	vec4 A;
	vec4 B;
};

struct HitRecord {
	vec4 p;
	vec4 normal;
	float t;
	float pads[3];
};

layout(std430, binding = 3) buffer cam 
{
	Camera camera;
};

layout(std430, binding = 4) buffer Dat
{
	vec4 point_input[];
};

layout (std430, binding = 5) buffer Col
{
	vec4 color_output[];
};

layout (std430, binding = 6) buffer Spheres
{
	SphereData sphere_input[];
};

layout (std430, binding = 7) buffer SphereCount
{
	int sphereCount;
};

layout (std430, binding = 8) buffer antiAliasingIterations
{
	int aaIterations;
};


layout (std430, binding = 9) buffer Time
{
	float time;
};


layout(local_size_x = 300, local_size_y = 1, local_size_z = 1) in;

const vec3 lowerLeftCorner = vec3(-2.0, -1.5, -1.0);
const vec3 horizontal = vec3(4.0, 0.0, 0.0);
const vec3 vertical = vec3(0.0, 3.0, 0.0);
const vec3 origin = vec3(0.0, 0.0, 0.0);
const int sphereStructSize = 12;
int seedIncrement = 0;

// Some magic numbers found at http://lukas-polok.cz/tutorial_sphere.htm
const vec4 unitSphereMagic = vec4(1111.1111, 3141.5926, 2718.2818, 0);

#define ONE 0.00390625
#define ONEHALF 0.001953125
// Quote Stefan Gustafsson: 
// "The numbers above are 1/256 and 0.5/256, change accordingly
// if you change the code to use another perm/grad texture size."

float NextSeed()
{
	return time + seedIncrement++;
}

// Random function found at https://www.shadertoy.com/view/4djSRW
float rand(vec2 pos){
    return fract(sin(dot(pos + NextSeed(), vec2(12.9898, 78.233))) * 43758.5453);
}

// Fade function by Stefan Gustafsson
float fade(const in float t) {
  return t*t*t*(t*(t*6.0-15.0)+10.0);
}

// 3D Noise function by Stefan Gustafsson (Input parameters changed)
float noise(const in float x, const in float y, const in float z)
{
  vec3 P = vec3(x, y, z);

  vec3 Pi = ONE*floor(P)+ONEHALF; // Integer part, scaled so +1 moves one texel
                                  // and offset 1/2 texel to sample texel centers
  vec3 Pf = fract(P);     // Fractional part for interpolation

  // Noise contributions from (x=0, y=0), z=0 and z=1
  float perm00 = texture2D(permTexture, Pi.xy).a ;
  vec3  grad000 = texture2D(permTexture, vec2(perm00, Pi.z)).rgb * 4.0 - 1.0;
  float n000 = dot(grad000, Pf);
  vec3  grad001 = texture2D(permTexture, vec2(perm00, Pi.z + ONE)).rgb * 4.0 - 1.0;
  float n001 = dot(grad001, Pf - vec3(0.0, 0.0, 1.0));

  // Noise contributions from (x=0, y=1), z=0 and z=1
  float perm01 = texture2D(permTexture, Pi.xy + vec2(0.0, ONE)).a ;
  vec3  grad010 = texture2D(permTexture, vec2(perm01, Pi.z)).rgb * 4.0 - 1.0;
  float n010 = dot(grad010, Pf - vec3(0.0, 1.0, 0.0));
  vec3  grad011 = texture2D(permTexture, vec2(perm01, Pi.z + ONE)).rgb * 4.0 - 1.0;
  float n011 = dot(grad011, Pf - vec3(0.0, 1.0, 1.0));

  // Noise contributions from (x=1, y=0), z=0 and z=1
  float perm10 = texture2D(permTexture, Pi.xy + vec2(ONE, 0.0)).a ;
  vec3  grad100 = texture2D(permTexture, vec2(perm10, Pi.z)).rgb * 4.0 - 1.0;
  float n100 = dot(grad100, Pf - vec3(1.0, 0.0, 0.0));
  vec3  grad101 = texture2D(permTexture, vec2(perm10, Pi.z + ONE)).rgb * 4.0 - 1.0;
  float n101 = dot(grad101, Pf - vec3(1.0, 0.0, 1.0));

  // Noise contributions from (x=1, y=1), z=0 and z=1
  float perm11 = texture2D(permTexture, Pi.xy + vec2(ONE, ONE)).a ;
  vec3  grad110 = texture2D(permTexture, vec2(perm11, Pi.z)).rgb * 4.0 - 1.0;
  float n110 = dot(grad110, Pf - vec3(1.0, 1.0, 0.0));
  vec3  grad111 = texture2D(permTexture, vec2(perm11, Pi.z + ONE)).rgb * 4.0 - 1.0;
  float n111 = dot(grad111, Pf - vec3(1.0, 1.0, 1.0));

  // Blend contributions along x
  vec4 n_x = mix(vec4(n000, n001, n010, n011),
                 vec4(n100, n101, n110, n111), fade(Pf.x));

  // Blend contributions along y
  vec2 n_xy = mix(n_x.xy, n_x.zw, fade(Pf.y));

  // Blend contributions along z
  float n_xyz = mix(n_xy.x, n_xy.y, fade(Pf.z));

  // We're done, return the final noise value.
  return n_xyz;
}


// Random direction found at http://lukas-polok.cz/tutorial_sphere.htm
vec3 RandomInUnitSphere(vec2 pos)
{
	vec2 tc = pos.xy * unitSphereMagic.xy;

	vec3 skewed_seed = vec3(NextSeed() * unitSphereMagic.z + tc.y - tc.x) + unitSphereMagic.yzw;
	vec3 velocity;
	velocity.x = noise(tc.x, tc.y, skewed_seed.x);
	velocity.y = noise(tc.y, skewed_seed.y, tc.x);
	velocity.z = noise(skewed_seed.z, tc.x, tc.y);

	velocity = normalize(velocity);

	return velocity;
}

vec4 RandomInUnitDisk(vec2 pos)
{
	vec4 p = vec4(RandomInUnitSphere(pos).xy, 0.0, 0.0);
	return normalize(p);
}

Ray GetRay(in float s, in float t)
{
	vec4 rd = camera.lensRadius * RandomInUnitDisk(vec2(s,t));
	vec4 offset = camera.u * rd.x + camera.v * rd.y;
	return Ray(camera.origin + offset, camera.lowerLeftCorner + s*camera.horizontal + t*camera.vertical - camera.origin - offset); 
}

vec4 PointAtParameter(in float t, in vec4 A, in vec4 B)
{
	return A + B * t;
}

bool Hit(in vec4 sphereCenter, in float sphereRadius, in Ray ray, in float tMin, in float tMax, out HitRecord rec)
{	
	vec4 oc = ray.A - sphereCenter;
	float a = dot(ray.B, ray.B);
	float b = dot(oc, ray.B);
	float c = dot(oc, oc) - sphereRadius * sphereRadius;
	float discriminant = b * b - a * c;
	if(discriminant > 0)
	{
		float temp = (-b - sqrt(discriminant)) / a;
		if(temp < tMax && temp > tMin)
		{
			rec.t = temp;
			rec.p = PointAtParameter(rec.t, ray.A, ray.B);
			rec.normal = (rec.p - sphereCenter) / sphereRadius;
			return true;
		}
		temp = (-b + sqrt(discriminant)) / a;
		if(temp < tMax && temp > tMin)
		{
			rec.t = temp;
			rec.p = PointAtParameter(rec.t, ray.A, ray.B);
			rec.normal = (rec.p - sphereCenter) / sphereRadius;
			return true;
		}	
	}
	return false;
}

vec4 Reflect(in vec4 v, in vec4 n) {
	return v - 2 * dot(v, n) * n;
}

bool Refract(in vec4 v, in vec4 n, float ni_over_nt, out vec4 refracted)
{
	vec4 uv = normalize(v);
	float dt = dot(uv, n);
	float discriminant = 1.0 - ni_over_nt * ni_over_nt * (1 - dt * dt);
	if(discriminant > 0)
	{
		refracted = ni_over_nt * (uv - n * dt) - n * sqrt(discriminant);
		return true;
	}
	false;
}

float Schlick(float cosine, float reflectionIndex) {
	float r0 = (1.0 - reflectionIndex) / (1.0 + reflectionIndex);
	r0 *= r0;
	return r0 + (1.0 - r0) * pow(1.0 - cosine, 5);
}

bool Scatter(in Ray inRay, in HitRecord rec, in int matType, in vec4 albedo, in float fuzz, in float reflectionIndex, out vec4 attenuation, out Ray scattered)
{
	if(matType == 0)
	{
		vec4 target = rec.p + rec.normal + vec4(RandomInUnitSphere(point_input[gl_GlobalInvocationID.x].xy), 0.0);
		scattered = Ray(rec.p, target - rec.p);
		attenuation = albedo;
		return true;
	}
	else if(matType == 1)
	{
		vec4 reflected = Reflect(normalize(inRay.B), rec.normal);
		scattered = Ray(rec.p, reflected + fuzz * vec4(RandomInUnitSphere(point_input[gl_GlobalInvocationID.x].xy), 0.0));
		attenuation = albedo;
		return dot(scattered.B, rec.normal) > 0;
	}
	else if(matType == 2)
	{
		vec4 outwardNormal;
		vec4 reflected = Reflect(inRay.B, rec.normal);
		float ni_over_nt;
		attenuation = vec4(1.0, 1.0, 1.0, 0.0);
		vec4 refracted;
		float reflectProb;
		float cosine;

		if(dot(inRay.B, rec.normal) > 0)
		{
			outwardNormal = -rec.normal;
			ni_over_nt = reflectionIndex;
			cosine = reflectionIndex * dot(inRay.B, rec.normal) / length(inRay.B);
		}
		else
		{
			outwardNormal = rec.normal;
			ni_over_nt = 1.0 / reflectionIndex;
			cosine = -dot(inRay.B, rec.normal) / length(inRay.B);
		}

		if(Refract(inRay.B, outwardNormal, ni_over_nt, refracted))
		{
			reflectProb = Schlick(cosine, reflectionIndex);
		}
		else {
			reflectProb = 1.0;
		}

		float random = rand(point_input[gl_GlobalInvocationID.x].xy);

		if(random < reflectProb) {
			scattered = Ray(rec.p, reflected);
		}
		else {
			scattered = Ray(rec.p, refracted);
		}
		return true;
	}
}

void GetClosestHit(in Ray ray, out int hitIndex, out bool hitAnything, out HitRecord rec)
{
	const float tMin = 0.00001;
	float closestSoFar = 10000000;
	hitAnything = false;

	HitRecord tempRec;
	for(int i = 0; i < sphereCount; i++)
	{
		vec4 sphereCenter = sphere_input[i].center;
		float sphereRadius = sphereCenter.w;
		sphereCenter.w = 0.0;

		if(Hit(sphereCenter, sphereRadius, ray, tMin, closestSoFar, tempRec))
		{	
			hitAnything = true;
			closestSoFar = tempRec.t;
			hitIndex = i;
			rec = tempRec;
		}	
	}
}


vec4 Color(Ray ray)
{
	HitRecord rec;
	bool hitAnything = false;
	int hitIndex = -1;

	GetClosestHit(ray, hitIndex, hitAnything, rec);


	if(hitAnything)
	{
		Ray scattered;
		vec4 attenuation;
		int matType = sphere_input[hitIndex].matType;
		vec4 albedo = sphere_input[hitIndex].albedo;
		float fuzz = sphere_input[hitIndex].fuzz;
		float reflectionIndex = sphere_input[hitIndex].reflectionIndex;

		vec4 returnColor = vec4(1,1,1,0.0);
		for(int i = 0; i < 10; i++)
		{
			if(hitAnything && Scatter(ray, rec, matType, albedo, fuzz, reflectionIndex, attenuation, scattered))
			{
				returnColor *= attenuation;

				ray = scattered;
				GetClosestHit(ray, hitIndex, hitAnything, rec);
				if(hitAnything)
				{
					matType = sphere_input[hitIndex].matType;
					albedo = sphere_input[hitIndex].albedo;
					fuzz = sphere_input[hitIndex].fuzz;
					reflectionIndex = sphere_input[hitIndex].reflectionIndex;
				}
			}
			else
			{
				
				vec4 unitDirection = normalize(ray.B);
				float lerpt = 0.5 * (unitDirection.y + 1.0);
				vec4 c = (1.0 - lerpt) * vec4(1.0, 1.0, 1.0, 0.0) + lerpt * vec4(0.5, 0.7, 1.0, 0.0);
				returnColor *= c;
				break;
			}
		}
		return returnColor;
	}
	else
	{
		vec4 unitDirection = normalize(ray.B);
		float lerpt = 0.5 * (unitDirection.y + 1.0);
		return (1.0 - lerpt) * vec4(1.0, 1.0, 1.0, 0.0) + lerpt * vec4(0.5, 0.7, 1.0, 0.0);
	}
}

void main()
{
	uint id = gl_GlobalInvocationID.x;

	float posx = point_input[id].x;
	float posy = point_input[id].y;
	float parameter = point_input[id].z;
	

	vec4 color = vec4(0, 0, 0, 0);
	for(int i = 0; i < aaIterations; i++)
	{
		float u = (posx + (rand(point_input[id].xy) * 2 - 1)) / 800.0;
		float v = (posy + (rand(point_input[id].xy) * 2 - 1)) / 600.0;
		Ray ray = GetRay(u, v);

		color += Color(ray);
	}
	color /= aaIterations;

	color_output[id] = color;
})";

#pragma endregion 




	void SetupCamera(Vec3 lookFrom)
	{
		Vec3 lookAt(0, 0, 0);
		float distToFocus = (lookFrom - lookAt).Magnitude();
		float aperture = 0.1f;

		cam = Camera(lookFrom, lookAt, Vec3(0, 1, 0), 20, static_cast<float>(WINDOW_WIDTH) / WINDOW_HEIGHT, aperture, distToFocus);
	}

void handle_key_event(GLFWwindow* window, int key, int scancode, int action, int modifiers)
{
	if(action != GLFW_PRESS)
		return;

	if(key == GLFW_KEY_ESCAPE)
		glfwDestroyWindow(window);
	else if (key == GLFW_KEY_D || key == GLFW_KEY_A || GLFW_KEY_W || GLFW_KEY_S)
	{
		Vec3 newLookFrom = cam.origin;
		if (key == GLFW_KEY_D)
		{
			newLookFrom += Vec3::Cross(-cam.origin, Vec3(0,1,0)).Normalized();
		}
		else if (key == GLFW_KEY_A)
		{
			newLookFrom += Vec3::Cross(cam.origin, Vec3(0,1,0)).Normalized();
		}
		else if(key == GLFW_KEY_W)
		{
			newLookFrom += -cam.origin.Normalized();
		}
		else if (key == GLFW_KEY_S)
		{
			newLookFrom += cam.origin.Normalized();
		}

		SetupCamera(newLookFrom);
		camerastruct camstruct{
		cam.origin,
		cam.lowerLeftCorner,
		cam.horizontal,
		cam.vertical,
		cam.u,
		cam.v,
		cam.w,
		cam.lensRadius
		};


		glGenBuffers(1, &camSSbo);
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, camSSbo);
		glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(camerastruct), &camstruct, GL_STATIC_DRAW);
	}
	else
		printf("key_event: %d\n", key);
}

void handle_mouse_event(GLFWwindow* window, int button, int action, int modifiers)
{
	if(action != GLFW_PRESS)
		return;

	if(button == GLFW_MOUSE_BUTTON_LEFT)
	{
		double mouse_x, mouse_y;
		glfwGetCursorPos(window, &mouse_x, &mouse_y);
		printf("bom (%f, %f)\n", mouse_x, mouse_y);
	}
}





Vec3 Color(const Ray& r, HittableList* world, int depth) {
	HitRecord rec;
	constexpr float floatMax = std::numeric_limits<float>::max();

	if (world->Hit(r, 0.001f, floatMax, rec))
	{
		Ray scattered;
		Vec3 attenuation;
		if (depth < 50 && rec.mat_ptr->Scatter(r, rec, attenuation, scattered)) {
			return attenuation * Color(scattered, world, depth + 1);
		}
		else {
			return Vec3(0,0,0);
		}
	}
	else
	{
		Vec3 unitDirection = r.Direction().Normalized();
		float t = 0.5f * (unitDirection.y + 1.0f);
		return (1.0f - t) * Vec3(1.0f, 1.0f, 1.0f) + t * Vec3(0.5f, 0.7f, 1.0f);
	}
}

HittableList* RandomScene() {
	const int n = 20;
	Hittable** list = new Hittable*[n];
	list[0] = new Sphere(Vec3(0, -1000, 0), 1000, new Lambertian(Vec3(0.5f, 0.5f, 0.5f)));
	int i = 1;

	for (int a = -2; a < 2; a++)
	{
		for(int b = -2; b < 2; b++)
		{
			float choose_mat = RandNext();
			Vec3 center(a+0.9f * RandNext(), 0.2f, b + 0.9f * RandNext());
			if ((center - Vec3(4.f, 0.2f, 0.f)).Magnitude() > 0.9f) {
				if (choose_mat < 0.8f) { // Diffuse
					list[i++] = new Sphere(center, 0.2f, new Lambertian(Vec3(RandNext() * RandNext(), RandNext() * RandNext(), RandNext() * RandNext())));
				}
				else if (choose_mat < 0.95f) { // Metal
					list[i++] = new Sphere(center, 0.2f, new Metal(Vec3(0.5f * (1 + RandNext()), 0.5f * (1 + RandNext()), 0.5f * (1 + RandNext())), 0.5f * RandNext()));
				}
				else {// Glass
					list[i++] = new Sphere(center, 0.2f, new Dielectric(1.5f));
				}
			}
		}
	}

	list[i++] = new Sphere(Vec3(0.f, 1.f, 0.f), 1.f, new Dielectric(1.5f));
	list[i++] = new Sphere(Vec3(-4, 1, 0), 1.f, new Lambertian(Vec3(0.4f, 0.2f, 0.1f)));
	list[i++] = new Sphere(Vec3(4, 1, 0), 1.f, new Metal(Vec3(0.7f, 0.6f, 0.5f), 0.0f));

	return new HittableList(list, i);
}

int main()
{

	glfwInit();
	GLFWwindow* window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "Graphics are cool", NULL, NULL);
	glfwMakeContextCurrent(window);

	glfwSetKeyCallback(window, handle_key_event);
	glfwSetMouseButtonCallback(window, handle_mouse_event);


	//glGenBuffers(1, NULL);

	glewInit();

	glViewport(0, 0, WINDOW_WIDTH, WINDOW_HEIGHT);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, WINDOW_WIDTH, 0, WINDOW_HEIGHT, -1, 1);

	#pragma region ComputeShaderTest
	GLuint compShader = glCreateShader(GL_COMPUTE_SHADER);
	glShaderSource(compShader, 1, &CALC_SRC, NULL);
	glCompileShader(compShader);

	static char LOG_BUFFER[1024];
	glGetShaderInfoLog(compShader, 1024, NULL, LOG_BUFFER);

	printf(LOG_BUFFER);


	GLuint program = glCreateProgram();
	glAttachShader(program, compShader);
	glLinkProgram(program);
	location_permTexture = glGetUniformLocation(program, "permTexture");

	glUseProgram(program);

	initPermTexture(&permTextureID);


	const int nx = WINDOW_WIDTH;
	const int ny = WINDOW_HEIGHT;

	SetupCamera(Vec3(3, 2, 13));

	camerastruct camstruct{
		cam.origin,
		cam.lowerLeftCorner,
		cam.horizontal,
		cam.vertical,
		cam.u,
		cam.v,
		cam.w,
		cam.lensRadius
	};
	

	glGenBuffers(1, &camSSbo);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, camSSbo);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(camerastruct), &camstruct, GL_STATIC_DRAW);

	computestruct* points = new computestruct[PIXEL_COUNT];

	for (int j = WINDOW_HEIGHT - 1; j >= 0; j--)
	{
		for (int i = 0; i < WINDOW_WIDTH; i++)
		{
			int idx = j * WINDOW_WIDTH + i;
			points[idx].u = i;
			points[idx].v = j;
			points[idx].parameter = 2.0f;
			points[idx].t = 2.f;
		}
	}


	GLuint posSSbo;
	glGenBuffers(1, &posSSbo);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, posSSbo);
	glBufferData(GL_SHADER_STORAGE_BUFFER, PIXEL_COUNT * sizeof(computestruct), points, GL_STATIC_DRAW);



	colorstruct* colors = new colorstruct[PIXEL_COUNT];


	GLuint colSSbo;
	glGenBuffers(1, &colSSbo);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, colSSbo);
	glBufferData(GL_SHADER_STORAGE_BUFFER, PIXEL_COUNT * sizeof(colorstruct), colors, GL_STATIC_DRAW);

	HittableList* world = static_cast<HittableList*>(RandomScene());

	int sphereCount = world->listSize;
	spherestruct* spheres = new spherestruct[sphereCount];

	for (int i = 0; i < sphereCount; i++)
	{
		Sphere* sphere = static_cast<Sphere*>(world->list[i]);
		spheres[i].center = sphere->center;
	//	printf("%f, %f, %f\n", spheres[i].center.x, spheres[i].center.y, spheres[i].center.z);
		spheres[i].center.w = sphere->radius;

		//printf("%f\n", spheres[i].center.w);

		if (Lambertian* lamb = dynamic_cast<Lambertian*>(sphere->material))
		{
			spheres[i].albedo = lamb->albedo;
			spheres[i].matType = 0;
		}
		else if (Metal* met = dynamic_cast<Metal*>(sphere->material))
		{
			spheres[i].albedo = met->albedo;
			spheres[i].fuzz = met->fuzz;
			spheres[i].matType = 1;
		}
		else if (Dielectric* die = dynamic_cast<Dielectric*>(sphere->material))
		{
			spheres[i].reflectionIndex = die->reflectionIndex;
			spheres[i].matType = 2;
		}
	}


	GLuint sphereSSbo;
	glGenBuffers(1, &sphereSSbo);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 6, sphereSSbo);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(spherestruct) * sphereCount, spheres, GL_STATIC_DRAW);

	//
	GLuint sphereCountSSbo;
	glGenBuffers(1, &sphereCountSSbo);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 7, sphereCountSSbo);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(GLint), &sphereCount, GL_STATIC_DRAW);


	GLuint aaIterationsSSbo;
	glGenBuffers(1, &aaIterationsSSbo);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 8, aaIterationsSSbo);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(GLint), &antiAliasingIterations, GL_STATIC_DRAW);

	float time = glfwGetTime();

	GLuint timeSSbo;
	glGenBuffers(1, &timeSSbo);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 9, timeSSbo);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(GLfloat), &time, GL_STATIC_DRAW);

	


	




	//printf("%f, %f, %f, %f", colors[0].r, colors[0].g, colors[0].b, colors[0].a);

#pragma endregion // Compute Shader test

	//const int nx = WINDOW_WIDTH;
	//const int ny = WINDOW_HEIGHT;
	//const int ns = 15;

	//float radius = cos(F_PI * 0.25f);

	///*Hittable* hittables[4];
	//hittables[0] = new Sphere(Vec3(0.f, 0.f, -1.f), radius, new Lambertian(Vec3(0.1f, 0.2f, 0.5f)));
	//hittables[1] = new Sphere(Vec3(0.f, -100 - radius, -1.f), 100, new Lambertian(Vec3(0.8f, 0.8f, 0.f)));
	//hittables[2] = new Sphere(Vec3(radius * 2.f, 0.f, -1.f), radius, new Metal(Vec3(0.8f, 0.6f, 0.2f), 0.3f));
	//hittables[3] = new Sphere(Vec3(-radius * 2.f, 0.f, -1.f), radius, new Dielectric(1.5f));

	//HittableList* world = new HittableList(hittables, 4);*/

	//Vec3 lookFrom(13, 2, 3);
	//Vec3 lookAt(0, 0, 0);
	//float distToFocus = (lookFrom - lookAt).Magnitude();
	//float aperture = 0.1f;

	//Camera cam(lookFrom, lookAt, Vec3(0, 1, 0), 20, float(nx) / float(ny), aperture, distToFocus);

	

	while (!glfwWindowShouldClose(window))
	{
		glClearColor(0.1f, 0.1f, 0.8f, 1.f);
		glClear(GL_COLOR_BUFFER_BIT);
		glClear(GL_DEPTH_BUFFER_BIT);

		/*glLoadIdentity();
		gluOrtho2D(0.0, 500.0, 500.0, 0.0);*/

		printf("begin");
		glDispatchCompute(PIXEL_COUNT / WORK_GROUP_SIZE, 1, 1);
		glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);



		glBindBuffer(GL_SHADER_STORAGE_BUFFER, colSSbo);
		colorstruct* ptr = (colorstruct*)glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_WRITE_ONLY);
		for (int i = 0; i < PIXEL_COUNT; i++)
		{
			//printf("%f, %f, %f", ptr[i].r, ptr[i].g, ptr[i].b);
			colors[i] = ptr[i];
		}
		glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);

		glBegin(GL_POINTS);

		#pragma region Compute Shader
		for (int j = ny - 1; j >= 0; j--)
		{
			for (int i = 0; i < nx; i++)
			{
				int idx = j * nx + i;

				//printf("%f, %f, %f\n", colors[idx].r, colors[idx].g, colors[idx].b);

				glColor3f(colors[idx].r, colors[idx].g, colors[idx].b);
				glVertex2i(i, j);
			}
		}
		#pragma endregion // Compute Shader Test

		//for (int j = ny - 1; j >= 0; j--)
		//{
		//	printf("row done %d\n", j);
		//	for (int i = 0; i < nx; i++)
		//	{
		//		Vec3 col(0.f, 0.f, 0.f);
		//		for (int s = 0; s < antiAliasingIterations; s++)
		//		{
		//			float r1 = RandNext();
		//			float r2 = RandNext();

		//			const float u = static_cast<float>(i + r1) / static_cast<float>(nx);
		//			const float v = static_cast<float>(j + r2) / static_cast<float>(ny);
		//			const Ray ray = cam.GetRay(u, v);
		//			//const Vec3 p = ray.PointAtParameter(2.0f);
		//			col += Color(ray, world, 0);
		//		}
		//		col /= float(antiAliasingIterations);
		//		col = Vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));
		//		
		//		glColor3f(col.r, col.g, col.b);
		//		glVertex2i(i, j);
		//	}
		//}

		printf("nice");

		glEnd();

		glfwSwapBuffers(window);

		glfwPollEvents();
	}

	//delete[] colors;
	//delete points;
	/*delete world;
	delete[] hittables;*/
}