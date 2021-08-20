#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <GL/GL.h>
#include "Math/MathUtility.h"
#include "Math/Vec3.h"
#include "Math/Vec4.h"
#include "Math/Ray.h"
#include "Geometry/Hittable.h"
#include "Geometry/Sphere.h"
#include "Framework/HittableList.h"
#include "Rendering/Camera.h"


/*
	Variables to control program
*/

const int WINDOW_WIDTH = 1200;
const int WINDOW_HEIGHT = 800;
const int PIXEL_COUNT = WINDOW_HEIGHT * WINDOW_WIDTH;
const int WORK_GROUP_SIZE = 300;
const int antiAliasingIterations = 30;
const Vec3 cameraLookAt = Vec3(0,0,0);
const Vec3 initialCameraPosition = Vec3(18, 4, 18);
#define USE_COMPUTE_SHADER


/*
	Global convenience variables
*/
GLuint camSSbo;
Camera cam;
bool wHeld;
bool aHeld;
bool sHeld;
bool dHeld;

/*
	Structs used to send data to compute shader
*/

struct WindowCoords {
	int x, y;
};

struct SphereData {
	Vec4 center;
	Vec4 albedo;
	float fuzz;
	float reflectionIndex;
	int matType;
	float padd;
};

struct CameraData {
	Vec4 origin;
	Vec4 lowerLeftCorner;
	Vec4 horizontal;
	Vec4 vertical;
	Vec4 u, v, w;
	float lensRadius;
};

/*
	The compute shader in its entirety
*/

#pragma region Compute Shader SRC
const char* CALC_SRC = R"(#version 430 compatibility

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
};

layout(rgba32f, binding = 0) uniform image2D imgOutput;

layout(std430, binding = 3) buffer cam 
{
	Camera camera;
};

layout(std430, binding = 4) buffer Dat
{
	ivec2 point_input[];
};

layout (std430, binding = 6) buffer Spheres
{
	SphereData sphere_input[];
};

uniform int sphereCount;
uniform int aaIterations;
uniform int windowWidth;
uniform int windowHeight;
uniform float time;

layout(local_size_x = 300, local_size_y = 1, local_size_z = 1) in;

#define PI 3.1415926538;
int seedIncrement = 0;

float NextSeed()
{
	seedIncrement += 5; // Arbritrary magic number preferably a bit larger than 1
	return time + seedIncrement;
}

// Random function found at https://www.shadertoy.com/view/4djSRW
float rand(vec2 pos){
    return fract(sin(dot(pos + NextSeed(), vec2(12.9898, 78.233))) * 43758.5453);
}

vec3 RandomInUnitSphere(vec2 pos)
{
	float z = (rand(pos) - 0.5) * 1.99;
	float rxy = sqrt(1 - z*z);
	float phi = rand(pos) * 2 * PI
	float x = rxy * cos(phi);
	float y = rxy * sin(phi);

	return normalize(vec3(x, y, z));
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
	float halfB = dot(oc, ray.B);
	float c = dot(oc, oc) - sphereRadius * sphereRadius;

	float discriminant = halfB * halfB - a * c;

	if(discriminant < 0) return false;

	float sqrtd = sqrt(discriminant);

	float root = (-halfB - sqrtd) / a;
	if(root < tMin || root > tMax) {
		root = (-halfB + sqrtd) / a;
		if(root < tMin || root > tMax)
			return false;
	}

	rec.t = root;
	rec.p = PointAtParameter(rec.t, ray.A, ray.B);
	rec.normal = (rec.p - sphereCenter) / sphereRadius;

	return true;
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

vec2 GetScreenPos()
{
	return point_input[gl_GlobalInvocationID.x];
}

bool Scatter(in Ray inRay, in HitRecord rec, in int matType, in vec4 albedo, in float fuzz, in float reflectionIndex, out vec4 attenuation, out Ray scattered)
{
	if(matType == 0)
	{
		vec3 randomDir = RandomInUnitSphere(GetScreenPos());
		vec4 target = rec.p + rec.normal + vec4(randomDir, 0.0);
		scattered = Ray(rec.p, target - rec.p);
		attenuation = albedo;
		return true;
	}
	else if(matType == 1)
	{
		vec4 reflected = Reflect(normalize(inRay.B), rec.normal);
		scattered = Ray(rec.p, reflected + fuzz * vec4(RandomInUnitSphere(GetScreenPos()), 0.0));
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

		float random = rand(GetScreenPos());

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
	const float tMin = 0.0001;
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

		vec4 returnColor = vec4(1,1,1,0);
		for(int i = 0; i < 50; i++)
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

	vec2 pos = point_input[id];
	

	vec4 color = vec4(0, 0, 0, 0);
	for(int i = 0; i < aaIterations; i++)
	{
		float u = (pos.x + (rand(pos) * 2 - 1)) / windowWidth;
		float v = (pos.y + (rand(pos) * 2 - 1)) / windowHeight;
		Ray ray = GetRay(u, v);

		color += Color(ray);
	}
	color /= aaIterations;

	ivec2 px = ivec2(pos);
	color.w = 1.0;
	imageStore(imgOutput, px, color);
})";

#pragma endregion 




void SetupCamera(Vec3 lookFrom)
{
	float distToFocus = (lookFrom - cameraLookAt).Magnitude();
	float aperture = 0.1f;

	cam = Camera(lookFrom, cameraLookAt, Vec3(0, 1, 0), 20, static_cast<float>(WINDOW_WIDTH) / WINDOW_HEIGHT, aperture, distToFocus);
}

void SendCameraDataToComputeShader() {
	CameraData camstruct{
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
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(CameraData), &camstruct, GL_STATIC_DRAW);
}

void handle_key_event(GLFWwindow* window, int key, int scancode, int action, int modifiers)
{

	if(key == GLFW_KEY_ESCAPE)
		glfwDestroyWindow(window);

	if (action == GLFW_PRESS)
	{
		if(key == GLFW_KEY_W)
			wHeld = true;
		else if(key == GLFW_KEY_S)
			sHeld = true;
		else if(key == GLFW_KEY_A)
			aHeld = true;
		else if(key == GLFW_KEY_D)
			dHeld = true;
	}
	else if(action == GLFW_RELEASE)
	{
		if(key == GLFW_KEY_W)
			wHeld = false;
		else if(key == GLFW_KEY_S)
			sHeld = false;
		else if(key == GLFW_KEY_A)
			aHeld = false;
		else if(key == GLFW_KEY_D)
			dHeld = false;
	}	
}


void UpdateInput()
{
	const float cameraSpeed = 100.f;
	if (wHeld || sHeld || aHeld || dHeld)
	{
		Vec3 newLookFrom = cam.origin;
		if (wHeld)
			newLookFrom += -cam.origin.Normalized();
		else if (sHeld)
			newLookFrom += cam.origin.Normalized();
		else if (aHeld)
			newLookFrom += Vec3::Cross(cam.origin, Vec3(0, 1, 0)).Normalized();
		else if (dHeld)
			newLookFrom += Vec3::Cross(-cam.origin, Vec3(0, 1, 0)).Normalized();

		SetupCamera(newLookFrom);
		SendCameraDataToComputeShader();
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
				else { // Glass
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
	GLFWwindow* window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "Path Tracer", NULL, NULL);
	glfwMakeContextCurrent(window);

	glfwSetKeyCallback(window, handle_key_event);

	glewInit();
	glViewport(0, 0, WINDOW_WIDTH, WINDOW_HEIGHT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, WINDOW_WIDTH, 0, WINDOW_HEIGHT, -1, 1);

	// Generate Spheres
	HittableList* world = static_cast<HittableList*>(RandomScene());

	SetupCamera(initialCameraPosition);

	#ifdef USE_COMPUTE_SHADER 


	//////////////////////// Setup shaders ////////////////////////
	GLuint compShader = glCreateShader(GL_COMPUTE_SHADER);
	glShaderSource(compShader, 1, &CALC_SRC, NULL);
	glCompileShader(compShader);

	static char LOG_BUFFER[1024];
	glGetShaderInfoLog(compShader, 1024, NULL, LOG_BUFFER);

	printf(LOG_BUFFER);


	//////////////////////// Setup program ////////////////////////
	GLuint program = glCreateProgram();
	glAttachShader(program, compShader);
	glLinkProgram(program);

	GLint u_Time = glGetUniformLocation(program, "time");
	GLint u_SphereCount = glGetUniformLocation(program, "sphereCount");
	GLint u_AAIterations = glGetUniformLocation(program, "aaIterations");
	GLint u_WindowWidth = glGetUniformLocation(program, "windowWidth");
	GLint u_WindowHeight = glGetUniformLocation(program, "windowHeight");

	glUseProgram(program);


	//////////////////////// Setup output texture //////////////////////// 
	GLuint texOutput;
	glGenTextures(1, &texOutput);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, texOutput);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, WINDOW_WIDTH, WINDOW_HEIGHT, 0, GL_RGBA, GL_FLOAT,
		NULL);
	glBindImageTexture(0, texOutput, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);

	
	glUniform1i(u_WindowWidth, WINDOW_WIDTH);
	glUniform1i(u_WindowHeight, WINDOW_HEIGHT);


	//////////////////////// Send camera data ////////////////////////
	SendCameraDataToComputeShader();


	//////////////////////// Fill a buffer with window coordinates ////////////////////////
	WindowCoords* points = new WindowCoords[PIXEL_COUNT];

	for (int j = WINDOW_HEIGHT - 1; j >= 0; j--)
	{
		for (int i = 0; i < WINDOW_WIDTH; i++)
		{
			int idx = j * WINDOW_WIDTH + i;
			points[idx].x = i;
			points[idx].y = j;
		}
	}
	GLuint posSSbo;
	glGenBuffers(1, &posSSbo);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, posSSbo);
	glBufferData(GL_SHADER_STORAGE_BUFFER, PIXEL_COUNT * sizeof(WindowCoords), points, GL_STATIC_DRAW);


	

	//////////////////////// Gather sphere data for compute shader ////////////////////////
	int sphereCount = world->listSize;
	SphereData* sphereDatas = new SphereData[sphereCount];

	for (int i = 0; i < sphereCount; i++)
	{
		Sphere* sphere = static_cast<Sphere*>(world->list[i]);
		sphereDatas[i].center = sphere->center;
		sphereDatas[i].center.w = sphere->radius;

		if (Lambertian* lamb = dynamic_cast<Lambertian*>(sphere->material))
		{
			sphereDatas[i].albedo = lamb->albedo;
			sphereDatas[i].matType = 0;
		}
		else if (Metal* met = dynamic_cast<Metal*>(sphere->material))
		{
			sphereDatas[i].albedo = met->albedo;
			sphereDatas[i].fuzz = met->fuzz;
			sphereDatas[i].matType = 1;
		}
		else if (Dielectric* die = dynamic_cast<Dielectric*>(sphere->material))
		{
			sphereDatas[i].reflectionIndex = die->reflectionIndex;
			sphereDatas[i].matType = 2;
		}
	}
	GLuint sphereSSbo;
	glGenBuffers(1, &sphereSSbo);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 6, sphereSSbo);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(SphereData) * sphereCount, sphereDatas, GL_STATIC_DRAW);


	glUniform1i(u_SphereCount, sphereCount);
	glUniform1i(u_AAIterations, antiAliasingIterations);

	#endif
	
	
	//////////////////////// Render Loop ////////////////////////
	while (!glfwWindowShouldClose(window))
	{
		double time = glfwGetTime();

		glClearColor(0.1f, 0.1f, 0.8f, 1.f);
		glClear(GL_COLOR_BUFFER_BIT);
		

#ifdef USE_COMPUTE_SHADER
		glUseProgram(program);
		glUniform1f(u_Time, static_cast<float>(time));

		glDispatchCompute(PIXEL_COUNT / WORK_GROUP_SIZE, 1, 1);
		glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

		


		glBindTexture(GL_TEXTURE_2D, texOutput);
		glEnable(GL_TEXTURE_2D);
		glBegin(GL_QUADS);
		glTexCoord2i(0, 1); glVertex2i(0, WINDOW_HEIGHT);
		glTexCoord2i(0, 0); glVertex2i(0, 0);
		glTexCoord2i(1, 0); glVertex2i(WINDOW_WIDTH, 0);
		glTexCoord2i(1, 1); glVertex2i(WINDOW_WIDTH, WINDOW_HEIGHT);
		glEnd();
		glDisable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 0);
		glFlush();
#else

		glBegin(GL_POINTS);
		for (int j = WINDOW_HEIGHT - 1; j >= 0; j--)
		{
			printf("row done %d\n", j);
			for (int i = 0; i < WINDOW_WIDTH; i++)
			{
				Vec3 col(0.f, 0.f, 0.f);
				for (int s = 0; s < antiAliasingIterations; s++)
				{
					float r1 = RandNext();
					float r2 = RandNext();

					const float u = static_cast<float>(i + r1) / WINDOW_WIDTH;
					const float v = static_cast<float>(j + r2) / WINDOW_HEIGHT;
					const Ray ray = cam.GetRay(u, v);

					col += Color(ray, world, 0);
				}
				col /= static_cast<float>(antiAliasingIterations);
				
				glColor3f(col.r, col.g, col.b);
				glVertex2i(i, j);
			}
		}
		glEnd();

#endif

		glfwSwapBuffers(window);
		double delta = glfwGetTime() - time;
		glfwPollEvents();

		UpdateInput();

		double fps = 1.f / delta;
		// Input causes severe deltatime spikes. Limit prints to 120 fps.
		if(fps < 120.f)
			printf("Fps: %f\n", 1.f / delta);
	}

	////////////////////////  Cleanup //////////////////////// 
	delete world;
	
	#ifdef USE_COMPUTE_SHADER
	delete[] sphereDatas;
	#endif
}