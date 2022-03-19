#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <GL/GL.h>
#include "Math/MathUtility.h"
#include "Math/Vec3.h"
#include "Math/Vec4.h"
#include "Geometry/Sphere.h"
#include "Rendering/Camera.h"
#include <vector>
#include <chrono>
#define _USE_MATH_DEFINES
#include <math.h>
#include <algorithm>


/*
	Variables to control program
*/
const int WINDOW_WIDTH = 1200;
const int WINDOW_HEIGHT = 800;
const int PIXEL_COUNT = WINDOW_HEIGHT * WINDOW_WIDTH;
const int WORK_GROUP_SIZE = 300;
const int AA_ITERATIONS = 20;
const int SPHERE_COUNT = 20; // One sphere is reserverd for the ground
const Vec3 CAMERA_LOOK_AT = Vec3(0,0.5f,0);
const Vec3 INITIAL_CAMERA_POSITION = Vec3(18, 4, 18);

const float CAMERA_ANGULAR_SPEED_HORIZONTAL = static_cast<float>(M_PI / 2);
const float CAMERA_ANGULAR_SPEED_VERTICAL = static_cast<float>(M_PI / 4);
const float CAMERA_FORWARD_MULTIPLIER = 0.2f; // Value between 0 and 1;


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
	float reflection_index;
	int mat_type;
	float padd;
};

struct CameraData {
	Vec4 origin;
	Vec4 lower_left_corner;
	Vec4 horizontal;
	Vec4 vertical;
	Vec4 u, v, w;
	float lens_radius;
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




void SetupCamera(Camera& cam, const Vec3& look_from)
{
	float dist_to_focus = (look_from - CAMERA_LOOK_AT).Magnitude();
	float aperture = 0.1f;

	cam.Set(look_from, CAMERA_LOOK_AT, Vec3(0, 1, 0), 20, static_cast<float>(WINDOW_WIDTH) / WINDOW_HEIGHT, aperture, dist_to_focus);
}

void SendCameraDataToComputeShader(const Camera& cam) {
	CameraData camstruct{
		cam.origin,
		cam.lowerLeftCorner,
		cam.horizontal,
		cam.vertical,
		cam.u,
		cam.v,
		cam.w,
		cam.lens_radius
	};

	GLuint camSSbo;
	glGenBuffers(1, &camSSbo);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, camSSbo);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(CameraData), &camstruct, GL_STATIC_DRAW);
}

// Camera controls
void UpdateCamera(GLFWwindow* window, Camera& cam, const float delta_time)
{
	Vec3 position = cam.origin;
	Vec3 look_at = CAMERA_LOOK_AT;

	Vec3 vec_from_lookat = position - look_at;

	//// Forward/backward movement
	const float min_distance = 5.f;
	float dist_from_lookat = vec_from_lookat.Magnitude();
	float forward_multiplier = std::clamp(CAMERA_FORWARD_MULTIPLIER * delta_time, 0.f, 0.9f);

	if (min_distance < dist_from_lookat && glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
		position -= vec_from_lookat.Normalized() * (1 - CAMERA_FORWARD_MULTIPLIER * delta_time);
	if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
		position += vec_from_lookat.Normalized() * (1 + CAMERA_FORWARD_MULTIPLIER * delta_time);

	
	
	// Left/right movement
	vec_from_lookat = position - look_at;
	float flat_angle = atan2(vec_from_lookat.z, vec_from_lookat.x);
	const float flat_dist_from_lookat = sqrt(vec_from_lookat.x * vec_from_lookat.x + vec_from_lookat.z * vec_from_lookat.z);

	if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
		flat_angle -= CAMERA_ANGULAR_SPEED_HORIZONTAL * delta_time;
	if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
		flat_angle += CAMERA_ANGULAR_SPEED_HORIZONTAL * delta_time;

	position.x = cos(flat_angle) * flat_dist_from_lookat + look_at.x;
	position.z = sin(flat_angle) * flat_dist_from_lookat + look_at.z;

	// Up/down movement
	vec_from_lookat = position - look_at;
	const float vertical_angle = atan2(vec_from_lookat.y, flat_dist_from_lookat);
	float new_vertical_angle = vertical_angle;
	if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
		new_vertical_angle += CAMERA_ANGULAR_SPEED_VERTICAL * delta_time;
	if (glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS)
		new_vertical_angle -= CAMERA_ANGULAR_SPEED_VERTICAL * delta_time;

	new_vertical_angle = std::clamp(new_vertical_angle, 0.f, 0.9f);

	float test_new_dist = sqrt(vec_from_lookat.x * vec_from_lookat.x + vec_from_lookat.z * vec_from_lookat.z);
	dist_from_lookat = vec_from_lookat.Magnitude();
	position.y = sin(new_vertical_angle) * dist_from_lookat + look_at.y;
	const float new_flat_dist = abs(cos(new_vertical_angle)) * dist_from_lookat;
	position.x *= new_flat_dist / flat_dist_from_lookat;
	position.z *= new_flat_dist / flat_dist_from_lookat;

	SetupCamera(cam, position);
	SendCameraDataToComputeShader(cam);
}

Material GenerateRandomMaterial()
{
	float choose_mat = RandNext();
	if (choose_mat < 0.8f) { // Diffuse
		return Material(Lambertian(Vec3(RandNext() * RandNext(), RandNext() * RandNext(), RandNext() * RandNext())));
	}
	else if (choose_mat < 0.95f) { // Metal
		return Material(Metal(Vec3(0.5f * (1 + RandNext()), 0.5f * (1 + RandNext()), 0.5f * (1 + RandNext())), 0.5f * RandNext()));
	}
	else { // Glass
		return Material(Dielectric(1.5f));
	}
}

// Generate scene with Spheres of various materials
std::vector<Sphere> RandomScene(const int sphere_count) {
	const float golden_ratio = 1.6180340f;
	auto spheres = std::vector<Sphere>(std::max(sphere_count, 4));

	// Sphere to make up the ground
	spheres[0] = Sphere(Vec3(0, -1000, 0), 1000, Lambertian(Vec3(0.5f, 0.5f, 0.5f)));
	// Three default spheres
	spheres[1] = Sphere(Vec3(0.f, 1.f, 0.f), 1.f, Dielectric(1.5f));
	spheres[2] = Sphere(Vec3(-4, 1, 0), 1.f, Lambertian(Vec3(0.4f, 0.2f, 0.1f)));
	spheres[3] = Sphere(Vec3(4, 1, 0), 1.f, Metal(Vec3(0.7f, 0.6f, 0.5f), 0.0f));

	// Create small spheres
	int i = 4;
	for (; i < sphere_count; i++)
	{
		const float angle = golden_ratio * i;
		const float dist = angle * 0.2f + 1;
		const Vec3 center(cos(angle) * dist, 0.2f, sin(angle) * dist);
		
		spheres[i] = Sphere(center, 0.2f, GenerateRandomMaterial());
	}

	return spheres;
}

///////////////////////////////////////////////////////////////
////////////////////////// MAIN LOOP //////////////////////////
///////////////////////////////////////////////////////////////

int main()
{
	///////////////////// Setup glfw and glew /////////////////////
	glfwInit();
	GLFWwindow* window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "Path Tracer", NULL, NULL);
	glfwMakeContextCurrent(window);

	glewInit();
	glViewport(0, 0, WINDOW_WIDTH, WINDOW_HEIGHT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, WINDOW_WIDTH, 0, WINDOW_HEIGHT, -1, 1);

	Camera cam;
	SetupCamera(cam, INITIAL_CAMERA_POSITION);


	////////////////////// Generate spheres ///////////////////////
	std::vector<Sphere> spheres = RandomScene(SPHERE_COUNT);


	//////////////////////// Setup shaders ////////////////////////
	GLuint comp_shader = glCreateShader(GL_COMPUTE_SHADER);
	glShaderSource(comp_shader, 1, &CALC_SRC, NULL);
	glCompileShader(comp_shader);

	static char LOG_BUFFER[1024];
	glGetShaderInfoLog(comp_shader, 1024, NULL, LOG_BUFFER);

	printf(LOG_BUFFER);

	//////////////////////// Setup program ////////////////////////
	GLuint program = glCreateProgram();
	glAttachShader(program, comp_shader);
	glLinkProgram(program);

	GLint u_time = glGetUniformLocation(program, "time");
	GLint u_sphereCount = glGetUniformLocation(program, "sphereCount");
	GLint u_aa_iterations = glGetUniformLocation(program, "aaIterations");
	GLint u_window_width = glGetUniformLocation(program, "windowWidth");
	GLint u_window_height = glGetUniformLocation(program, "windowHeight");

	glUseProgram(program);


	//////////////////////// Setup output texture //////////////////////// 
	GLuint tex_output;
	glGenTextures(1, &tex_output);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, tex_output);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, WINDOW_WIDTH, WINDOW_HEIGHT, 0, GL_RGBA, GL_FLOAT,
		NULL);
	glBindImageTexture(0, tex_output, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);

	
	glUniform1i(u_window_width, WINDOW_WIDTH);
	glUniform1i(u_window_height, WINDOW_HEIGHT);


	//////////////////////// Send camera data ////////////////////////
	SendCameraDataToComputeShader(cam);


	///////////// Fill a buffer with window coordinates /////////////
	std::vector<WindowCoords> points(PIXEL_COUNT);

	for (int j = WINDOW_HEIGHT - 1; j >= 0; j--)
	{
		for (int i = 0; i < WINDOW_WIDTH; i++)
		{
			int idx = j * WINDOW_WIDTH + i;
			points[idx].x = i;
			points[idx].y = j;
		}
	}
	GLuint pos_ssbo;
	glGenBuffers(1, &pos_ssbo);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, pos_ssbo);
	glBufferData(GL_SHADER_STORAGE_BUFFER, PIXEL_COUNT * sizeof(WindowCoords), points.data(), GL_STATIC_DRAW);


	

	////////////// Gather sphere data for compute shader ////////////////
	std::vector<SphereData> sphere_datas(SPHERE_COUNT);

	for (int i = 0; i < SPHERE_COUNT; i++)
	{
		SphereData& sphere_data = sphere_datas[i];

		const Sphere& sphere = spheres[i];
		sphere_data.center = sphere.center;
		sphere_data.center.w = sphere.radius;


		const Material& material = sphere.material;
		switch (material.matType)
		{
		case MaterialType::Dielectric:
			sphere_data.reflection_index = material.dielectric.reflection_index;
			break;
		case MaterialType::Lambertian:
			sphere_data.albedo = material.lambertian.albedo;
			break;
		case MaterialType::Metal:
			sphere_data.albedo = material.metal.albedo;
			sphere_data.fuzz = material.metal.fuzz;
			break;
		}
		sphere_data.mat_type = static_cast<int>(material.matType);
	}
	GLuint sphere_ssbo;
	glGenBuffers(1, &sphere_ssbo);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 6, sphere_ssbo);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(SphereData) * SPHERE_COUNT, sphere_datas.data(), GL_STATIC_DRAW);


	glUniform1i(u_sphereCount, SPHERE_COUNT);
	glUniform1i(u_aa_iterations, AA_ITERATIONS);
	
	
	auto frame_start = std::chrono::steady_clock::now();

	//////////////////////// Render Loop ////////////////////////
	while (!glfwWindowShouldClose(window))
	{
		glClearColor(0.1f, 0.1f, 0.8f, 1.f);
		glClear(GL_COLOR_BUFFER_BIT);
		
		const double time = glfwGetTime();
		glUseProgram(program);
		glUniform1f(u_time, static_cast<float>(time));

		glDispatchCompute(PIXEL_COUNT / WORK_GROUP_SIZE, 1, 1);
		glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);


		glBindTexture(GL_TEXTURE_2D, tex_output);
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

		glfwSwapBuffers(window);
		const auto frame_end = std::chrono::steady_clock::now();
		const auto frame_delta = std::chrono::duration_cast<std::chrono::nanoseconds>(frame_end - frame_start).count();
		frame_start = std::chrono::steady_clock::now();
		const float delta_time = static_cast<float>(frame_delta / 1000000000.0);
		
		glfwPollEvents();
		UpdateCamera(window, cam, delta_time);

		const float fps = 1.f / delta_time;
		printf("Fps: %f\n", fps);
	}
}