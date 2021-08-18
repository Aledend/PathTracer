#include <stdio.h>
#include <stdlib.h>
#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <GL/GL.h>
#include <limits>
#include "Math/Vec3.h"
#include "Math/Ray.h"
#include "Physics/Collision.h"
#include "Geometry/Sphere.h"
#include "Geometry/Hittable.h"
#include "Framework/HittableList.h"
#include "Camera.h"


const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 600;
const int PIXEL_COUNT = WINDOW_HEIGHT * WINDOW_WIDTH;
const int WORK_GROUP_SIZE = 100;

struct computestruct {
	float u, v, parameter, t;
};
struct colorstruct {
	float r, g, b, a;
};
struct spherestruct {
	Vec3 center;
	float radius;
};


const char* CALC_SRC = R"(#version 430 compatibility
//#extension GL_ARB_compute_shader: enable
//#extension GL_ARB_shader_storage_buffer_object: enable

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
	vec4 sphere_input[];
};

layout (std430, binding = 7) buffer SphereCount
{
	int sphereCount;
};

layout (std430, binding = 8) buffer antiAliasingIterations
{
	int aaIterations;
};

layout (std430, binding = 9) buffer antiAliasingRandoms
{
	float aaRandoms[];
};

layout(local_size_x = 100, local_size_y = 1, local_size_z = 1) in;

const vec3 lowerLeftCorner = vec3(-2.0, -1.5, -1.0);
const vec3 horizontal = vec3(4.0, 0.0, 0.0);
const vec3 vertical = vec3(0.0, 3.0, 0.0);
const vec3 origin = vec3(0.0, 0.0, 0.0);


void GetRay(in float u, in float v, out vec3 A, out vec3 B)
{
	A = origin;
	B = lowerLeftCorner + u*horizontal + v*vertical - origin;
}

vec3 PointAtParameter(in float t, in vec3 A, in vec3 B)
{
	return A + B * t;
}

bool Hit(in vec3 sphereCenter, in float sphereRadius, in vec3 rayA, in vec3 rayB, in float tMin, in float tMax, out float t, out vec3 p, out vec3 normal)
{
	vec3 oc = rayA - sphereCenter;
	float a = dot(rayB, rayB);
	float b = dot(oc, rayB);
	float c = dot(oc, oc) - sphereRadius * sphereRadius;
	float discriminant = b * b - a * c;
	if(discriminant > 0)
	{
		float temp = (-b - sqrt(discriminant)) / a;
		if(temp < tMax && temp > tMin)
		{
			t = temp;
			p = PointAtParameter(t, rayA, rayB);
			normal = (p - sphereCenter) / sphereRadius;
			return true;
		}
		temp = (-b + sqrt(discriminant)) / a;
		if(temp < tMax && temp > tMin)
		{
			t = temp;
			p = PointAtParameter(t, rayA, rayB);
			normal = (p - sphereCenter) / a;
			return true;
		}	
	}
	return false;
}

vec3 Color(in vec3 rayA, in vec3 rayB)
{
	bool hitAnything = false;
	float closestSoFar = 1000000;
	float t = 0.0;
	vec3 p;
	vec3 normal;


	for(int i = 0; i < sphereCount; i++)
	{

		vec3 sphereCenter = vec3(sphere_input[i].x, sphere_input[i].y, sphere_input[i].z);
		float sphereRadius = sphere_input[i].w;

		float tempt;
		vec3 tempp;
		vec3 tempnormal;

		if(Hit(sphereCenter, sphereRadius, rayA, rayB, 0.0, 1000000, tempt, tempp, tempnormal))
		{
			hitAnything = true;
			closestSoFar = tempt;
			t = tempt;
			p = tempp;
			normal = tempnormal;
		}	
	}

	if(hitAnything)
	{
		return 0.5 * vec3(normal.x + 1.0, normal.y + 1.0, normal.z + 1.0);
	}
	else
	{
		vec3 unitDirection = normalize(rayB);
		float somet = 0.5 * (unitDirection.y + 1.0);
		return (1.0 - somet) * vec3(1.0, 1.0, 1.0) + somet * vec3(0.5, 0.7, 1.0);
	}
}

void main()
{
	uint id = gl_GlobalInvocationID.x;
	float posx = point_input[id].x;
	float posy = point_input[id].y;
	float parameter = point_input[id].z;

	//color_output[id] = vec4(posx / 800.0, posy / 600.0, 0.2, 1);
	//return;
	vec3 A;
	vec3 B;
	

	vec3 color = vec3(0, 0, 0);
	for(int i = 0; i < aaIterations; i++)
	{
		uint idx = id * aaIterations * 2 + i;
		float u = (posx + aaRandoms[idx]) / 800.0;
		float v = (posy + aaRandoms[idx+1]) / 600.0;

		GetRay(u, v, A, B);

		color += Color(A, B);
	}
	color /= aaIterations;
	
	color_output[id] = vec4(color, 1);
})";


Vec3 color(const Ray& r, HittableList* world) {
	HitRecord rec;
	if (world->Hit(r, 0.0f, std::numeric_limits<float>::max(), rec))
	{
		return 0.5f * Vec3(rec.normal.x() + 1.f, rec.normal.y() + 1.f, rec.normal.z() + 1.f);
	}
	else
	{
		Vec3 unitDirection = r.Direction().Normalized();
		float t = 0.5f * (unitDirection.y() + 1.0f);
		return (1.0f - t) * Vec3(1.0f, 1.0f, 1.0f) + t * Vec3(0.5f, 0.7f, 1.0f);
	}
}


void handle_key_event(GLFWwindow* window, int key, int scancode, int action, int modifiers)
{
	if(action != GLFW_PRESS)
		return;

	if(key == GLFW_KEY_ESCAPE)
		glfwDestroyWindow(window);
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
	//GLuint compShader = glCreateShader(GL_COMPUTE_SHADER);
	//glShaderSource(compShader, 1, &CALC_SRC, NULL);
	//glCompileShader(compShader);

	//static char LOG_BUFFER[1024];
	//glGetShaderInfoLog(compShader, 1024, NULL, LOG_BUFFER);

	//printf(LOG_BUFFER);


	//GLuint program = glCreateProgram();
	//glAttachShader(program, compShader);
	//glLinkProgram(program);


	//glUseProgram(program);




	//computestruct* points = new computestruct[PIXEL_COUNT];

	//for (int j = WINDOW_HEIGHT - 1; j >= 0; j--)
	//{
	//	for (int i = 0; i < WINDOW_WIDTH; i++)
	//	{
	//		int idx = j * WINDOW_WIDTH + i;
	//		points[idx].u = i;// + static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	//		points[idx].v = j;// + static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	//		points[idx].parameter = 2.0f;
	//		points[idx].t = 2.f;
	//	}
	//}

	//GLuint posSSbo;
	//glGenBuffers(1, &posSSbo);
	//glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, posSSbo);
	//glBufferData(GL_SHADER_STORAGE_BUFFER, PIXEL_COUNT * sizeof(computestruct), points, GL_DYNAMIC_DRAW);



	//colorstruct* colors = new colorstruct[PIXEL_COUNT];


	//GLuint colSSbo;
	//glGenBuffers(1, &colSSbo);
	//glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, colSSbo);
	//glBufferData(GL_SHADER_STORAGE_BUFFER, PIXEL_COUNT * sizeof(colorstruct), colors, GL_DYNAMIC_DRAW);

	//int sphereCount = 2;
	//spherestruct* spheres = new spherestruct[sphereCount];
	//spheres[0].center = Vec3(0.f, 0.f, -1.f);
	//spheres[0].radius = 0.5f;
	//spheres[1].center = Vec3(-1.5f, 0.f, -1.f);
	//spheres[1].radius = 0.3f;



	//GLuint sphereSSbo;
	//glGenBuffers(1, &sphereSSbo);
	//glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 6, sphereSSbo);
	//glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(spherestruct) * sphereCount, spheres, GL_DYNAMIC_DRAW);

	//
	//GLuint sphereCountSSbo;
	//glGenBuffers(1, &sphereCountSSbo);
	//glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 7, sphereCountSSbo);
	//glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(GLint), &sphereCount, GL_DYNAMIC_DRAW);


	//const int antiAliasingIterations = 100;
	//GLuint aaIterationsSSbo;
	//glGenBuffers(1, &aaIterationsSSbo);
	//glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 8, aaIterationsSSbo);
	//glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(GLint), &antiAliasingIterations, GL_DYNAMIC_DRAW);


	//float* aaRandoms = new float[PIXEL_COUNT * antiAliasingIterations];
	//for (int i = 0; i < PIXEL_COUNT * antiAliasingIterations; i++)
	//{
	//	aaRandoms[i] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	//}

	//GLuint aaSSbo;
	//glGenBuffers(1, &aaSSbo);
	//glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 9, aaSSbo);
	//glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(GLfloat) * PIXEL_COUNT * antiAliasingIterations, aaRandoms, GL_DYNAMIC_DRAW);



	//glDispatchCompute(PIXEL_COUNT / WORK_GROUP_SIZE, 1, 1);
	//glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
	//glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, 0);
	//glBindBuffer(GL_SHADER_STORAGE_BUFFER, colSSbo);


	//colorstruct* ptr = (colorstruct*)glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_WRITE_ONLY);
	//for (int i = 0; i < PIXEL_COUNT; i++)
	//{
	//	colors[i] = ptr[i];
	//}
	//glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);




	//printf("%f, %f, %f", colors[0].r, colors[0].g, colors[0].b);

#pragma endregion // Compute Shader test

	const int nx = WINDOW_WIDTH;
	const int ny = WINDOW_HEIGHT;

	/*Hittable* hittables[2];
	hittables[0] = new Sphere(Vec3(0.f, 0.f, -1.f), 0.5f);
	hittables[1] = new Sphere(Vec3(0.f, -100.5f, -1.f), 100.f);
	HittableList* world = new HittableList(hittables, 2);*/

	//Camera cam;

	while (!glfwWindowShouldClose(window))
	{
		glClearColor(0.1f, 0.1f, 0.8f, 1.f);
		glClear(GL_COLOR_BUFFER_BIT);
		glClear(GL_DEPTH_BUFFER_BIT);

		/*glLoadIdentity();
		gluOrtho2D(0.0, 500.0, 500.0, 0.0);*/



		glBegin(GL_POINTS);
		
		#pragma region Compute Shader
		//for (int j = ny - 1; j >= 0; j--)
		//{
		//	for (int i = 0; i < nx; i++)
		//	{
		//		int idx = j * nx + i;

		//		//printf("%f, %f, %f\n", colors[idx].r, colors[idx].g, colors[idx].b);

		//		glColor3f(colors[idx].r, colors[idx].g, colors[idx].b);
		//		glVertex2i(i, j);
		//	}
		//}
		#pragma endregion // Compute Shader Test


		//for (int j = ny - 1; j >= 0; j--)
		//{
		//	for (int i = 0; i < nx; i++)
		//	{
		//		Vec3 col(0.f, 0.f, 0.f);
		//		for (int s = 0; s < ns; s++)
		//		{
		//			float r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		//			float r2 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

		//			const float u = static_cast<float>(i + r1) / static_cast<float>(nx);
		//			const float v = static_cast<float>(j + r2) / static_cast<float>(ny);
		//			const Ray ray = cam.GetRay(u, v);
		//			const Vec3 p = ray.PointAtParameter(2.0f);
		//			col += color(ray, world);
		//		}
		//		col /= float(ns);

		//		
		//		//glColor3f(col.r(), col.g(), col.b());
		//		//glVertex2i(i, j);
		//	}
		//}


		glEnd();

		glfwSwapBuffers(window);

		glfwPollEvents();
	}

	//delete[] colors;
	//delete points;
	/*delete world;
	delete[] hittables;*/
}