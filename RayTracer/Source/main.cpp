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


const int WINDOW_WIDTH = 50;
const int WINDOW_HEIGHT = 10;
const int PIXEL_COUNT = WINDOW_HEIGHT * WINDOW_WIDTH;
const int WORK_GROUP_SIZE = 10;

struct computestruct {
	float u, v, parameter, t;
};
struct colorstruct {
	float r, g, b, a;
};


const char* CALC_SRC = R"(#version 430 compatibility
//#extension GL_ARB_compute_shader: enable
//#extension GL_ARB_shader_storage_buffer_object: enable

//layout(std430, binding = 4) buffer Dat
//{
//	vec4 input_data[];
//};

layout (std430, binding = 5) buffer Col
{
	vec4 output_data[];
};

layout(local_size_x = 100, local_size_y = 1, local_size_z = 1) in;

const vec3 lowerLeftCorner = vec3(-2.0, -1.5, -1.0);
const vec3 horizontal = vec3(4.0, 0.0, 0.0);
const vec3 vertical = vec3(0.0, 3.0, 0.0);
const vec3 origin = vec3(0.0, 0.0, 0.0);
void main()
{
	uint id = gl_GlobalInvocationID.x;
	//float posx = input_data[id].x;
	//float posy = input_data[id].y;
	//float parameter = input_data[id].z;
	//vec3 ray = origin + (lowerLeftCorner + posx * horizontal + posy * vertical - origin) * parameter;
	output_data[0] = vec4(1.0, 1.0, 1.0, 1.0);
})";


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


	GLuint compShader = glCreateShader(GL_COMPUTE_SHADER);
	glShaderSource(compShader, 1, &CALC_SRC, NULL);
	glCompileShader(compShader);

	static char LOG_BUFFER[1024];
	glGetShaderInfoLog(compShader, 1024, NULL, LOG_BUFFER);

	GLuint program = glCreateProgram();
	glAttachShader(program, compShader);
	glLinkProgram(program);


	glUseProgram(program);


	//GLuint buffer;
	//glGenBuffers(1, &buffer);
	//glBindBuffer(GL_ARRAY_BUFFER, buffer);

	//posstruct posData[WINDOW_HEIGHT * WINDOW_WIDTH];

	

	//glBufferData(GL_ARRAY_BUFFER, sizeof(colorData), colorData, GL_STATIC_DRAW);

	//GLuint vao;
	//glGenVertexArrays(1, &vao);
	//glBindVertexArray(vao);

	//glEnableVertexAttribArray(0);
	//glVertexAttribPointer(0, 3, GL_FLOAT, false, 0, 0);






	computestruct* points = new computestruct[PIXEL_COUNT];

	for (int j = WINDOW_HEIGHT - 1; j >= 0; j--)
	{
		for (int i = 0; i < WINDOW_WIDTH; i++)
		{
			points[i].u = i + static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
			points[i].v = j + static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
			points[i].parameter = 2.0f;
			points[i].t = 2.f;
		}
	}

	GLuint posSSbo;
	glGenBuffers(1, &posSSbo);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, posSSbo);
	glBufferData(GL_SHADER_STORAGE_BUFFER, PIXEL_COUNT * sizeof(computestruct), points, GL_DYNAMIC_DRAW);

	

	colorstruct* colors = new colorstruct[PIXEL_COUNT];

	for (int i = 0; i < PIXEL_COUNT; i++)
	{
		colors[i] = {0.f, 0.f, 0.f, 0.f};
	}

	GLuint colSSbo;
	glGenBuffers(1, &colSSbo);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, colSSbo);
	glBufferData(GL_SHADER_STORAGE_BUFFER, PIXEL_COUNT * sizeof(colorstruct), colors, GL_DYNAMIC_DRAW);




	glDispatchCompute(PIXEL_COUNT / WORK_GROUP_SIZE, 1, 1);
	glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, 0);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, colSSbo);


	colorstruct* ptr = (colorstruct*)glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_WRITE_ONLY);
	for (int i = 0; i < 1; i++)
	{
		colors[i] = ptr[i];
	}
	glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);



	printf("%f", colors[0].r);

	printf(LOG_BUFFER);
	return 1;

	const int nx = WINDOW_WIDTH;
	const int ny = WINDOW_HEIGHT;
	const int ns = 10;

	Hittable* hittables[2];
	hittables[0] = new Sphere(Vec3(0.f, 0.f, -1.f), 0.5f);
	hittables[1] = new Sphere(Vec3(0.f, -100.5f, -1.f), 100.f);
	HittableList* world = new HittableList(hittables, 2);

	Camera cam;

	while (!glfwWindowShouldClose(window))
	{
		glClearColor(0.1f, 0.1f, 0.8f, 1.f);
		glClear(GL_COLOR_BUFFER_BIT);
		glClear(GL_DEPTH_BUFFER_BIT);

		/*glLoadIdentity();
		gluOrtho2D(0.0, 500.0, 500.0, 0.0);*/



		glBegin(GL_POINTS);
		
		for (int i = 0; i < PIXEL_COUNT; i++)
		{
			glColor3f(colors[i].r, colors[i].g, colors[i].b);
			float x = i & WINDOW_WIDTH;
			float y = i / WINDOW_WIDTH;
			glVertex2i(x, y);
		}
		


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

		printf("lol");

		glEnd();

		glfwSwapBuffers(window);

		glfwPollEvents();
	}

	delete[] colors;
	//delete points;
	delete world;
	delete[] hittables;
}