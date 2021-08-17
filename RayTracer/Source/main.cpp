#include <stdio.h>
#include <stdlib.h>
#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <GL/GL.h>
#include "Math/Vec3.h"


int WINDOW_WIDTH = 800;
int WINDOW_HEIGHT = 600;


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


	while (!glfwWindowShouldClose(window))
	{
		glClearColor(0.1f, 0.1f, 0.8f, 1.f);
		glClear(GL_COLOR_BUFFER_BIT);

		Vec3 newVec;
		//glLoadIdentity();
		//gluOrtho2D(0.0, 500.0, 500.0, 0.0);

		glBegin(GL_POINTS);

		const int nx = WINDOW_WIDTH;
		const int ny = WINDOW_HEIGHT;

		for (int j = ny - 1; j >= 0; j--)
		{
			for (int i = 0; i < nx; i++)
			{
				float r = float(i) / float(nx);
				float g = float(j) / float(ny);
				float b = 0.2;
				/*int ir = int(255.99*r);
				int ig = int(255.99*g);
				int ib = int(255.99*b); */
				
				glColor3f(r, g, b);
				glVertex2i(i, j);
			}
		}

		glEnd();

		glfwSwapBuffers(window);

		glfwPollEvents();
	}
}