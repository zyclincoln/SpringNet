#include <iostream>
#include <unistd.h>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <Eigen/Dense>

#include "shader.hpp"

#include "elastic_plane.h"
#include "net_reader.h"

using namespace std;
using namespace Eigen;

int main(){
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);

	GLFWwindow* window = glfwCreateWindow(800, 600, "Framework Test", nullptr, nullptr);

	if(window == nullptr){
		cerr << "Failed to create GLFW window" <<endl;
		glfwTerminate();
		return -1;
	}

	glfwMakeContextCurrent(window);

	glewExperimental = GL_TRUE;
	if(glewInit() != GLEW_OK){
		cout << "Failed to initialize GLEW" << endl;
		return -1;
	}

	glViewport(0, 0, 800, 600);
	glEnable(GL_DEPTH_TEST);

	Net_Reader net_reader("./net.txt");
	net_reader.read_file();
	ElasticPlane elastic_plane(net_reader.length(), net_reader.width());
	elastic_plane.add_points(net_reader.mass(), net_reader.points_position(), net_reader.points_speed());
	elastic_plane.add_static_points(net_reader.index_of_static_points());
	elastic_plane.set_time_step_ms(0.05);
	
	GLuint VBO, VAO, EBO;
	glGenBuffers(1, &VBO);
	glGenBuffers(1, &EBO);
	glGenVertexArrays(1, &VAO);
	
	glBindVertexArray(VAO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);

	elastic_plane.generate_draw_line();
	glBufferData(GL_ARRAY_BUFFER,
		sizeof(double)*(elastic_plane.draw_line().rows()),
		elastic_plane.draw_line().data(),
		GL_DYNAMIC_DRAW);
	// glBufferData(GL_ARRAY_BUFFER, 
	// 	sizeof(double)*(elastic_plane.points_position().rows()), 
	// 	elastic_plane.points_position().data(),
	// 	GL_DYNAMIC_DRAW);

	// glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
 //    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 
 //    	sizeof(unsigned int)*elastic_plane.draw_line_index().rows(), 
 //    	elastic_plane.draw_line_index().data(), 
 //    	GL_DYNAMIC_DRAW);

	glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 3*sizeof(GLdouble), (GLvoid*)0);
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

	Shader point_shader("./PointShader.vertex", "./PointShader.fragment");
	point_shader.Use();

	while(!glfwWindowShouldClose(window)){
		glfwPollEvents();
		glClearColor(0.2f, 0.3f, 0.4f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		elastic_plane.generate_draw_line();
		glBindVertexArray(VAO);
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glBufferData(GL_ARRAY_BUFFER, 
			sizeof(double)*(elastic_plane.draw_line().rows()),
			elastic_plane.draw_line().data(),
			GL_DYNAMIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindVertexArray(0);

		glBindVertexArray(VAO);
		glDrawArrays(GL_LINES, 0, elastic_plane.draw_line_index().rows());
		// glDrawElements(GL_LINES, elastic_plane.draw_line_index().rows(), GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);

		usleep(10000);
		elastic_plane.next_frame();

		//getchar();

		glfwSwapBuffers(window);
	}

	glfwTerminate();
	return 0;
}