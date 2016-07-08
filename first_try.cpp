#include <iostream>
#include <unistd.h>
#include <math.h>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <Eigen/Dense>

#include "Shader.hpp"

//#define DEBUG_DETAIL

using namespace std;
using namespace Eigen;

MatrixXd mass = MatrixXd::Identity(9, 9);
double time_step = 0.1;
double length = 1;
double min_y0,max_y0;

void Refresh(VectorXd &point, VectorXd &speed);

int main(){
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR,3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR,3);
	glfwWindowHint(GLFW_OPENGL_PROFILE,GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_RESIZABLE,GL_FALSE);

	GLFWwindow* window = glfwCreateWindow(800, 600, "Spring Net First Try", nullptr, nullptr);
	if(window == nullptr){
		cout<<"Failed to create GLFW window"<<endl;
		glfwTerminate();
		return -1;
	}

	glfwMakeContextCurrent(window);

	glewExperimental = GL_TRUE;
	if(glewInit() != GLEW_OK){
		cout<<"Failed to initialize GLEW"<<endl;
		return -1;
	}

	glViewport(0, 0, 800, 600);
	glEnable(GL_DEPTH_TEST);

	VectorXd point;
	point.resize(9,1);
	point << 0,-0.4,0, -0.3,0,0, 0.3,0,0;

	min_y0 = -0.4;
	max_y0 = 0;

	VectorXd speed;
	speed.resize(9,1);
	speed.setZero();

	GLuint VBO;
	glGenBuffers(1, &VBO);
	GLuint VAO;
	glGenVertexArrays(1, &VAO);
	glBindVertexArray(VAO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(double)*point.rows()*point.cols(), point.data(), GL_DYNAMIC_DRAW);

	glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 3*sizeof(GLdouble), (GLvoid*)0);
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

	Shader PointShader("./PointShader.vertex","./PointShader.fragment");
	PointShader.Use();

	while(!glfwWindowShouldClose(window)){
		glfwPollEvents();
		glClearColor(0.2f,0.3f,0.4f,1.0f);
		glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

		Refresh(point, speed);
		
		#ifdef DEBUG_DETAIL
		getchar();
		#endif

		usleep(10000);

		glBindVertexArray(VAO);
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glBufferData(GL_ARRAY_BUFFER, sizeof(double)*point.rows()*point.cols(), point.data(), GL_DYNAMIC_DRAW);
		glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 3*sizeof(GLdouble), (GLvoid*)0);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindVertexArray(0);

		glBindVertexArray(VAO);
		glDrawArrays(GL_TRIANGLES, 0, 6);
		glBindVertexArray(0);

		glfwSwapBuffers(window);
	}

	glfwTerminate();
	return 0;
}

void Refresh(VectorXd &point, VectorXd &speed){
	MatrixXd A = MatrixXd::Zero(9, 9), stiff = MatrixXd::Zero(9, 9);
	VectorXd b = VectorXd::Zero(9, 1), delta_v = VectorXd::Zero(9, 1), offset_p = VectorXd::Zero(9, 1);

	for(int i=0; i<3; i++){
		int index0 = i, index1 = (i+1)%3;
		Vector3d p0 = point.block(index0*3, 0, 3, 1);
		Vector3d p1 = point.block(index1*3, 0, 3, 1);
		
		Matrix3d sub_stiff = (p1 - p0)*(p1 - p0).transpose()/pow(length, 3) + Matrix3d::Identity()*(1 - 1/length);
		
		stiff.block(index0*3, index0*3, 3, 3) += sub_stiff;
		stiff.block(index0*3, index1*3, 3, 3) += -sub_stiff;
		stiff.block(index1*3, index0*3, 3, 3) += -sub_stiff;
		stiff.block(index1*3, index1*3, 3, 3) += sub_stiff;
		
		offset_p.block(index0*3, 0, 3, 1) += ((p1 - p0).norm() - length)*(p1 - p0).normalized();
		offset_p.block(index1*3, 0, 3, 1) += ((p0 - p1).norm() - length)*(p0 - p1).normalized();
		
		#ifdef DEBUG_DETAIL
		cout<<"===vector p===\n"<<p0<<endl<<"===\n"<<p1<<endl;
		cout<<"===substiff===\n"<<sub_stiff<<endl;
		cout<<"===stiff===\n"<<stiff<<endl;
		#endif

	}

	for(int i=3; i<9; i++){
		stiff.col(i).setZero();
		stiff.row(i).setZero();
	}

	A = mass + pow(time_step,2)*stiff;
	b = time_step*(offset_p-stiff*time_step*speed);

        delta_v.head(3) = A.topLeftCorner(3, 3).colPivHouseholderQr().solve(b.head(3));
	
	speed = speed + delta_v;
	point = point + time_step*speed;

	#ifdef DEBUG_DETAIL
	cout<<"===stiff===\n"<<stiff<<endl;
	cout<<"===speed===\n"<<speed<<endl;
	cout<<"===point===\n"<<point<<endl;
	cout<<"===spring length===\n"<<(point.block(0, 0, 3, 1) - point.block(3, 0, 3, 1)).norm()<<endl;

	#endif

	if(point(1,0)<min_y0){
		min_y0=point(1,0);
		cout<<"new min y0: "<<endl<<point(1,0)<<endl;
	}
	if(point(1,0)>max_y0){
		max_y0=point(1,0);
		cout<<"new max y0: "<<endl<<point(1,0)<<endl;
	}
	cout<<"===spring length===\n"<<(point.block(0, 0, 3, 1) - point.block(3, 0, 3, 1)).norm()<<endl;
}
