#include <iostream>
#include <unistd.h>
#include <string.h>
#include <sstream>
#include <sys/time.h>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <Eigen/Dense>

#include "view/shader.h"
#include "utility/vtk_reader.h"
#include "utility/tetrahedron_moduler.h"
#include "model/fem_system.h"
#include "algorithm/implicit_euler_integrator.h"
#include "utility/vtk_generator.h"

using namespace std;
using namespace Eigen;
using namespace zyclincoln;

void read_file(VTKReader &vtk_reader, TetrahedronModuler &moduler, const char* path){
  string vtk_path = string(path) + ".vtk";
  string constraint_path = string(path) + ".txt";

  vtk_reader.set_file_name(vtk_path);
  vtk_reader.ParseModel();

  moduler.AddPoint(vtk_reader.points_position());
  moduler.AddTetrahedron(vtk_reader.tetrahedrons_point());
  moduler.AddConstraint(constraint_path);
}

GLFWwindow* view_init(){
  glfwInit();
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);

  GLFWwindow* window = glfwCreateWindow(800, 600, "Linear Tetrahedron Test", nullptr, nullptr);

  if(window == nullptr){
    cerr << "Failed to create GLFW window" <<endl;
    glfwTerminate();
    return nullptr;
  }

  glfwMakeContextCurrent(window);

  glewExperimental = GL_TRUE;
  if(glewInit() != GLEW_OK){
    cout << "Failed to initialize GLEW" << endl;
    return nullptr;
  }

  glViewport(0, 0, 800, 600);
  glEnable(GL_DEPTH_TEST);

  return window;
}

void build_tetrahedron_system(TetrahedronModuler &moduler, LinearFEMSystem &system){
  system.set_time_step_ms(0.3);
  system.add_points(moduler.point());
  system.add_tetrahedrons(moduler.tetrahedron());
  system.add_static_points(moduler.index_of_static_points());
//  system.update_draw_line();
}

int run(LinearFEMSystem &system, GLFWwindow* window){
  ImplicitEulerIntegrator integrator;
  GLuint vbo, vao;
  glGenBuffers(1, &vbo);
  glGenVertexArrays(1, &vao);

  glBindVertexArray(vao);
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glBufferData(GL_ARRAY_BUFFER, sizeof(double)*(system.draw_line().rows()),
    system.draw_line().data(), GL_DYNAMIC_DRAW);
  glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 3*sizeof(GLdouble), (GLvoid*)0);
  glEnableVertexAttribArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);

  Shader draw_shader("PointShader.vertex", "PointShader.fragment");
  draw_shader.Use();

  while(!glfwWindowShouldClose(window)){
    glfwPollEvents();
    glClearColor(0.2f, 0.3f, 0.4f, 0.1f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    system.update_draw_line();
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(double)*(system.draw_line().rows()),
      system.draw_line().data(), GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    
    glBindVertexArray(vao);
    glDrawArrays(GL_LINES, 0, system.draw_line().rows());
    glBindVertexArray(0);
    integrator.next_frame(system);
    glfwSwapBuffers(window);
  }

  glfwTerminate();
  return 0;
}

int main(int argc, char** argv){
  VTKReader vtk_reader;
  TetrahedronModuler moduler;
  read_file(vtk_reader, moduler, argv[1]);

//  GLFWwindow* window = view_init();
  LinearFEMSystem system(moduler.point().size(), moduler.tetrahedron().size());
  build_tetrahedron_system(moduler, system);
//  return run(system, window);

  ImplicitEulerIntegrator integrator;
  integrator.set_reduce();
  // VectorXd old_position = system.position_vector();
  for(unsigned int i = 0; i < 300; i++){ 
    timeval starttime, endtime;
    gettimeofday(&starttime,0);

    cout <<" calculate " << i << " frame : " << endl;
    stringstream ss;
    string output(argv[1]);
    ss << output << "_reduce_result_" << i+1 << ".vtk";
    ss >> output;

    integrator.next_frame(system);
    linear_fem_system_to_vtk(system, output);
    // system.update_position_vector(old_position);
    
    gettimeofday(&endtime,0);
    double timeuse = 1000000*(endtime.tv_sec - starttime.tv_sec) + endtime.tv_usec - starttime.tv_usec;

    cout << "total time: " << timeuse/1000 << "ms" << endl;
  }
  return 0;
}