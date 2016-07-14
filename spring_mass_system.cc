#include <iostream>
#include <unistd.h>
#include <string.h>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <Eigen/Dense>

#include "view/shader.h"
#include "utility/obj_reader.h"
#include "utility/spring_mass_moduler.h"
#include "model/elastic_system.h"
#include "algorithm/implicit_euler_integrator.h"

using namespace std;
using namespace Eigen;
using namespace zyclincoln;

void read_file(OBJReader &obj_reader, SpringMassModuler &moduler, const char* path){
  string obj_path = string(path) + ".obj";
  string constraint_path = string(path) + ".txt";

  obj_reader.set_file_name(obj_path);
  obj_reader.ParseModel();

  moduler.AddPoint(obj_reader.points_position());
  moduler.AddSpring(obj_reader.unique_pairs(), obj_reader.edges_length());
  moduler.AddConstraint(constraint_path);
}

GLFWwindow* view_init(){
  glfwInit();
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);

  GLFWwindow* window = glfwCreateWindow(800, 600, "Framework Test", nullptr, nullptr);

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

void build_elastic_system(SpringMassModuler &moduler, ElasticSystem &system){
  system.set_time_step_ms(0.2);
  system.add_points(moduler.point());
  system.add_springs(moduler.spring());
  system.add_static_points(moduler.index_of_static_points());
}

int run(ElasticSystem &system, GLFWwindow* window){
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

int main(int argc, char* argv[]){
  OBJReader obj_reader;
  SpringMassModuler moduler;
  read_file(obj_reader, moduler, argv[1]);

  GLFWwindow* window = view_init();

  ElasticSystem system(moduler.point().size(), moduler.spring().size());
  build_elastic_system(moduler, system);

  return run(system, window);
}