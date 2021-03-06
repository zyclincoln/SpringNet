project(SpringNet)

cmake_minimum_required(VERSION 2.8.11)
set(CMAKE_CXX_FLAGS "-std=c++11")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModule")


find_package(glfw3 3.2 REQUIRED)

find_package(GLEW)
include_directories(${GLEW_INCLUDE_DIRS})
link_directories(${GLEW_LIBRARY_DIRS})
add_definitions(${GLEW_DEFINITIONS})

if(NOT GLEW_FOUND)
	message(Error " GLEW not found")
endif(NOT GLEW_FOUND)

find_package( OpenMP REQUIRED)  
if(OPENMP_FOUND)  
message("OPENMP FOUND")  
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")  
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")  
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")  
endif()  

find_package(OpenGL REQUIRED)
include_directories(${OPENGL_INCLUDE_DIR})

find_package(Eigen3 REQUIRED)
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})

find_package(Boost COMPONENTS timer REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

SET(LIB_SRC
  algorithm/implicit_euler_integrator.cc
  algorithm/point_potential_energy.cc
  algorithm/spring_potential_energy.cc
  algorithm/tetrahedron_potential_energy.cc
#  model/elastic_system.cc
  model/fem_system.cc
  model/point.cc
  model/spring.cc
  model/matrix_op.cc
  model/tetrahedron.cc
  utility/obj_reader.cc
  utility/spring_mass_moduler.cc
  utility/vtk_reader.cc
  utility/tetrahedron_moduler.cc
)

add_library(spring_mass_system STATIC ${LIB_SRC})

# add_executable(FrameWorkTest framework_test.cc elastic_plane.cc matrix_op.cc obj_reader.cc)
# add_executable(SpringMassSystem spring_mass_system.cc)
# add_executable(Test ./test/vtk_reader_test.cc)
add_executable(LinearFEMSystem tetrahedron_linear_system.cc)
add_executable(CorotateFEMSystem tetrahedron_corotate_system.cc)
add_executable(vol_to_vtk vol_to_vtk.cc)
add_executable(LinearModelReduce linear_model_reduce.cc)

# target_link_libraries(FrameWorkTest
# 	glfw
# 	GLEW
# 	${OPENGL_LIBRARIES}
#   ${Boost_LIBRARIES}
# )

# target_link_libraries(SpringMassSystem
#   glfw
#   GLEW
#   ${OPENGL_LIBRARIES}
#   spring_mass_system
#   ${Boost_LIBRARIES}
# )

# target_link_libraries(Test
#   spring_mass_system
#   ${Boost_LIBRARIES}
# )

target_link_libraries(LinearFEMSystem
  glfw
  GLEW
  ${OPENGL_LIBRARIES}
  spring_mass_system
)

target_link_libraries(CorotateFEMSystem
  glfw
  GLEW
  ${OPENGL_LIBRARIES}
  spring_mass_system
)

target_link_libraries(vol_to_vtk
  spring_mass_system
)

target_link_libraries(LinearModelReduce
  glfw
  GLEW
  ${OPENGL_LIBRARIES}
  spring_mass_system
)
