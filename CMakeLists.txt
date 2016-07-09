project(SpringNet)

cmake_minimum_required(VERSION 2.8.11)
set(CMAKE_CXX_FLAGS "-std=c++11")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModule")


# find_package(glfw3 3.2 REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_search_module(GLFW3 REQUIRED glfw3)

find_package(GLEW)
include_directories(${GLEW_INCLUDE_DIRS})
link_directories(${GLEW_LIBRARY_DIRS})
add_definitions(${GLEW_DEFINITIONS})

if(NOT GLEW_FOUND)
	message(Error " GLEW not found")
endif(NOT GLEW_FOUND)

find_package(OpenGL REQUIRED)
include_directories(${OPENGL_INCLUDE_DIR})

find_package(Eigen3 REQUIRED)
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})

find_package(Boost COMPONENTS system thread REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})

add_executable(FirstTry first_try.cpp)
add_executable(FrameWorkTest framework_test.cc elastic_plane.cc matrix_op.cc obj_reader.cc)

target_link_libraries(FirstTry
#	glfw
	GLEW
	${OPENGL_LIBRARIES}
	${GLFW3_STATIC_LIBRARIES}
)

target_link_libraries(FrameWorkTest
#	glfw
	GLEW
	${OPENGL_LIBRARIES}
	${GLFW3_STATIC_LIBRARIES}
)
