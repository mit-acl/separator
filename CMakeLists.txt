cmake_minimum_required(VERSION 2.8.3)
project(separator)

## Compile as C++11, supported in ROS Kinetic and newer
add_compile_options(-std=c++11)

find_package(catkin REQUIRED)
find_package( Eigen3 REQUIRED )
include_directories(${EIGEN3_INCLUDE_DIR})

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#
# Download GLPK
#

set(GLPK_DIR ${CMAKE_CURRENT_BINARY_DIR}/glpk-download)

# Download GLPK
configure_file(
  "${CMAKE_MODULE_PATH}/glpk.cmake.in"
  "${GLPK_DIR}/CMakeLists.txt"
  IMMEDIATE @ONLY)

execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" . WORKING_DIRECTORY ${GLPK_DIR})
execute_process(COMMAND ${CMAKE_COMMAND} --build . WORKING_DIRECTORY ${GLPK_DIR})

# For convenience (link to static library)
set(GLPK_INCLUDE_DIRS ${GLPK_DIR}/install/include)
set(GLPK_LIBS -L${GLPK_DIR}/install/lib;libglpk.a)

###############################################################################

catkin_package(
    INCLUDE_DIRS include ${GLPK_INCLUDE_DIRS}
    LIBRARIES separator
)

include_directories(
	include
    ${GLPK_INCLUDE_DIRS}
	${catkin_INCLUDE_DIRS}
)

add_library(separator SHARED src/separator.cpp)
target_link_libraries(separator ${catkin_LIBRARIES} ${GLPK_LIBS} m)
add_dependencies(separator ${separator_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})



add_executable(test_separator src/test_separator.cpp)
target_link_libraries(test_separator PUBLIC separator)
add_dependencies(test_separator separator)
