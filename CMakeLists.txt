# Decide which version of CMake is needed to process this project.
# Actually not sure which version is needed - older versions probably work, too.
cmake_minimum_required(VERSION 3.0)

# The name of this project. It is possible to add several projects and compile
# them separately.
project(ann2)

# Some compiler specific things need to be done.
# If your using neither GNU nor Visual Studio compiler, you might want/have to
# add a matching case to control specific aspects.
if (${CMAKE_CXX_COMPILER_ID} STREQUAL GNU)
  # The GNU compiler needs to be told to use C++11. Also, warnings are nice.
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")
elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
  # The Visual Studio compiler needs to be identified in some places of the C++
  # code. Define a preprocessor variable arbitrarily called ANN_USE_MSVC.
  add_definitions(-D_USE_MATH_DEFINES)
endif (${CMAKE_CXX_COMPILER_ID} STREQUAL GNU)

# Files belonging to the project are listed here and collected in a variable.
set(SOURCE_FILES
  main.cpp
  floating.h
  matrix.h
  matrix.cpp
  x_factory.h
  x_factory.cpp
  y_factory.h
  y_factory.cpp
)

# Create an executable file with the same name as this project and use the
# collected source files as code. Only one of the files in the list may contain
# a main function.
add_executable(${PROJECT_NAME} ${SOURCE_FILES})
