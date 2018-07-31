cmake_minimum_required(VERSION 2.6)
project(xfoil)

ADD_DEFINITIONS(-std=c++11)

find_package(Qt5Core REQUIRED)
find_package(pybind11 REQUIRED)
find_package(Eigen3 REQUIRED)

EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} -c
"from distutils.sysconfig import get_python_lib; import os; print(os.path.join(get_python_lib()))"
OUTPUT_VARIABLE _PYTHON_INSTALL_DIR OUTPUT_STRIP_TRAILING_WHITESPACE )
string(REPLACE "\\" "/" PYTHON_INSTALL_DIR ${_PYTHON_INSTALL_DIR})


include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${Qt5Core_INCLUDE_DIRS})
INCLUDE_DIRECTORIES( ${pybind11_INCLUDE_DIRS})

pybind11_add_module(xfoil XFoil.cpp xfoil.cpp)
target_link_libraries(xfoil PUBLIC Qt5::Core)

install(TARGETS xfoil DESTINATION ${PYTHON_INSTALL_DIR})