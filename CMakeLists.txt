cmake_minimum_required(VERSION 2.8.3)

project(pyIbex)
set(PACKAGE_VERSION "1.1.4")
SET(PYTHON_INSTALL /usr/local/lib/python3.4/dist-packages)
SET(PYTHON_EXECUTABLE /usr/bin/python3.4)
# build type
IF(NOT CMAKE_BUILD_TYPE)
	#SET(CMAKE_BUILD_TYPE "Release")
	SET(CMAKE_BUILD_TYPE "DEBUG")
	#SET(CMAKE_BUILD_TYPE "RELWITHDEBINFO")
	#SET(CMAKE_BUILD_TYPE "MINSIZEREL")
ENDIF()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

IF (False)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
ENDIF()

# option declaration
option(BUILD_TESTS "Build tests" ON)
option(BUILD_EXAMPLE "Build example" OFF)
option(BUILD_STATIC "make standalone lib" OFF)
set(PYTHON_VERSION "3.4" CACHE STRING "target python version")

message(STATUS "Try to use python : ${PYTHON_VERSION}")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/)

SET(PYTHON_EXECUTABLE /usr/bin/python3.4)
SET(PYTHON_INSTALL /usr/local/lib/python3.4/dist-packages)
if(NOT WIN32)
  #SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC"  )
	set(IbexLib_USE_STATIC 1)
	#SET(CMAKE_EXE_LINKER_FLAGS "-static")
endif()
SET(PYTHON_EXECUTABLE /usr/bin/python3.4)
SET(PYTHON_INSTALL /usr/local/lib/python3.4/dist-packages)
SET(PYTHON_LIBRARY /usr/lib/x86_64-linux-gnu/libpython3.4m.so)
SET(PYTHON_INCLUDE_DIR /usr/include/python3.4)
include(${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/check_dependancies.cmake)
# Find python install path (from morse install script)
SET(PYTHON_EXECUTABLE /usr/bin/python3.4)
SET(PYTHON_INSTALL /usr/local/lib/python3.4/dist-packages)
EXECUTE_PROCESS(COMMAND
					${PYTHON_EXECUTABLE} -c "import distutils.sysconfig, sys; sys.stdout.write(distutils.sysconfig.get_python_lib(1,0,\"\"))"
					OUTPUT_VARIABLE PYTHON_INSTALL
					ERROR_VARIABLE PYTHON_STDERR
					RESULT_VARIABLE PYTHON_ERR
				)
SET(PYTHON_INSTALL /usr/local/lib/python3.4/dist-packages)
if(NOT PYTHON_INSTDIR)
	SET(PYTHON_INSTDIR ${PYTHON_INSTALL}) #${CMAKE_INSTALL_PREFIX}/
endif(NOT PYTHON_INSTDIR)

if(WIN32)
    STRING(REGEX REPLACE "\\\\" "/" PYTHON_INSTDIR ${PYTHON_INSTDIR})
endif(WIN32)
MESSAGE(STATUS "will install python files in ${PYTHON_INSTDIR}")

if(WIN32)
	set(CompilerFlags
        CMAKE_CXX_FLAGS
        CMAKE_CXX_FLAGS_DEBUG
        CMAKE_CXX_FLAGS_RELEASE
        CMAKE_C_FLAGS
        CMAKE_C_FLAGS_DEBUG
        CMAKE_C_FLAGS_RELEASE
        )
foreach(CompilerFlag ${CompilerFlags})
  string(REPLACE "/MD" "/MT" ${CompilerFlag} "${${CompilerFlag}}")
endforeach()
endif()

add_subdirectory(src)

if (BUILD_TESTS)
	ENABLE_TESTING()
	include( CTest )
	add_subdirectory(tests)
endif(BUILD_TESTS)

if (BUILD_EXAMPLE)
	add_subdirectory(examples)
endif(BUILD_EXAMPLE)




#########################################################################
### Uninstall ###
#########################################################################
configure_file("${CMAKE_HOME_DIRECTORY}/CMakeModules/cmake_uninstall.cmake.in" "cmake_uninstall.cmake" IMMEDIATE @ONLY)
add_custom_target(uninstall "${CMAKE_COMMAND}" -P "cmake_uninstall.cmake")


#########################################################################
### Package configuration ###
#########################################################################


set(CPACK_PACKAGE_NAME ${CMAKE_PROJECT_NAME})
set(CPACK_PACKAGE_VENDOR "B. Desrochers")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Python wrapper for Ibex Library")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "1")
set(CPACK_PACKAGE_VERSION_PATCH "3")
set(CPACK_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
#set(CPACK_PACKAGE_INSTALL_DIRECTORY "pyIbex")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Benoit Desrochers <ben.ensta@gmail.com>")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libboost-python1.58.0")
SET(CPACK_GENERATOR "DEB")
set (CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
INCLUDE(CPack)


message(STATUS "Using python lib ${PYTHON_LIBRARY}")
