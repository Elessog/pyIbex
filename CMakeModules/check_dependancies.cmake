#check_dependancy.cmake


#set(PYTHON_MIN_VERSION 3.4)
set(PythonInterp_FIND_VERSION "${PYTHON_MIN_VERSION}")
set(PythonLibs_FIND_VERSION "${PYTHON_MIN_VERSION}")
find_package(PythonInterp REQUIRED)
find_package(PythonLibs REQUIRED)

IF(PYTHONLIBS_FOUND)
  INCLUDE_DIRECTORIES("${PYTHON_INCLUDE_DIRS}")
  find_package(NumPy REQUIRED)
  INCLUDE_DIRECTORIES(${NUMPY_INCLUDE_DIRS})
ELSE()
  MESSAGE(FATAL_ERROR "Unable to find PythonLibs.")
ENDIF()

# To statically ling boost with pyIbex on Windows
if(BUILD_STATIC OR WIN32)
  message( "HERE")
  SET(Boost_USE_STATIC_LIBS     ON)
  add_definitions(-DBOOST_PYTHON_STATIC_LIB)
  find_package(Boost COMPONENTS python)
else()
  SET(Boost_USE_STATIC_LIBS     OFF)
  message(STATUS "looking for boost python-py34")
  find_package(Boost COMPONENTS python-py34)
  IF(NOT Boost_FOUND)
    find_package(Boost COMPONENTS python)
  ENDIF()
endif()


IF(Boost_FOUND)
  INCLUDE_DIRECTORIES("${Boost_INCLUDE_DIRS}")
  SET(Boost_USE_MULTITHREADED    ON)
  SET(Boost_USE_STATIC_RUNTIME     ON)
  SET(LIBS ${LIBS} ${Boost_LIBRARIES})
ELSEIF(NOT Boost_FOUND)
  MESSAGE(FATAL_ERROR "Unable to find Boost.")
ENDIF()

FIND_PACKAGE(IbexLib)
if(IBEX_FOUND)
  INCLUDE_DIRECTORIES(${IBEX_INCLUDE_DIRS})
  SET(LIBS ${LIBS} ${IBEX_LIBRARIES} )
else()
  MESSAGE(FATAL_ERROR "Unable to find IbexLib. Need to set IBEX_ROOT ${IBEX_ROOT}")
endif()
