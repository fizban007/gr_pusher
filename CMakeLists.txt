cmake_minimum_required(VERSION 2.8)
project(GR)

option(build_tests "Build the test suite." ON)

# Set CXX flags
set(CMAKE_CXX_STANDARD 14)
set(CXX_EXTRA_FLAGS "-DBOOST_ALL_NO_LIB=1 -DBOOST_NO_CXX14_CONSTEXPR")
if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  set(CMAKE_CXX_FLAGS_RELEASE "-std=c++14 -pthread -O3 -g -DNDEBUG -Wall -Wextra -fPIC -xHost -no-prec-div ${CXX_EXTRA_FLAGS}")
else()
  # set(CXX_EXTRA_FLAGS "-Wno-misleading-indentation -Wno-unused-parameter")
  set(CXX_EXTRA_FLAGS "-Wno-unused-parameter")
  set(CMAKE_CXX_FLAGS_RELEASE "-std=c++14 -pthread -O3 -g -DNDEBUG -Wall -Wextra -fPIC -march=native ${CXX_EXTRA_FLAGS}")
endif()
set(CMAKE_CXX_FLAGS_DEBUG "-std=c++14 -pthread -g3 -O0 -Wall -Wextra -fPIC ${CXX_EXTRA_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE}")

# Set CMake parameters
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# Find external dependencies
find_package(Eigen3)
if (EIGEN3_FOUND)
  message(${EIGEN3_INCLUDE_DIR})
endif()

# Set include directories
set(CUDALE_INCLUDE_DIR deps/CudaLE/include)
set(CATCH_INCLUDE_DIR deps/catch/single_include)
set(FADBAD_INCLUDE_DIR deps/FADBAD++)

include_directories(${EIGEN3_INCLUDE_DIR} ${CUDALE_INCLUDE_DIR} ${CATCH_INCLUDE_DIR} ${FADBAD_INCLUDE_DIR} include)

# Add src directory to the build process
add_subdirectory(src)

# Include unit tests if build_tests is ON
if (build_tests)
  add_subdirectory(tests)
endif()

# This is to get all header files recognized by QtCreator
file(GLOB_RECURSE headers "include/*.h" "include/*/*.h" "include/*.hpp" "include/*/*.hpp")
add_custom_target(myIncludes SOURCES ${headers})
