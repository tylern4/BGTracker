cmake_minimum_required(VERSION 3.4...3.18)
project(BGE CUDA CXX)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)


# Add include directory and cuda include directories
include_directories(${BGE_SOURCE_DIR}/include ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

# Create a library (statatic by default)
add_library(BGTracker
		${BGE_SOURCE_DIR}/src/BGCuda.cu
		${BGE_SOURCE_DIR}/src/Bend.cc
		${BGE_SOURCE_DIR}/src/Drif.cc
		${BGE_SOURCE_DIR}/src/Quad.cu
		${BGE_SOURCE_DIR}/src/Sext.cc
		${BGE_SOURCE_DIR}/src/util.cc)
install(TARGETS BGTracker DESTINATION lib)

option(BUILD_PYBIND11 "Build with BUILD_ support" true)
if(BUILD_PYBIND11)
	# Include all the parts for building pybind11 from the directory
	add_subdirectory(pybind11)
	pybind11_add_module(pyGBCuda python/wrapper.cc)
	target_link_libraries(pyGBCuda PUBLIC BGTracker)
	target_compile_definitions(pyGBCuda PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
	install(TARGETS pyGBCuda DESTINATION lib)
endif()

# Get all the executables from exe folder, compile
file(GLOB exes ${BGE_SOURCE_DIR}/exe/*.c*)
foreach(executable ${exes})
	get_filename_component(exe_name ${executable} NAME_WE)	
	add_executable(${exe_name} ${executable})
	target_link_libraries(${exe_name} PUBLIC BGTracker)
	install(TARGETS ${exe_name} DESTINATION bin)
endforeach()

