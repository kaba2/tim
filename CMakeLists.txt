# Description: Main CMake build file for the TIM library
# Documentation: building.txt

cmake_minimum_required (VERSION 3.16)

if(NOT "${CMAKE_CXX_STANDARD}")
  set(CMAKE_CXX_STANDARD 20)
endif()

# Set up VCPKG
# ------------

set(CMAKE_TOOLCHAIN_FILE ${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake
  CACHE STRING "Vcpkg toolchain file")

if(DEFINED ENV{VCPKG_ROOT} AND NOT DEFINED CMAKE_TOOLCHAIN_FILE)
  set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake"
      CACHE STRING "")
endif()

if(DEFINED ENV{VCPKG_DEFAULT_TRIPLET} AND NOT DEFINED VCPKG_TARGET_TRIPLET)
  set(VCPKG_TARGET_TRIPLET "$ENV{VCPKG_DEFAULT_TRIPLET}" CACHE STRING "")
endif()

project (Tim)

# CMake build options
# -------------------

option (BuildLibraries "Build Tim's main libraries." ON)
option (BuildMatlab "Build Tim's Matlab-libraries." ON)
option (BuildTests "Build Tim's tests." ON)

# ECMake
# ------

set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Initialize ECMake.
include ("ecmake")

# Setup compilers.
include ("SetupCompilers")

message(STATUS 
	"Current compiler is ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")

if (CMAKE_BUILD_TYPE MATCHES Debug)
	add_definitions(
		-DPASTEL_ENABLE_ASSERTS
	)
endif()

# Configure external libraries
# ----------------------------

# Directory for CMake to search for CMake-included() files.
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set (PastelDirectory "${ProjectDirectory}/../pastel")
set (Pastel_DIR "${PastelDirectory}/cmake")

find_package(Pastel CONFIG REQUIRED)
include_directories (${PastelIncludeDirectory})
link_directories (${PastelLibraryDirectory})

find_package(range-v3 CONFIG REQUIRED)
include_directories (range-v3)
add_definitions(-DRANGES_NO_STD_FORWARD_DECLARATIONS)
get_target_property(RangesIncludeDirectory range-v3 INTERFACE_INCLUDE_DIRECTORIES)

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_libraries(${Boost_LIBRARIES})

add_definitions(
	-DBOOST_BIND_GLOBAL_PLACEHOLDERS
)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIRS})
link_libraries(Eigen3::Eigen)

add_definitions (
	# std::result_of was removed in C++20.
	-DEIGEN_HAS_STD_RESULT_OF=0
)

find_package(TBB CONFIG REQUIRED)
include_directories(TBB::tbb)
link_libraries(TBB::tbb)
get_target_property(TbbIncludeDirectory TBB::tbb INTERFACE_INCLUDE_DIRECTORIES)

if (CMAKE_BUILD_TYPE MATCHES "Debug")
	get_target_property(TbbLibraryPath TBB::tbb IMPORTED_LOCATION_DEBUG)
	get_target_property(TbbMallocLibraryPath TBB::tbbmalloc IMPORTED_LOCATION_DEBUG)
else()
	get_target_property(TbbLibraryPath TBB::tbb IMPORTED_LOCATION_RELEASE)
	get_target_property(TbbMallocLibraryPath TBB::tbbmalloc IMPORTED_LOCATION_RELEASE)
endif()

get_filename_component(TbbLibraryDirectory ${TbbLibraryPath} DIRECTORY)
get_filename_component(TbbLibraryName ${TbbLibraryPath} NAME)

get_filename_component(TbbMallocLibraryDirectory ${TbbMallocLibraryPath} DIRECTORY)
get_filename_component(TbbMallocLibraryName ${TbbMallocLibraryPath} NAME)

if (UNIX)
	# FIX: This is a hack to get the library name correct under Linux.
	string (REPLACE lib "" TbbLibraryName ${TbbLibraryName})
	string (REPLACE .a "" TbbLibraryName ${TbbLibraryName})

	string (REPLACE lib "" TbbMallocLibraryName ${TbbMallocLibraryName})
	string (REPLACE .a "" TbbMallocLibraryName ${TbbMallocLibraryName})
endif()

if (BuildMatlab)
	set (MatlabDirectoryWindows "C:/Program Files/Polyspace/R2020b")
	set (MatlabDirectoryMac "/Applications/MATLAB_R2020b.app")
	set (MatlabDirectoryLinux "/usr/local/Polyspace/R2020b")

	if (WIN32)
		set (MatlabDirectory ${MatlabDirectoryWindows})
	elseif(UNIX)
		if (APPLE)
			set (MatlabDirectory ${MatlabDirectoryMac})
		else()
			set (MatlabDirectory ${MatlabDirectoryLinux})
		endif()
	endif()

	set (MatlabIncludeDirectory "${MatlabDirectory}/extern/include")

	EcCheckPathExists("Matlab (include)" "${MatlabIncludeDirectory}")

	include_directories (${MatlabIncludeDirectory})
endif()

message(STATUS "Range-v3: ${RangesIncludeDirectory}")
message(STATUS "Boost: ${Boost_INCLUDE_DIRS} ${Boost_LIBRARIES}")
message(STATUS "Eigen: ${EIGEN3_INCLUDE_DIRS}")
message(STATUS "TBB: ${TbbIncludeDirectory} ${TbbLibraryPath}")
message(STATUS "Tbbmalloc: ${TbbMallocIncludeDirectory} ${TbbMallocLibraryPath}")

# Source file globs
# -----------------

# The source files to include in a C++ project.
set (TimSourceGlobSet *.cpp *.h *.hpp *.txt)

# The source files to include in a Matlab project.
set (TimMatlabSourceGlobSet *.m *.m.cmake)

# The documentation files of the project.
set (TimDocumentationGlobSet *.txt)

# Recurse to sub-projects
# -----------------------

add_subdirectory (tim)

if (BuildTests)
	add_subdirectory (test)
endif()

if (BuildMatlab)
	add_subdirectory (matlab)
endif()

if (BuildExamples)
	#add_subdirectory (example)
endif()

