# Description: Main CMake build file for the TIM library
# Documentation: building.txt

cmake_minimum_required (VERSION 2.8)

project (Tim)

# Paths
# -----

# Change the following directories to reflect your own
# build environment. If a library is not used, the path
# is ignored.

# The directory of the Pastel library's source code.
# The includes are of the form 'pastel/sys/array.h'
set (PastelIncludeDirectory ../pastel)
set (PastelLibraryDirectory ../pastel/lib)

# The directory of the Boost library's source code.
# The includes are of the form 'boost/static_assert.hpp'.
set (BoostIncludeDirectory ../boost_1_53_0)

# The directory of the Matlab header files.
# The includes are of the form 'mex.h'.
set (MatlabIncludeDirectory /Applications/MATLAB_R2013a.app/extern/include)

# Note: To succesfully _compile_ TIM's Matlab-libraries, 
# you only need the header files for Matlab. 
# This allows you to try the compilation
# even if you do not have the external library
# binaries. It is only in the _linking_ phase of 
# executables and shared libraries (i.e. tests and 
# examples) where the binaries are needed .

# No need to give a library path for Matlab:
# Mex files are built from within Matlab.

# For TIM, the includes are of the form 'tim/core/signal.h'.
include_directories (.)
include_directories (${PastelIncludeDirectory})
include_directories (${BoostIncludeDirectory})
include_directories (${MatlabIncludeDirectory})

# CMake build options
# -------------------

option (BuildLibraries "Build TIM's main libraries." ON)
option (BuildMatlab "Build TIM's Matlab-libraries." ON)
option (BuildTests "Build TIM's tests." ON)

# Force to use an out-of-source build
# -----------------------------------

if ("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
	message (SEND_ERROR 
   		"TIM does not allow in-source builds (e.g. 'cmake .'); you should do an "
   		"out-of-source build instead (e.g. 'cmake ..' in 'release/' directory). "
   		"This call produced the file 'CMakeCache.txt' and the 'CMakeFiles' directory "
   		"in the TIM's source directory. You must remove them for the out-of-source "
   		" build to work; otherwise CMake attempts an in-source build again."
	)

   return()
endif()

# Global variables for the build-script
# -------------------------------------

# Do-it-yourself clang detection
# Under Mac Os X the clang compiler is aliased to c++
# (or something like that; I don't know why cmake
# finds my compiler as /usr/bin/c++).
if (CMAKE_CXX_COMPILER MATCHES "(.*clang)|(.*c\\+\\+)")
	set (CMAKE_COMPILER_IS_CLANGXX 1)
	set (CMAKE_COMPILER_ID Clang)
endif ()

# The libraries are built static.
set (TimLibraryType STATIC)

# The source files to include in the build.
set (TimSourceGlobSet *.cpp *.h *.hpp)

# Set a default build type if none was specified.
# This only applies to single-configuration tool-sets,
# such as Unix Makefiles.
if (NOT CMAKE_BUILD_TYPE AND 
	NOT CMAKE_CONFIGURATION_TYPES)
	message(STATUS 
		"Setting build type to 'Release' as none was specified.")
	set(CMAKE_BUILD_TYPE Release CACHE STRING 
		"Choose the type of build." FORCE)
	# Set the possible values of build types for cmake-gui.
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS 
 		"Debug" 
  		"Release")
endif()

string (TOLOWER "${CMAKE_BUILD_TYPE}" LOWER_CMAKE_BUILD_TYPE)

# The directory to place the static libraries.
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY 
	${CMAKE_SOURCE_DIR}/lib/${LOWER_CMAKE_BUILD_TYPE})

# The directory to place the shared libraries.
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY 
	${CMAKE_SOURCE_DIR}/lib/${LOWER_CMAKE_BUILD_TYPE})

# The directory to place the built executables.
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY 
	${CMAKE_SOURCE_DIR}/bin/${LOWER_CMAKE_BUILD_TYPE})

# This is for the multi-configuration build-scripts
# (such as Visual Studio or XCode).
foreach (OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES})
    string (TOUPPER ${OUTPUTCONFIG} UPPER_OUTPUTCONFIG)
    string (TOLOWER ${OUTPUTCONFIG} LOWER_OUTPUTCONFIG)
    set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${UPPER_OUTPUTCONFIG} 
    	${CMAKE_SOURCE_DIR}/lib/${LOWER_OUTPUTCONFIG})
    set (CMAKE_LIBRARY_OUTPUT_DIRECTORY_${UPPER_OUTPUTCONFIG} 
    	${CMAKE_SOURCE_DIR}/lib/${LOWER_OUTPUTCONFIG})
    set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_${UPPER_OUTPUTCONFIG} 
    	${CMAKE_SOURCE_DIR}/bin/${LOWER_OUTPUTCONFIG})
endforeach()

# TIM's global definitions
# ---------------------------

# Debug-mode definitions.
set (CMAKE_CXX_FLAGS_DEBUG 
	"${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG -D_DEBUG"
)

# clang and g++
# -------------

if (CMAKE_COMPILER_IS_GNUGXX OR CMAKE_COMPILER_IS_CLANGXX) 
	# Same flags apply to both g++ and clang.

	add_definitions (
		# Enables C++11 compiler support.
		-std=c++11 
		-stdlib=libc++
		# Enables some additional warnings.
		-Wall 
	)

	# Enables C++11 linker support.
	set (CMAKE_EXE_LINKER_FLAGS 
		"${CMAKE_EXE_LINKER_FLAGS} -stdlib=libc++")
	set (CMAKE_SHARED_LINKER_FLAGS 
		"${CMAKE_SHARED_LINKER_FLAGS} -stdlib=libc++")

	# Disable some warnings.
	add_definitions (
		# Pragma warnings caused by OpenMP support not being enabled.
		-Wno-unknown-pragmas
		# Comparison between an unsigned and a signed integer.
		-Wno-sign-compare
		# Conversion between an unsigned and a signed integer.
		-Wno-sign-conversion
		# Unused variables.
		-Wno-unused-variable
		# Unused values.
		-Wno-unused-value
		# Unused but set variable.
		# Commented out; Clang does not recognize this.
		# "-Wno-unused-but-set-variable",
		# Unused functions.
		-Wno-unused-function
		# Breaking strict aliasing rules.
		-Wno-strict-aliasing
		# Compiler warns that it optimizes code based on the 
		# assumption that signed integer overflows do not occur.
		-Wno-strict-overflow
		# Compiler warns 'that >= 0' is always true for an 
		# unsigned integer.
		-Wno-tautological-compare
	)
endif()

# Visual Studio
# -------------

if (MSVC)
	add_definitions (
		# Disable Microsoft's Secure STL
		_ITERATOR_DEBUG_LEVEL=0
	)

	# Disable some warnings.
	add_definitions (
		# "'identifier' : truncation from 'type1' to 'type2'"
		/wd4305
		# "'argument' : conversion from 'type1' to 'type2', possible loss of data."
		/wd4244
		# "'identifier' : class 'type' needs to have dll-interface to be used by clients of class 'type2'"
		/wd4251 
		# "new behavior: elements of array 'array' will be default initialized"
		/wd4351
		# "'function': was declared deprecated" (referring to STL functions)
		/wd4996
		# "'var' : conversion from 'size_t' to 'type', possible loss of data"
		/wd4267
		# "'expression' : signed/unsigned mismatch"
		/wd4018 
		# "'operation' : conversion from 'type1' to 'type2' of greater size"
		/wd4312
		# "nonstandard extension used : formal parameter 'identifier' was previously defined as a type"
		/wd4224
		# "qualifier applied to function type has no meaning; ignored"
		/wd4180
		# "'type' : forcing value to bool 'true' or 'false' (performance warning)"
		/wd4800 
		# "'operation' : unsafe use of type 'bool' in operation"
   		/wd4804
	)

endif()

# Recurse to sub-projects
# -----------------------

add_subdirectory (tim)

if (BuildTests)
	#add_subdirectory (test)
endif()
