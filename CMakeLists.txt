# Description: Main CMake build file for the TIM library
# Documentation: building.txt

cmake_minimum_required (VERSION 2.8)

project (Tim)

# CMake build options
# -------------------

option (BuildLibraries "Build Tim's main libraries." ON)
option (BuildMatlab "Build Tim's Matlab-libraries." ON)
option (BuildTests "Build Tim's tests." ON)

# ECMake
# ------

# Directory for CMake to search for CMake-included() files.
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Initialize ECMake.
include ("ecmake")

# Setup compilers.
include ("SetupCompilers")

# Configure external libraries
# ----------------------------

# Manual override for external library paths.
include ("Paths")

find_package("Ranges_")
include_directories (${RangesIncludeDirectory})

find_package("Pastel_")
include_directories (${PastelIncludeDirectory})

find_package("Boost_")
include_directories (${BoostIncludeDirectory})

find_package("Armadillo_")
include_directories (${ArmadilloIncludeDirectory})

find_package("Tbb_")
include_directories (${TbbIncludeDirectory})

find_package("Blas_")
find_package("Lapack_")

if (BuildMatlab)
	find_package("Matlab_")
	include_directories (${MatlabIncludeDirectory})
endif()

# Copy files
# ----------

if (WIN32)
	# Copy dll-libraries aside executables for
	# Windows to find them, when running an
	# executable.

	set (DllPathSet 
		"${BlasDllPath}"
		"${LapackDllPath}")

	foreach (DllPath ${DllPathSet})
		EcCopyAsideExecutables("${DllPath}")
	endforeach()
endif()

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
	#add_subdirectory (test)
endif()

if (BuildMatlab)
	add_subdirectory (matlab)
endif()

if (BuildExamples)
	#add_subdirectory (example)
endif()

