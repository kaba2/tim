# Description: ECMake

# Turn on solution folders.
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Set a default build type if none was specified
# ----------------------------------------------

if (NOT CMAKE_BUILD_TYPE)
	message(STATUS "Setting build type to 'Release' as none was specified.")
	set(CMAKE_BUILD_TYPE Release)
endif()

if (CMAKE_CONFIGURATION_TYPES)
	# Restrict multi-configuration generators to the current build-type.
	set (CMAKE_CONFIGURATION_TYPES ${CMAKE_BUILD_TYPE})
endif()

option (BuildMatlabMex
	"Make libraries usable for Matlab mex (force release-mode C and C++ standard libraries)." 
	ON)

# Form tool-set string to differentiate builds
# --------------------------------------------

# Find out whether the generator is 32-bit or 64-bit.
math(EXPR GENERATOR_BITS "8*${CMAKE_SIZEOF_VOID_P}")

# Find out the compiler-id in lower-case.
# For example: msvc, gnu, clang
string (TOLOWER ${CMAKE_CXX_COMPILER_ID} CompilerId)
string (TOLOWER "${CMAKE_BUILD_TYPE}" LOWER_CMAKE_BUILD_TYPE)
string (TOUPPER "${CMAKE_BUILD_TYPE}" UPPER_CMAKE_BUILD_TYPE)

# We use a tool-set id to separate the outputs of 
# different compilers to different directories.
# The tool-set id consists of a compiler-id,
# the bitness of the generator, and the build-type. 
# For example:
# msvc64-release: Visual Studio, 64 bits, release
# gnu32-debug: GCC, 32 bits, debug
# clang64-debug: Clang, 64 bits, debug
set (ToolSet "${CompilerId}${GENERATOR_BITS}-${LOWER_CMAKE_BUILD_TYPE}")

# Force to use an out-of-source build
# -----------------------------------

if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
	# This is an in-source build. Report an error.
	message (SEND_ERROR 
   		"${CMAKE_PROJECT_NAME} does not allow in-source builds (e.g. 'cmake .'); you should do an "
   		"out-of-source build instead (e.g. 'cmake ..' in 'build_${ToolSet}/' directory). "
   		"This call produced the file 'CMakeCache.txt' and the 'CMakeFiles' directory "
   		"in the ${CMAKE_PROJECT_NAME}'s source directory. You must remove them for the out-of-source "
   		" build to work; otherwise CMake attempts an in-source build again."
	)

   return()
endif()

# Define output directories
# -------------------------

set (ProjectDirectory "${CMAKE_SOURCE_DIR}")
set (ProjectIncludeDirectory "${ProjectDirectory}")
set (ProjectLibraryDirectory "${ProjectDirectory}/lib/${ToolSet}")
set (ProjectExecutableDirectory "${ProjectDirectory}/bin/${ToolSet}")
set (ProjectMatlabDirectory "${ProjectExecutableDirectory}/matlab")

set (${CMAKE_PROJECT_NAME}Directory ${ProjectDirectory})
set (${CMAKE_PROJECT_NAME}IncludeDirectory ${ProjectIncludeDirectory})
set (${CMAKE_PROJECT_NAME}LibraryDirectory ${ProjectLibraryDirectory})
set (${CMAKE_PROJECT_NAME}ExecutableDirectory ${ProjectExecutableDirectory})
set (${CMAKE_PROJECT_NAME}MatlabDirectory ${ProjectMatlabDirectory})

include_directories (${ProjectIncludeDirectory})

# Set output directories
# ----------------------

# The directory to place the static libraries (e.g. lib/msvc64-release).
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${ProjectLibraryDirectory})

# The directory to place the shared libraries (e.g. lib/msvc64-release).
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${ProjectLibraryDirectory})

# The directory to place the built executables (e.g. bin/msvc64-release).
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${ProjectExecutableDirectory})

# This is for the multi-configuration build-scripts
# (such as Visual Studio and XCode).
foreach (OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES})
    string (TOUPPER ${OUTPUTCONFIG} UPPER_OUTPUTCONFIG)
    string (TOLOWER ${OUTPUTCONFIG} LOWER_OUTPUTCONFIG)

	# The library output directory is of the form "lib/msvc64-release".
    set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${UPPER_OUTPUTCONFIG} 
    	"${ProjectLibraryDirectory}")
    set (CMAKE_LIBRARY_OUTPUT_DIRECTORY_${UPPER_OUTPUTCONFIG} 
    	"${ProjectLibraryDirectory}")

	# The executable output directory is of the form "bin/msvc64-release".
    set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_${UPPER_OUTPUTCONFIG} 
    	"${ProjectExecutableDirectory}")
endforeach()

# Set some options
# ----------------

if (MSVC)
	# Do not add ZERO_CHECK project into the Visual Studio solution.
	# From VS2008 to VS2013, the ZERO_CHECK project is always out
	# of date, which causes the Visual Studio to ask at every build
	# whether ZERO_CHECK should be built, although nothing was changed.
	# The purpose of the ZERO_CHECK project is to check whether there
	# are changes to the CMake files themselves, and to regenerate the
	# project files if so. But even then the projects are regenerated
	# during the build, and they need to be reloaded, and that is not 
	# very smooth. It is better to suppress this feature and to 
	# regenerate the project files manually whenever the CMake files
	# are changed.
	set(CMAKE_SUPPRESS_REGENERATION TRUE)
endif()

# Helper macros
# -------------

# Checks if the given paths exist.
macro(EcCheckPathExists Name PathSet)
	foreach(Path ${PathSet})
		if(EXISTS ${Path})
			message (STATUS "${Name}: ${Path}")
		else()
			set (Tried "")
			if (NOT ("${Path}" STREQUAL ""))
				set (Tried " (tried ${Path})")
			endif()
			message (SEND_ERROR "Cannot find ${Name}${Tried}. Either install ${Name}, or correct the path in Pastel's root CMakeLists.txt.")
			return()
		endif()
	endforeach()
endmacro()

# Copies a file aside executables.
macro (EcCopyAsideExecutables FilePath)
	if (CMAKE_CONFIGURATION_TYPES)
		# This is a multi-configuration generator,
		# such as Visual Studio or XCode.
		foreach (OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES})
	    	# Copy the file to where the executables are,
	    	# for each configuration.
			file (COPY "${FilePath}" 
				DESTINATION "${ProjectExecutableDirectory}")
		endforeach()
	else()
		# This is a single-configuration generator,
		# such as Unix Makefiles.

		# Copy the file to where the executables are.
		file (COPY "${FilePath}" 
			DESTINATION "${ProjectExecutableDirectory}")
	endif()

	# Copy the file to where the Matlab interface is.
	file (COPY "${FilePath}" 
		DESTINATION "${ProjectExecutableDirectory}/matlab")
endmacro()

# Creates source-groups for files based on the physical directory tree.
macro (EcCreateSourceGroups SourceSet)
foreach (FilePath ${SourceSet})
	# message (STATUS ${FilePath})

	# Get the path to the source file, relative to the current directory.
	file (RELATIVE_PATH FileRelativePath ${CMAKE_CURRENT_LIST_DIR} ${FilePath})

	# Append / to the beginning, so that the regex-replacement
	# works also in the current directory.
	set (FileRelativePath "/${FileRelativePath}")

	# Get the directory-part of the path.
	# I could not find a way for specifying a non-capturing group, 
	# so I opted to append the / to the beginning, and then do
	# the following.
	string (REGEX REPLACE "(.*/)[^/]*$" "\\1" DirectoryRelativePath ${FileRelativePath})

	# Replace / with \.
	string (REPLACE "/" "\\" SourceGroupName ${DirectoryRelativePath})

	# message (STATUS ${FileRelativePath})
	# message (STATUS ${DirectoryRelativePath})
	# message (STATUS ${SourceGroupName})

	# Create a source group.
	source_group(${SourceGroupName} FILES ${FilePath})
endforeach()
endmacro()

# Adds a library, or an executable, and creates source-groups based on 
# the physical directory tree.
macro (EcAddLibrary Type LibraryName SourceGlobSet)
	file (GLOB_RECURSE SourceSet ${SourceGlobSet})

	EcCreateSourceGroups("${SourceSet}")

	#message (STATUS "${LibraryName} is ${Type}" )

	if ("${Type}" STREQUAL "library")
		add_library(${LibraryName} ${SourceSet})
		add_library(${LibraryName}::${LibraryName} ALIAS ${LibraryName})
	elseif ("${Type}" STREQUAL "executable")
		add_executable (${LibraryName} ${SourceSet})
	else ()
		message (FATAL_ERROR "Unknown library type ${Type}.")
	endif()
endmacro()

# Configures a Pastel Matlab library.
#
# Copies each file 'path/name.ext' in the library 
# into 'bin/matlab/path/name.ext'. When the
# file-name is of the form 'path/name.template.ext', 
# where ext is the file-name extension, the
# CMake macros in the file are substituted,
# and the file is renamed to 'path/name.ext'
# before copying.
#
# SourceGlobSet:
# A set of file-globs which determine the
# files to be included in the library.
macro (EcAddMatlabLibrary SourceGlobSet)
	file (GLOB_RECURSE SourceSet RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${SourceGlobSet})
	foreach(FilePath ${SourceSet})
		set (OutputFilePath ${ProjectMatlabDirectory}/${FilePath})
		set (Options "")

		get_filename_component(FileExtension ${FilePath} EXT)

		if (${FileExtension} MATCHES "template.(.+)$")
			string (REGEX REPLACE "(.+).template.(.+)$" "\\1.\\2" OutputFilePath ${OutputFilePath})
		else()
			set (Options COPYONLY)
		endif()
		
		configure_file(${FilePath} ${OutputFilePath} ${Options})
		#message (STATUS "Configured ${FilePath} to ${OutputFilePath}.")
		endforeach()
endmacro()

# Extracts library name from library filename.
macro (EcLibraryNameFromFilename LibraryName LibraryFilename)
	if (UNIX)
		string (REGEX REPLACE "lib(.*)\\.(a|so|dylib)$" "\\1" ${LibraryName} "${LibraryFilename}")
	elseif(WIN32)
		string (REGEX REPLACE "(.*)\\.lib$" "\\1" ${LibraryName} "${LibraryFilename}")
	endif()
endmacro()
