-- This is a Premake script for producing
-- build files for the TIM library.

-- Available external libraries
-- ============================

-- Whether you have Pastel 1.2 include files and binaries.
gotPastel = true

-- Whether you have Boost 1.45 include files.
gotBoost = true

-- Whether you have SDL 1.2 include files and binaries.
gotSdl = true

-- Whether you have Matlab (say, 2008a or never) include files.
gotMatlab = true

-- Note: To succesfully _compile_ the libraries, 
-- you only need the header files for the external
-- libraries. This allows you to try the compilation
-- even if you did not have the external library
-- binaries. It is only in the _linking_ phase of 
-- executables and shared libraries (i.e. tests and 
-- examples) where the binaries are needed .

-- Main build switches
-- ===================

-- Whether to build the libraries.
buildLibraries = true

-- Whether to build the tests.
buildTests = true

-- Detailed build switches
-- =======================

-- Requirements
basicRequirements = gotPastel and gotBoost
libraryRequirements = basicRequirements and buildLibraries
testRequirements = basicRequirements and buildTests and gotSdl

-- Whether to build the TIM Core library.
buildTimCore = true and libraryRequirements

-- Whether to build the TIM Console interface.
buildTimConsole = true and libraryRequirements

-- Whether to build the TIM Matlab interface.
buildTimMatlab = true and libraryRequirements and gotMatlab

-- Whether to build the test project.
buildTimCoreTest = true and testRequirements

-- Paths
-- =====

-- Change the following directories to reflect your own
-- build environment.

-- The directory of the Pastel library's source code.
-- The includes are of the form 'pastel/sys/array.h'
pastelIncludeDir = "../pastel"
pastelLibraryDir = "../pastel/build/vs2008/lib"

-- The directory of the Boost library's source code.
-- The includes are of the form 'boost/static_assert.hpp'.
boostIncludeDir = "../../external/boost_1_45_0"

-- The directory of the SDL library's header files.
-- The includes are of the form 'SDL.h'.
sdlIncludeDir = "../../external/SDL-1.2.14/include"
sdlLibraryDir = "../../external/SDL-1.2.14/lib"

-- The directory of the Matlab header files.
-- The includes are of the form 'mex.h'.
matlabIncludeDir = "C:/Program Files/MATLAB/R2008a/extern/include"

-- No need to give a library path for Matlab:
-- Mex files are built from within Matlab.

outputDirectory = "build/" .. _ACTION

solution "Tim"

	language "C++"
	platforms
	{
		"native",
		"x32",
		"x64"
	}

	location(outputDirectory)
	targetdir(outputDirectory .. "/lib")

	configurations 
	{
		-- Debug information, with ASSERTs and PENSURES, without OpenMP.
		"debug",
		-- With ASSERTs and PENSUREs, without OpenMP.
		"develop",
		-- Without ASSERTs and PENSUREs, without OpenMP.
		"release-without-openmp",
		-- Without ASSERTs and PENSUREs, with OpenMP.
		"release"
	}

	flags
	{
		-- Treat wchar_t as a native type.
		"NativeWChar", 
		-- Do not give warnings from 64-bit checks.
		"No64BitChecks", 
		-- Do not use precompiled headers.
		"NoPCH"
	}
	
	configuration "debug"
		-- Debug libraries are suffixed with 'd'.
		targetsuffix "_d"
		-- Enable debug information.
		flags {"Symbols"}
		-- Enable ASSERTs and PENSUREs.
		defines 
		{
			"PASTEL_DEBUG_MODE", 
			"PASTEL_ENABLE_PENSURES"
		}
		
	configuration "develop"
		-- Developer libraries are suffixed with 'v'.
		targetsuffix "_v"
		-- Enable optimizations.
		flags {"Optimize"}
		-- Enable ASSERTs and PENSUREs.
		defines
		{
			"PASTEL_DEBUG_MODE", 
			"PASTEL_ENABLE_PENSURES",
		}

	configuration "release*"
		-- Release libraries do not have a suffix.
		-- Enable optimizations.
		flags {"Optimize"}

	-- Determine the SDL library name.	
	sdlLibrary = "SDL"
	
	includeDirectorySet = 
	{
		"./",
		boostIncludeDir,
		sdlIncludeDir,
		pastelIncludeDir,
		matlabIncludeDir
	}
	
	libraryDirectorySet =
	{
		sdlLibraryDir,
		pastelLibraryDir
	}
	
	fileSet = 
	{
		"*.cpp",
		"*.hpp",
		"*.h"
	}

	-- Visual Studio specific build options
	configuration "vs*"
		-- Disable warnings.
		buildoptions
		{
			-- "'identifier' : truncation from 'type1' to 'type2'"
			"/wd4305",
			-- "'argument' : conversion from 'type1' to 'type2', possible loss of data."
			"/wd4244",
			-- "'identifier' : class 'type' needs to have dll-interface to be used by clients of class 'type2'"
			"/wd4251", 
			-- "new behavior: elements of array 'array' will be default initialized"
			"/wd4351",
			-- "'function': was declared deprecated" (referring to STL functions)
			"/wd4996",
			-- "'var' : conversion from 'size_t' to 'type', possible loss of data"
			"/wd4267",
			-- "'expression' : signed/unsigned mismatch"
			"/wd4018", 
			-- "'operation' : conversion from 'type1' to 'type2' of greater size"
			"/wd4312",
			-- "nonstandard extension used : formal parameter 'identifier' was previously defined as a type"
			"/wd4224",
			-- "qualifier applied to function type has no meaning; ignored"
			"/wd4180",
			-- "'type' : forcing value to bool 'true' or 'false' (performance warning)"
			"/wd4800", 
			-- "'operation' : unsafe use of type 'bool' in operation"
			"/wd4804"
		}

		-- Disable Microsoft's Secure STL
		defines
		{
			"_SECURE_SCL=0",
			"_HAS_ITERATOR_DEBUGGING=0",
			"PASTEL_VISUAL_STUDIO"
		}
		
		-- Disable language extensions
		buildoptions
		{
			"/Za"			
		}
		
	-- GCC specific build options.
	configuration "gmake"
		-- Enables some additional warnings.
		buildoptions { "-Wall" }
		-- Disable some warnings.
		buildoptions 
		{ 
			-- Pragma warnings caused by OpenMP support not being enabled.
			"-Wno-unknown-pragmas", 
			-- Comparison between an unsigned and a signed integer.
			"-Wno-sign-compare", 
			-- Conversion between an unsigned and a signed integer.
			"-Wno-sign-conversion",
			-- Breaking strict aliasing rules.
			"-Wno-strict-aliasing",
			-- Compiler warns that it optimizes code based on the 
			-- assumption that signed integer overflows do not occur.
			"-Wno-strict-overflow"
		}
		
	-- Enable OpenMP if requested

	configuration "debug"
		defines 
		{ 
			"PASTEL_ENABLE_PENSURES",
			"PASTEL_DEBUG_MODE"
		}

	configuration "develop"
		defines 
		{ 
			"PASTEL_ENABLE_PENSURES",
			"PASTEL_DEBUG_MODE"
		}

	configuration "release"
		defines { "PASTEL_ENABLE_OMP" }

	configuration { "vs*",  "release" }
		buildoptions { "/openmp" }		

	configuration { "gmake",  "release" }
		buildoptions { "-fopenmp" }		
		links { "gomp" }

	function addPrefix(prefix, stringSet)
		resultSet = {}
		for i, name in pairs(stringSet)
		do
			resultSet[#resultSet + 1] = prefix .. name
		end
		return resultSet
	end

	libKind = "StaticLib"

	if buildTimCore
	then
		project "TimCore"
			kind(libKind)
			includedirs(includeDirectorySet)
			libdirs(libraryDirectorySet)
			files(addPrefix("tim/core/", fileSet))
	end

	if buildTimMatlab
	then	
		project "TimMatlab"
			kind(libKind)
			includedirs(includeDirectorySet)
			libdirs(libraryDirectorySet)
			files(addPrefix("tim/matlab/", fileSet))
	end

	if buildTimConsole
	then
		project "TimConsole"
			kind("ConsoleApp")
			includedirs(includeDirectorySet)
			libdirs(libraryDirectorySet)
			files(addPrefix("tim/console/", fileSet))
			links
			{
				"PastelGeometry",
				"PastelMath",
				"PastelSys",
				"TimCore"
			}
	end

	if buildTimCoreTest
	then
		project "TimCoreTest"
			kind("ConsoleApp")
			includedirs(includeDirectorySet)
			libdirs(libraryDirectorySet)
			files(addPrefix("test/coretest/", fileSet))
			links
			{
				"SDL",
				"PastelDevice",
				"PastelGfx",
				"PastelGeometry",
				"PastelDsp",
				"PastelMath",
				"PastelSys",
				"TimCore"
			}
	end