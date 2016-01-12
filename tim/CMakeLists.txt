if (BuildLibraries)
	add_subdirectory ("core")
	#add_subdirectory ("console")
endif()

if (BuildMatlab)
	add_subdirectory ("corematlab")
endif()
