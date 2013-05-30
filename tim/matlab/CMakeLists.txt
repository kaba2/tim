project (TimMatlab)

file (GLOB SourceSet ${TimSourceGlobSet})
add_library (timmatlab ${TimLibraryType} ${SourceSet})

target_link_libraries (
	timmatlab
	timcore
	${PastelLibraryDirectory}/pastelsys 
	${PastelLibraryDirectory}/pastelmath
	${PastelLibraryDirectory}/pastelgeometry 
	${PastelLibraryDirectory}/pastelmatlab 
)
