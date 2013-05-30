project (TimCore)

file (GLOB SourceSet ${TimSourceGlobSet})
add_library (timcore ${TimLibraryType} ${SourceSet})

target_link_libraries (
	timcore
	${PastelLibraryDirectory}/pastelsys 
	${PastelLibraryDirectory}/pastelmath
	${PastelLibraryDirectory}/pastelgeometry 
)
