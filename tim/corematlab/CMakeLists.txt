project (TimMatlab)

file (GLOB SourceSet ${TimSourceGlobSet})
source_group (Code FILES ${SourceSet})

file (GLOB DocumentationSet ${TimDocumentationGlobSet})
source_group (Documentation FILES ${DocumentationSet})

add_library (timmatlab ${TimLibraryType} ${SourceSet})

target_link_libraries (
	timmatlab
	timcore
	${PastelLibraryDirectory}/pastelsys 
	${PastelLibraryDirectory}/pastelmath
	${PastelLibraryDirectory}/pastelgeometry 
	${PastelLibraryDirectory}/pastelmatlab 
)
