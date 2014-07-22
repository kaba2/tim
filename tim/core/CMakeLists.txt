project (TimCore)

file (GLOB SourceSet ${TimSourceGlobSet})
source_group (Code FILES ${SourceSet})

file (GLOB DocumentationSet ${TimDocumentationGlobSet})
source_group (Documentation FILES ${DocumentationSet})

add_library (timcore ${TimLibraryType} ${SourceSet})

target_link_libraries (
	timcore
	${PastelLibraryDirectory}/pastelsys 
	${PastelLibraryDirectory}/pastelmath
	${PastelLibraryDirectory}/pastelgeometry 
)
