project (TimCore)

EcAddLibrary (library timcore "${TimSourceGlobSet}")

target_link_libraries (
	timcore
	${PastelLibraryDirectory}/pastelsys 
	${PastelLibraryDirectory}/pastelmath
	${PastelLibraryDirectory}/pastelgeometry 
)
