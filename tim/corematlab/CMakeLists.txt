project (TimCoreMatlab)

EcAddLibrary (library timcorematlab "${TimSourceGlobSet}")

target_link_libraries (
	timcorematlab
	timcore
	${PastelLibraryDirectory}/pastelsys 
	${PastelLibraryDirectory}/pastelmath
	${PastelLibraryDirectory}/pastelgeometry 
	${PastelLibraryDirectory}/pastelmatlab 
)
