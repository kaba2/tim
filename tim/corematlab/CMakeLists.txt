project (TimCoreMatlab)

EcAddLibrary (library timcorematlab "${TimSourceGlobSet}")

target_link_libraries (
	timcorematlab
	timcore
	pastelsys 
	pastelmath
	pastelgeometry 
	pastelmatlab 
)
