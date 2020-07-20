project (TimCore)

EcAddLibrary (library timcore "${TimSourceGlobSet}")

target_link_libraries (
	timcore
	pastel
)
