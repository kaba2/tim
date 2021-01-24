project (AnotherTest)

EcAddLibrary (executable anothertest "${TimSourceGlobSet}")

target_link_libraries (
    anothertest
	timcore
	pastel
)

