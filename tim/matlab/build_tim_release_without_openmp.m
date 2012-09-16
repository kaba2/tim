% Description: TIM Matlab interface release-without-openmp build script
% Documentation: building_matlab.txt

mex tim_matlab.cpp ...
	-D_HAS_ITERATOR_DEBUGGING=0 ...
	-I'../../../pastel' ...
	-I'../../../tim' ...
	-I'../../../../external/boost_1_49_0' ...
	-L'../../../pastel/build/vs2010/lib/release_without_openmp' ...
	-L'../../../tim/build/vs2010/lib/release_without_openmp' ...
	-lTimMatlab ...
	-lTimCore ...
	-lPastelGeometry ...
	-lPastelMath ...
	-lPastelMatlab ...
	-lPastelSys

