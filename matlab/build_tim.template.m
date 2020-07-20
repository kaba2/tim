% BUILD_TIM
% Builds Tim Matlab mex-libraries.
%
% build_tim()
% build_tim('key', value, ...)
%
% Optional arguments
% ------------------
%
% VERBOSE ('verbose') is a string which specifies whether to print 
% additional information about the build process. Must be either 'on' 
% or 'off'.
% Default: off

% Description: Builds Tim Matlab mex-libraries
% Documentation: building_timmatlab.txt

function build_tim(varargin)

% Add Pastel to Matlab path so that
% process_options is available.
addpath(['${PastelExecutableDirectory}', '/matlab']);

% Optional input arguments
mode = '${LOWER_CMAKE_BUILD_TYPE}';
verbose = 'off';
eval(pastelmatlab.process_options(...
    {'libraryName', 'mode', 'verbose'}, ...
    varargin));

% Determine the bitness of the running Matlab version
% (not of the operating system or the computer).
[ignore, maxArraySize] = computer;
bits = 32;
if maxArraySize >= 2^32
    bits = 64;
end

modeSet = {'debug', 'release', 'relwithdebinfo'};
if ~ismember(mode, modeSet)
    error(['MODE must be one of debug, release, or relwithdebinfo.']);
end

verboseSet = {'on', 'off'};
if ~ismember(verbose, verboseSet)
    error(['VERBOSE must be one on or off.']);
end

libraryName = 'core';
completeLibraryName = ['tim_matlab'];

inputDirectory = ['${TimDirectory}/tim/', libraryName, 'matlab'];
outputDirectory = '+tim';

% Directories
% -----------

defineSet = {};
includeDirectorySet = ...
{...
    '${TimIncludeDirectory}', ...
    '${PastelIncludeDirectory}', ...
    '${RangesIncludeDirectory}', ...
    '${Boost_INCLUDE_DIRS}', ...
    '${TbbIncludeDirectory}', ...
    '${EIGEN3_INCLUDE_DIRS}'...
};

libraryDirectorySet = ...
{...
    '${PastelLibraryDirectory}', ...
    '${TimLibraryDirectory}', ...
    '${TbbLibraryDirectory}' ...
};

% Libraries

librarySet = ...
{...
    '${TbbLibraryName}' ...
};

if strcmp(libraryName, 'core')
    librarySet = [...
        librarySet, ...
        'timcorematlab', ...
        'timcore', ...
        'pastel', ...
        'pastelmatlab' ...
    ];
end

fileSet = {[inputDirectory, '/tim_matlab.cpp']};

commandSet = form_build_command(...
    fileSet, ...
    'outputLibraryName', completeLibraryName, ...
    'includeDirectorySet', includeDirectorySet, ...
    'libraryDirectorySet', libraryDirectorySet, ...
    'librarySet', librarySet, ...
    'defineSet', defineSet, ...
    'outputDirectory', outputDirectory, ...
    'mode', mode, ...
    'verbose', verbose, ...
    'run', true);

end