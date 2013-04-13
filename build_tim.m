% BUILD_TIM
% Builds a Matlab mex-library for TIM.
%
% build_pastel(mode)
%
% where
%
% MODE is a string which specifies the configuration to use for building TIM. 
% It must be one of 
%     debug
%     develop
%     release (default)
%     release-without-openmp

% Description: Builds a Matlab mex-library for TIM.
% Documentation: building.txt

function build_tim(mode)

if nargin < 1
    mode = 'release';
end

modeSet = {'debug', 'develop', 'release', 'release-without-openmp'};
if ~ismember(mode, modeSet)
    error(['MODE must be one of debug, develop, release, ', ...
        'or release-without-openmp.']);
end

disp(' ');
disp(['Building a mex file for TIM in ', mode, ' mode.']);
disp(' ');

timIncludePath = '.';
pastelIncludePath = '..\pastel';
boostIncludePath = '..\external\boost_1_51_0';

timLibraryPath = ['build\vs2010\lib\', mode];
pastelLibraryPath = ['..\pastel\build\vs2010\lib\', mode];

inputPath = ['tim\matlab'];
outputPath = ['tim\+tim'];

% Paths
% -----

defineSet = {};
includePathSet = {};
libraryPathSet = {};

includePathSet{end + 1} = timIncludePath;
includePathSet{end + 1} = pastelIncludePath;
includePathSet{end + 1} = boostIncludePath;

libraryPathSet{end + 1} = timLibraryPath;
libraryPathSet{end + 1} = pastelLibraryPath;

% Libraries

librarySet = { ...
	'TimMatlab', ...
	'TimCore', ...
	'PastelGeometry', ...
	'PastelMath', ...
	'PastelMatlab', ...
	'PastelSys' ...
};

% Preprocessor definitions 
% ------------------------

defineSet{end + 1} = '_ITERATOR_DEBUG_LEVEL=0';

if strcmp(mode, 'release')
    defineSet{end + 1} = 'PASTEL_ENABLE_OMP';
end

if strcmp(mode, 'develop') || strcmp(mode, 'debug')
    defineSet{end + 1} = 'PASTEL_ENABLE_PENSURES';
end

% Form the build-command
% ----------------------

commandSet = {};
commandSet{end + 1} = 'mex';
commandSet{end + 1} = [' ', inputPath, '\tim_matlab.cpp'];

% Add preprocessor definitions.
for i = 1 : numel(defineSet)
    commandSet{end + 1} = [' -D', defineSet{i}];
end

% Add include paths.
for i = 1 : numel(includePathSet)
    commandSet{end + 1} = [' -I''', includePathSet{i}, ''''];
end

% Add library paths.
for i = 1 : numel(libraryPathSet)
    commandSet{end + 1} = [' -L''', libraryPathSet{i}, ''''];
end

% Add output path.
commandSet{end + 1} = [' -outdir ''', outputPath, ''''];

% Add libraries.
for i = 1 : numel(librarySet)
    commandSet{end + 1} = [' -l', librarySet{i}];
end

% Other flags.
if strcmp(mode, 'debug')
    commandSet{end + 1} = ' -g';
end

% Run the build-command
% ---------------------

% Print the executed build-command in pretty form.
for i = 1 : numel(commandSet)
    disp(commandSet{i});
end

% Run the command.
buildCommand = [commandSet{:}];
eval(buildCommand);