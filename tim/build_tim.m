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
%     release
%     relwithdebinfo
% Default: release

% Description: Builds a Matlab mex-library for TIM.
% Documentation: building.txt

function build_tim(mode)

if nargin < 1
    mode = 'release';
end

modeSet = {'debug', 'release', 'relwithdebinfo'};
if ~ismember(mode, modeSet)
    error(['MODE must be either debug, release, or relwithdebinfo.']);
end

disp(' ');
disp(['Building a mex file for TIM in ', mode, ' mode.']);
disp(' ');

timIncludePath = '..';
pastelIncludePath = '../../pastel';
boostIncludePath = '../../boost';
tbbIncludePath = '../../tbb42/include';

timLibraryPath = ['../lib/', mode];
pastelLibraryPath = ['../../pastel/lib/', mode];
tbbLibraryPath = ['../../tbb42/lib/intel64/vc11'];

inputPath = ['matlab'];
outputPath = ['+tim'];

tbbName = 'tbb';
if strcmp(mode, 'debug')
    tbbName = 'tbb_debug';
end

% Paths
% -----

defineSet = {};
includePathSet = {...
    timIncludePath, ...
    pastelIncludePath, ...
    boostIncludePath, ...
    tbbIncludePath};
libraryPathSet = {...
    timLibraryPath, ...
    pastelLibraryPath, ...
    tbbLibraryPath};

% Libraries

librarySet = { ...
    tbbName, ...
	'timmatlab', ...
	'timcore', ...
	'pastelgeometry', ...
	'pastelmath', ...
	'pastelmatlab', ...
	'pastelsys' ...
};

% Preprocessor definitions 
% ------------------------

defineSet{end + 1} = '_ITERATOR_DEBUG_LEVEL=0';

% Form the build-command
% ----------------------

commandSet = {};
commandSet{end + 1} = 'mex';
commandSet{end + 1} = [' ', inputPath, '/tim_matlab.cpp'];

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
