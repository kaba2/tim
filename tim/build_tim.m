% BUILD_TIM
% Builds TIM Matlab mex-library.
%
% build_tim()
% build_tim('key', value, ...)
%
% Optional arguments
% ------------------
%
% MODE ('mode') is a string which specifies the configuration to use for 
% building the TIM library. It must be one of 
%     debug
%     release
%     relwithdebinfo
% Default: release
%
% VERBOSE ('verbose') is a string which specifies whether to print 
% additional information about the build process. Must be either 'on' 
% or 'off'.
% Default: off

% Description: Builds a Matlab mex-library for TIM.
% Documentation: building.txt

function build_tim(varargin)

% Path to Pastel.
pastelPath = '../../pastel';

% Path to Boost.
boostPath = '../../boost_1_56_0';

% Path to Threading Building Blocks.
tbbPath = '../../tbb42';

% Add Pastel to Matlab path so that
% process_options is available.
addpath([pastelPath, '/pastel']);

% Optional input arguments
mode = 'release';
verbose = 'off';
eval(pastelsys.process_options(...
    {'mode', 'verbose'}, ...
    varargin));

modeSet = {'debug', 'release', 'relwithdebinfo'};
if ~ismember(mode, modeSet)
    error(['MODE must be one of debug, release, or relwithdebinfo.']);
end

verboseSet = {'on', 'off'};
if ~ismember(verbose, verboseSet)
    error(['VERBOSE must be one on or off.']);
end

% Determine the bitness of the running Matlab version
% (not of the operating system or the computer).
[ignore, maxArraySize] = computer;
bits = 32;
if maxArraySize >= 2^32
    bits = 64;
end

disp(' ');
disp(['Building a mex file for TIM in ', ...
    num2str(bits), '-bit ', ...
    mode, ' mode.']);
disp(' ');

% Compute the actual paths.

pastelIncludePath = pastelPath;
pastelLibraryPath = [...
    '../../pastel/lib', ...
    '/', num2str(bits), ...
    '/', mode];

timIncludePath = '..';
timLibraryPath = [...
    '../lib', ...
    '/', num2str(bits), ...
    '/', mode];

boostIncludePath = boostPath;

tbbIncludePath = [tbbPath, '/include'];
tbbLibraryPath = [tbbPath, '/lib'];

if ispc()
    % Windows
    if bits == 64
        tbbLibraryPath = [tbbLibraryPath, '/intel64/vc11'];
    else
        tbbLibraryPath = [tbbLibraryPath, '/ia32/vc11'];
    end
elseif ismac()
    % Mac Os X
    tbbLibraryPath = [tbbLibraryPath, '/libc++'];
else
    % Linux and others
    tbbLibraryPath = [tbbLibraryPath, ''];
end

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

if strcmp(verbose, 'on')
    commandSet{end + 1} = ' -v';
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
