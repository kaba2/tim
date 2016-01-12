% BUILD_TIM
% Builds Tim Matlab mex-libraries.
%
% build_tim()
% build_tim('key', value, ...)
%
% Optional arguments
% ------------------
%
% MODE ('mode') is a string which specifies the configuration to use for 
% building the Tim sub-library. It must be one of 
%     debug
%     release
%     relwithdebinfo
% Default: release
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
addpath(['${PastelDirectory}', '/matlab']);

% Optional input arguments
mode = 'release';
verbose = 'off';
eval(pastelsys.process_options(...
    {'libraryName', 'mode', 'verbose'}, ...
    varargin));

% Determine the bitness of the running Matlab version
% (not of the operating system or the computer).
[ignore, maxArraySize] = computer;
bits = 32;
if maxArraySize >= 2^32
    bits = 64;
end

generatedBits = ${GENERATOR_BITS};
if bits ~= generatedBits
    error(['Build file is ', num2str(generatedBits), '-bit, but ', ...
        'Matlab is ', num2str(bits), '-bit. ', ...
        'Switch to ', num2str(bits), '-bit Matlab, or rerun Tim CMake with a ', ...
        num2str(bits), '-bit generator.'])
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

disp(' ');
disp(['Building a mex file for ', ...
    'tim', libraryName, 'matlab', ...
    ' in ', ...
    num2str(bits), '-bit ', ...
    mode, ' mode.']);
disp(' ');

inputDirectory = ['${TimDirectory}/tim/', libraryName, 'matlab'];
outputDirectory = '+tim'

% Directories
% -----------

defineSet = {};
includeDirectorySet = ...
{...
    '${TimIncludeDirectory}', ...
    '${PastelIncludeDirectory}', ...
    '${BoostIncludeDirectory}', ...
    '${TbbIncludeDirectory}', ...
    '${ArmadilloIncludeDirectory}'...
};

libraryDirectorySet = ...
{...
    ['${PastelLibraryDirectory}/', mode], ...
    ['${TimLibraryDirectory}/', mode], ...
    '${TbbLibraryDirectory}', ...
    '${ArmadilloLibraryDirectory}', ...
};

if ~ismac()
    libraryDirectorySet = ...
    [...
        libraryDirectorySet, ...
        '${LapackLibraryDirectory}', ...
        '${BlasLibraryDirectory}'
    ];
end

% Libraries

librarySet = ...
{...
    '${TbbLibraryName}', ...
    '${ArmadilloLibraryName}' ...
};

if ~ismac()
    librarySet = ...
    [...
        librarySet, ...
        '${BlasLibraryName}', ...
        '${LapackLibraryName}' ...
    ];
end

if strcmp(libraryName, 'core')
    librarySet = [...
        librarySet, ...
        'timcorematlab', ...
        'timcore', ...
        'pastelgeometry', ...
        'pastelmath', ...
        'pastelmatlab', ...
        'pastelsys' ...
    ];
end

% Handle semicolon-separated lists
% --------------------------------

function resultSet = breakSemicolons(stringSet)
    resultSet = {};
    for i = 1 : numel(stringSet)
        string = stringSet{i};
        splitSet = regexp(string, ';', 'split');
        resultSet = [resultSet, splitSet];
    end
end

includeDirectorySet = breakSemicolons(includeDirectorySet);
libraryDirectorySet = breakSemicolons(libraryDirectorySet);
librarySet = breakSemicolons(librarySet);

% Preprocessor definitions 
% ------------------------

defineSet{end + 1} = '_ITERATOR_DEBUG_LEVEL=0';

% Form the build-command
% ----------------------

commandSet = {};
commandSet{end + 1} = 'mex';
commandSet{end + 1} = [' ', inputDirectory, '/tim_matlab.cpp'];

% Add preprocessor definitions.
for i = 1 : numel(defineSet)
    commandSet{end + 1} = [' -D', defineSet{i}];
end

% Add include paths.
for i = 1 : numel(includeDirectorySet)
    commandSet{end + 1} = [' -I''', includeDirectorySet{i}, ''''];
end

% Add library directories.
for i = 1 : numel(libraryDirectorySet)
    libraryDirectory = libraryDirectorySet{i};
    if isempty(libraryDirectory)
        % No directory was specified; skip it.
        continue;
    end

    [iStart, iEnd] = regexp(libraryDirectory, '\.framework$');
    if ~isempty(iStart) && ismac()
        % This is a framework directory.
        commandSet{end + 1} = [' LDFLAGS=''$LDFLAGS -F ', libraryDirectory, ''''];
    else
        % This is a normal library directory.
        commandSet{end + 1} = [' -L''', libraryDirectory, ''''];
    end
end

% Add output path.
commandSet{end + 1} = [' -outdir ''', outputDirectory, ''''];

% Add libraries.
for i = 1 : numel(librarySet)
    library = librarySet{i};
    if isempty(library)
        % No library was specified; skip it.
        continue;
    end

    [iStart, iEnd] = regexp(library, '\.framework$');
    if ~isempty(iStart) && ismac()
        % This is a framework (Mac OS X).
        
        % Get the framework name by stripping out the
        % trailing '.framework'.
        frameworkName = library(1 : iStart - 1);

        % Mex does not have a direct flag to specify
        % frameworks. Instead, we embed it directly
        % into compiler flags.
        commandSet{end + 1} = [' CXXLIBS=''$CXXLIBS -framework ', frameworkName, ''''];
    else
        % This is a normal library.
        commandSet{end + 1} = [' -l', librarySet{i}];
    end
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

end