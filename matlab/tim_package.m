% TIM_PACKAGE
% Matlab-package the calling m-file is contained in.
%
% packageName = tim_package(...)
%
% where
%
% PACKAGENAME is a string which gives the name of the 
% package in which the calling m-file is contained in.
% If NAME is given, then [packagename, '.', name].
%
% Optional arguments
% ------------------
%
% A string which is concatenated to the package
% name in the form packageName.name.

% Description: Matlab-package the calling m-file is contained in.
% Documentation: tim_matlab.txt

function packageName = tim_package(varargin)

name = '';
if nargin >= 1
	name = varargin{1};
end

% Find out the m-file of the calling function.
callStack = dbstack('-completenames');
callerPath = callStack(2).file;

% Note: the following does not work. Don't know why exactly.
%callerPath = evalin('caller', 'mfilename(''fullpath'')')

packageNameSet = regexpi(callerPath, ...
    '\+([^\\/]*)', 'tokens', 'once');

packageName = '';
if ~isempty(packageNameSet)
	packageName = packageNameSet{1};
end

if ~isempty(name)
	packageName = strcat([packageName, '.', name]);
end
