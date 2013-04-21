% TIM_PACKAGE
% Matlab-package the calling m-file is contained in.
%
% tim_package()
%
% where
%
% PACKAGENAME is a string which gives the name of the 
% package in which the calling m-file is contained in.

% Description: Matlab-package the calling m-file is contained in.
% Documentation: tim_matlab.txt

function packageName = tim_package()

% Find out the m-file of the calling function.
callStack = dbstack('-completenames');
callerPath = callStack(2).file;

% Note: the following does not work. Don't know why exactly.
%callerPath = evalin('caller', 'mfilename(''fullpath'')')

packageNameSet = regexpi(callerPath, ...
    '\+(tim[^\\/]*)', 'tokens', 'once');

packageName = '';
if ~isempty(packageNameSet)
	packageName = packageNameSet{1};
end
