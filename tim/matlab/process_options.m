% PROCESS_OPTIONS
% A function pre-preprocessing optional input options. This function is
% required by multiple functions of TIM Matlab interface
%
% [command, remove] = process_options(keySet, optionSet)
%
% where
%
% KEYSET is a cell array with the names of the optional input options, e.g.
% {'threads', 'lags'}
%
% OPTIONSET is the varargin argument at the calling function, e.g.
% {'threads', 4, 'lags', [1 10 20]}
%
% COMMAND is the string that should be evaluated by the calling function in
% order to assign the provided values to the corresponding argument names,
% e.g. COMMAND would be 'threads=optionSet{4};lags={optionSet{2};' if
% OPTIONSET = {'lags', [1 10 20], 'threads', 4}
%
% Type 'help tim' for more documentation

% Description: Input options pre-processing
% Documentation: tim_matlab_impl.txt

function [command, remove] = process_options(keySet, optionSet)

command = '';

if nargin < 2,
    remove = [];
    return;
end

options = length(optionSet);
remove = false(1, options);
i = 1;
while i < options
    if ~ischar(optionSet{i}),
        error('MISC:process_options:invalidInput', ...
            'Optional input arguments must be given in 'key'-value pairs.');
    end
    
    [~, loc] = ismember(optionSet{i}, keySet);
    if loc > 0 && ~all(isempty(optionSet{i + 1}))
        command = [command keySet{loc} '=optionSet{' num2str(i+1)  '};']; %#ok<AGROW>
        remove(i:i + 1) = true;
    end
    
    i = i + 2;
end
