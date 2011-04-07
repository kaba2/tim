% MUTUAL_INFORMATION_T
% A temporal mutual information estimate from samples.
%
% I = mutual_information_t(X, Y, timeWindowRadius, 
%       xLag, yLag, k, filter)
%
% where
%
% X and Y are signal sets.
%
% Type 'help tim' for more documentation.

% Description: Temporal mutual information estimation
% Documentation: mutual_information.txt

function I = mutual_information_t(X, Y, timeWindowRadius, ...
    xLag, yLag, k, filter)

concept_check(nargin, 'inputs', [3, 5, 6, 7]);
concept_check(nargout, 'outputs', 0 : 1);

if nargin < 4
    xLag = 0;
    yLag = 0;
end

if nargin < 6
    k = 1;
end

if nargin < 7
    filter = 1;
end

if isnumeric(X)
    X = {X};
end

if isnumeric(Y)
    Y = {Y};
end

if ~iscell(X) || ~iscell(Y)
    error('X or Y is not a cell-array.');
end

if numel(X) ~= numel(Y)
    error('The number of trials in X and Y differ.');
end

% Pass parameter error checking to entropy_combination.

I = entropy_combination_t(...
    [X(:)'; Y(:)'], ...
    [1, 1, 1; 2, 2, 1], timeWindowRadius, ...
    {xLag, yLag}, ...
    k, filter);
