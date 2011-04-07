% MUTUAL_INFORMATION 
% A mutual information estimate from samples.
%
% I = mutual_information(X, Y, xLag, yLag, k)
%
% where
%
% X and Y are signal sets.
%
% Type 'help tim' for more documentation.

% Description: Mutual information estimation
% Documentation: mutual_information.txt

function I = mutual_information(X, Y, xLag, yLag, k)

concept_check(nargin, 'inputs', [2, 4, 5]);
concept_check(nargout, 'outputs', 0 : 1);

if nargin < 3
    xLag = 0;
    yLag = 0;
end

if nargin < 5
    k = 1;
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

I = entropy_combination(...
    [X(:)'; Y(:)'], ...
    [1, 1, 1; 2, 2, 1], ...
    {xLag, yLag}, ...
    k);
