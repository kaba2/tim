% MUTUAL_INFORMATION_PT
% A temporal partial mutual information estimate from samples.
%
% I = mutual_information_pt(
%         X, Y, Z, timeWindowRadius, 
%         xLag, yLag, zLag, k, filter)
%
% where
%
% X, Y, and Z are signal sets.
%
% Type 'help tim' for more documentation.

% Description: Temporal partial mutual information estimation
% Documentation: mutual_information.txt

function I = mutual_information_pt(...
    X, Y, Z, timeWindowRadius, xLag, yLag, zLag, k, filter)

check(nargin, 'inputs', [4, 7, 8, 9]);
check(nargout, 'outputs', 0 : 1);

if nargin < 5
    xLag = 0;
    yLag = 0;
    zLag = 0;
end

if nargin < 8
    k = 1;
end

if nargin < 9
    filter = 1;
end

if isnumeric(X)
    X = {X};
end

if isnumeric(Y)
    Y = {Y};
end

if isnumeric(Z)
    Z = {Z};
end

if ~iscell(X) || ~iscell(Y) || ~iscell(Z)
    error('X, Y, or Z is not a cell-array.');
end

if numel(X) ~= numel(Y) || numel(X) ~= numel(Z)
    error('The number of trials in X, Y, and Z differ.');
end

% Pass parameter error checking to entropy_combination.

I = entropy_combination_t(...
    [X(:)'; Z(:)'; Y(:)'], ...
    [1, 2, 1; 2, 3, 1; 2, 2, -1], timeWindowRadius, ...
    {xLag, zLag, yLag}, ...
    k, filter);
