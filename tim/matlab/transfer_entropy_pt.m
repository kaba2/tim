% TRANSFER_ENTROPY_PT
% A temporal partial transfer entropy estimate from samples.
%
% I = transfer_entropy_pt(X, Y, Z, W, 
%       timeWindowRadius, xLag, yLag, zLag, wLag, k, filter)
%
% where
%
% X, Y, Z, and W are signal sets.
%
% Type 'help tim' for more documentation.

% Description: Temporal partial transfer entropy estimation
% Documentation: tim_matlab_matlab.txt

function I = transfer_entropy_pt(X, Y, Z, W, ...
    timeWindowRadius, xLag, yLag, zLag, wLag, k, filter)

check(nargin, 'inputs', [5, 9, 10, 11]);
check(nargout, 'outputs', 0 : 1);

if nargin < 6
    xLag = 0;
    yLag = 0;
    zLag = 0;
    wLag = 0;
end

if nargin < 10
    k = 1;
end

if nargin < 11
    filter = [1];
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

if isnumeric(W)
    W = {W};
end

if ~iscell(X) || ~iscell(Y) || ~iscell(Z) || ~iscell(W)
    error('X, Y, Z, or W is not a cell-array.');
end

if numel(X) ~= numel(Y) || numel(X) ~= numel(Z) || ...
    numel(X) ~= numel(W)
    error('The number of trials in X, Y, Z, and W differ.');
end

% Pass parameter error checking to entropy_combination.

I = entropy_combination_t(...
    [W(:)'; X(:)'; Z(:)'; Y(:)'], ...
    [1, 3, 1; 2, 4, 1; 2, 3, -1], ...
    timeWindowRadius, ...
    {wLag, xLag, zLag, yLag}, ...
    k, filter);
