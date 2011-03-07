% TRANSFER_ENTROPY_T
% A temporal transfer entropy estimate from samples.
%
% I = transfer_entropy_t(X, Y, W, 
%       timeWindowRadius, xLag, yLag, wLag, k, filter)
%
% where
%
% X, Y, and W are signal sets.
%
% Type 'help tim' for more documentation.

% Description: Temporal transfer entropy estimation
% Documentation: transfer_entropy.txt

function I = transfer_entropy_t(X, Y, W, ...
    timeWindowRadius, xLag, yLag, wLag, k, filter)

check(nargin, 'inputs', [4, 7, 8, 9]);
check(nargout, 'outputs', 0 : 1);

if nargin < 5
    xLag = 0;
    yLag = 0;
    wLag = 0;
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

if isnumeric(W)
    W = {W};
end

if ~iscell(X) || ~iscell(Y) || ~iscell(W)
    error('X, Y, or W is not a cell-array.');
end

if numel(X) ~= numel(Y) || numel(X) ~= numel(W)
    error('The number of trials in X, Y, and W differ.');
end

I = entropy_combination_t(...
    [W(:)'; X(:)'; Y(:)'], ...
    [1, 2, 1; 2, 3, 1; 2, 2, -1], ...
    timeWindowRadius, ...
    {wLag, xLag, yLag}, k, filter);
