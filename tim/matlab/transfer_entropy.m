% TRANSFER_ENTROPY 
% A transfer entropy estimate from samples.
%
% I = transfer_entropy(X, Y, W, 
%       xLag, yLag, wLag, k)
%
% where
%
% X, Y, and W are signal sets.
%
% Type 'help tim' for more documentation.

% Description: Transfer entropy estimation
% Documentation: transfer_entropy.txt

function I = transfer_entropy(X, Y, W, ...
    xLag, yLag, wLag, k)

check(nargin, 'inputs', [3, 6, 7]);
check(nargout, 'outputs', 0 : 1);

if nargin < 4
    xLag = 0;
    yLag = 0;
    wLag = 0;
end

if nargin < 7
    k = 1;
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

% Pass parameter error checking to entropy_combination.

I = entropy_combination(...
    [W(:)'; X(:)'; Y(:)'], ...
    [1, 2, 1; 2, 3, 1; 2, 2, -1], ...
    {wLag, xLag, yLag}, ...
    k);
