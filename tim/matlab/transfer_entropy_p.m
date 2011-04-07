% TRANSFER_ENTROPY_P
% A partial transfer entropy estimate from samples.
%
% I = transfer_entropy_p(X, Y, Z, W, 
%       xLag, yLag, zLag, wLag, k)
%
% where
%
% X, Y, Z, and W are signal sets.
%
% Type 'help tim' for more documentation.

% Description: Partial transfer entropy estimation
% Documentation: transfer_entropy.txt

function I = transfer_entropy_p(X, Y, Z, W, ...
    xLag, yLag, zLag, wLag, k)

concept_check(nargin, 'inputs', [4, 8, 9]);
concept_check(nargout, 'outputs', 0 : 1);

if nargin < 5
    xLag = 0;
    yLag = 0;
    zLag = 0;
    wLag = 0;
end

if nargin < 9
    k = 1;
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

I = entropy_combination(...
    [W(:)'; X(:)'; Z(:)'; Y(:)'], ...
    [1, 3, 1; 2, 4, 1; 2, 3, -1], ...
    {wLag, xLag, zLag, yLag}, ...
    k);
