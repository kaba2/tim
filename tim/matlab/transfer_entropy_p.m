% TRANSFER_ENTROPY_P
% A partial transfer entropy estimate from samples.
%
% I = transfer_entropy_p(X, Y, Z, W, 
%       xLag, yLag, zLag, wLag, k, threads)
%
% where
%
% X, Y, Z, and W are signal sets.
%
% Type 'help tim' for more documentation.

% Description: Partial transfer entropy estimation
% Documentation: tim_matlab_matlab.txt

function I = transfer_entropy_p(X, Y, Z, W, ...
    xLag, yLag, zLag, wLag, k, threads)

if nargin < 4
    error('Not enough input arguments.');
end

if nargin > 10
    error('Too many input arguments.');
end

if nargin >= 5 && nargin < 8
    error('Lags must be specified all at once to avoid errors.');
end

if nargin < 5
    xLag = 0;
    yLag = 0;
    zLag = 0;
    wLag = 0;
end

if nargin < 9
    k = 1;
end

if nargin < 10
    threads = maxNumCompThreads;
end

if isnumeric(X)
    I = transfer_entropy_p({X}, Y, Z, W, ...
        xLag, yLag, zLag, wLag, k, threads);
    return
end

if isnumeric(Y)
    I = transfer_entropy_p(X, {Y}, Z, W, ...
        xLag, yLag, zLag, wLag, k, threads);
    return
end

if isnumeric(Z)
    I = transfer_entropy_p(X, Y, {Z}, W, ...
        xLag, yLag, zLag, wLag, k, threads);
    return
end

if isnumeric(W)
    I = transfer_entropy_p(X, Y, Z, {W}, ...
        xLag, yLag, zLag, wLag, k, threads);
    return
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
    k, threads);
