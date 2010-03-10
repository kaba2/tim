% TRANSFER_ENTROPY 
% A transfer entropy estimate from samples.
%
% I = transfer_entropy(X, Y, W, 
%       xLag, yLag, wLag, k, threads)
%
% where
%
% X, Y, and W are signal sets.
%
% Type 'help tim_matlab' for more documentation.

% Description: Transfer entropy estimation
% Documentation: tim_matlab_matlab.txt

function I = transfer_entropy(X, Y, W, ...
    xLag, yLag, wLag, k, threads)

if nargin < 3
    error('Not enough input arguments.');
end

if nargin >= 4 && nargin < 6
    error('Lags must be specified all at once to avoid errors.');
end

if nargin < 4
    xLag = 0;
    yLag = 0;
    wLag = 0;
end

if nargin < 7
    k = 1;
end

if nargin < 8
    threads = maxNumCompThreads;
end

if isnumeric(X)
    I = transfer_entropy({X}, Y, W, ...
        xLag, yLag, wLag, k, threads);
    return
end

if isnumeric(Y)
    I = transfer_entropy(X, {Y}, W, ...
        xLag, yLag, wLag, k, threads);
    return
end

if isnumeric(W)
    I = transfer_entropy(X, Y, {W}, ...
        xLag, yLag, wLag, k, threads);
    return
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
    k, threads);
