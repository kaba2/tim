% TRANSFER_ENTROPY_T
% A temporal transfer entropy estimate from samples.
%
% I = transfer_entropy_t(X, Y, W, 
%       timeWindowRadius, xLag, yLag, wLag, k, filter, threads)
%
% where
%
% X, Y, and W are signal sets.
%
% Type 'help tim' for more documentation.

% Description: Temporal transfer entropy estimation
% Documentation: tim_matlab_matlab.txt

function I = transfer_entropy_t(X, Y, W, ...
    timeWindowRadius, xLag, yLag, wLag, k, filter, threads)

if nargin < 4
    error('Not enough input arguments.');
end

if nargin > 10
    error('Too many input arguments.');
end

if nargin >= 5 && nargin < 7
    error('Lags must be specified all at once to avoid errors.');
end

if nargin < 5
    xLag = 0;
    yLag = 0;
    wLag = 0;
end

if nargin < 8
    k = 1;
end

if nargin < 9
    filter = [1];
end

if nargin < 10
    threads = maxNumCompThreads;
end

if isnumeric(X)
    I = transfer_entropy_t({X}, Y, W, timeWindowRadius, ...
        xLag, yLag, wLag, k, filter, threads);
    return
end

if isnumeric(Y)
    I = transfer_entropy_t(X, {Y}, W, timeWindowRadius, ...
        xLag, yLag, wLag, k, filter, threads);
    return
end

if isnumeric(W)
    I = transfer_entropy_t(X, Y, {W}, timeWindowRadius, ...
        xLag, yLag, wLag, k, filter, threads);
    return
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
    {wLag, xLag, yLag}, k, filter, threads);
