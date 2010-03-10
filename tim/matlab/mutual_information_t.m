% MUTUAL_INFORMATION_T
% A temporal mutual information estimate from samples.
%
% I = mutual_information_t(X, Y, timeWindowRadius, 
%       xLag, yLag, k, filter, threads)
%
% where
%
% X and Y are signal sets.
%
% Type 'help tim_matlab' for more documentation.

% Description: Temporal mutual information estimation
% Documentation: tim_matlab_matlab.txt

function I = mutual_information_t(X, Y, timeWindowRadius, ...
    xLag, yLag, k, threads)

if nargin < 3
    error('Not enough input arguments.');
end

if nargin >= 4 && nargin < 5
    error('Lags must be specified all at once to avoid errors.');
end

if nargin < 4
    xLag = 0;
    yLag = 0;
end

if nargin < 6
    k = 1;
end

if nargin < 7
    threads = maxNumCompThreads;
end

if isnumeric(X)
    I = mutual_information_t({X}, Y, timeWindowRadius, ...
        xLag, yLag, k, threads);
    return
end

if isnumeric(Y)
    I = mutual_information_t(X, {Y}, timeWindowRadius, ...
        xLag, yLag, k, threads);
    return
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
    k, threads);
