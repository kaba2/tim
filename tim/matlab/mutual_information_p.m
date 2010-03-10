% MUTUAL_INFORMATION_P
% A partial mutual information estimate from samples.
%
% I = mutual_information_p(X, Y, Z, xLag, yLag, zLag, k, threads)
%
% where
%
% X, Y, and Z are signal sets.
%
% Type 'help tim_matlab' for more documentation.

% Description: Partial mutual information estimation
% Documentation: tim_matlab_matlab.txt

function I = mutual_information_p(...
    X, Y, Z, xLag, yLag, zLag, k, threads)

if nargin < 3
    error('Not enough input arguments.');
end

if nargin >= 4 && nargin < 6
    error('Lags must be specified all at once to avoid errors.');
end

if nargin < 4
    xLag = 0;
    yLag = 0;
    zLag = 0;
end

if nargin < 7
    k = 1;
end

if nargin < 8
    threads = maxNumCompThreads;
end

if isnumeric(X)
    I = mutual_information_p({X}, Y, Z, xLag, yLag, zLag, k, threads);
    return
end

if isnumeric(Y)
    I = mutual_information_p(X, {Y}, Z, xLag, yLag, zLag, k, threads);
    return
end

if isnumeric(Z)
    I = mutual_information_p(X, Y, {Z}, xLag, yLag, zLag, k, threads);
    return
end

if ~iscell(X) || ~iscell(Y) || ~iscell(Z)
    error('X, Y, or Z is not a cell-array.');
end

if numel(X) ~= numel(Y) || numel(X) ~= numel(Z)
    error('The number of trials in X, Y, and Z differ.');
end

% Pass parameter error checking to entropy_combination.

I = entropy_combination(...
    [X(:)'; Z(:)'; Y(:)'], ...
    [1, 2, 1; 2, 3, 1; 2, 2, -1], ...
    {xLag, zLag, yLag}, ...
    k, threads);
