% MUTUAL_INFORMATION 
% A mutual information estimate from samples.
%
% I = mutual_information(X, Y, xLag, yLag, k, threads)
%
% where
%
% X and Y are signal sets.
%
% Type 'help tim' for more documentation.

% Description: Mutual information estimation
% Documentation: tim_matlab_matlab.txt

function I = mutual_information(X, Y, xLag, yLag, k, threads)

if nargin < 2
    error('Not enough input arguments.');
end

if nargin > 6
    error('Too many input arguments.');
end
 
if nargin >= 3 && nargin < 4
    error('Lags must be specified all at once to avoid errors.');
end

if nargin < 3
    xLag = 0;
    yLag = 0;
end

if nargin < 5
    k = 1;
end

if nargin < 6
    threads = maxNumCompThreads;
end

if isnumeric(X)
    I = mutual_information({X}, Y, xLag, yLag, k, threads);
    return
end

if isnumeric(Y)
    I = mutual_information(X, {Y}, xLag, yLag, k, threads);
    return
end

if ~iscell(X) || ~iscell(Y)
    error('X or Y is not a cell-array.');
end

if numel(X) ~= numel(Y)
    error('The number of trials in X and Y differ.');
end

% Pass parameter error checking to entropy_combination.

I = entropy_combination(...
    [X(:)'; Y(:)'], ...
    [1, 1, 1; 2, 2, 1], ...
    {xLag, yLag}, ...
    k, threads);
