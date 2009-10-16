% MUTUAL_INFORMATION_P
% A partial mutual information estimate from samples.
%
% I = mutual_information_p(X, Y, Z, xLag, yLag, zLag, k, threads)
%
% where
%
% X, Y, and Z are arbitrary-dimensional cell-arrays whose 
% linearizations contain q trials of signal X, Y, and Z, 
% respectively.
%
% XLAG, YLAG and ZLAG are the lags in samples applied to 
% signals X, Y and Z, respectively.
%
% K determines which k:th nearest neighbor the algorithm
% uses for estimation. Default 1.
%
% THREADS determines the number of threads to use for parallelization.
% To fully take advantage of multiple cores in your machine, set this
% to the number of cores in your machine. Note however that this makes 
% your computer unresponsive to other tasks. When you need responsiveness, 
% spare one core for other work. Default maxNumCompThreads.
%
% Each signal is a real (m x n)-matrix that contains n samples of an
% m-dimensional signal. The signals contained in X (or Y or Z) must all 
% have equal dimensionality, but their number of samples may vary. 
% If the number of samples varies with trials, the function uses 
% the minimum sample count among the trials of X, Y, and Z.
% The number of trials in X, Y, and Z must be equal.

% Description: Partial mutual information estimation
% Documentation: tim_matlab.txt

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

if ~iscell(X) || ~iscell(Y) || ~iscell(Z)
    error('X, Y, or Z is not a cell-array.');
end

% Pass parameter error checking to entropy_combination.

I = entropy_combination(...
    [X(:), Z(:), Y(:)]', ...
    [1, 2, 1; 2, 3, 1; 2, 2, -1], [xLag, zLag, yLag], k, threads);
