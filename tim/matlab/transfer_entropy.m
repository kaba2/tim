% TRANSFER_ENTROPY 
% A transfer entropy estimate from samples.
%
% I = transfer_entropy(X, Y, W, 
%       xLag, yLag, wLag, k, threads)
%
% where
%
% X, Y, and W are cell-arrays of arbitrary dimension whose linearization
% contains q trials of the signals X, Y, and W, respectively. 
%
% XLAG, YLAG, and WLAG are the lags in samples applied to 
% signals X, Y, and W, respectively. Each can be given either as a scalar
% or as an array. In case some of the lags are given as arrays, those 
% arrays must have the same number of elements, and a scalar lag is 
% interpreted as an array of the same size with the given value as 
% elements. Default 0.
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
% I is a real (L x 1)-matrix of computed transfer entropies, where L is 
% the number of specified lags. The I(i) corresponds to the transfer
% entropy estimate using the lags XLAG(i), YLAG(i), and ZLAG(i).
%
% Each signal is a real (m x n)-matrix that contains n samples of an
% m-dimensional signal. The signals contained in X (or Y or W) 
% must all have equal dimensionality, but their number of samples may vary. 
% If the number of samples varies with trials, the function uses 
% the minimum sample count among the trials of X, Y, and W.
% The number of trials in X, Y, and W must be equal.

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
