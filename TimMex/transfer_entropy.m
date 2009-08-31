% TRANSFER_ENTROPY 
% A transfer entropy estimate from samples.
%
% I = transfer_entropy(X, Y, W, 
%       xLag, yLag, wLag, k, threads)
%
% where
%
% X, Y, and W are cell arrays of arbitrary dimension whose linearization
% contains q trials of the signals X, Y, and W, respectively. 
%
% XLAG, YLAG, and WLAG are the lags in samples applied to 
% signals X, Y, and W, respectively.
%
% K determines which k:th nearest neighbor the algorithm
% uses for estimation. Default 1.
%
% THREADS determines the number of threads to use for parallelization.
% To fully take advantage of multiple cores in your machine, set this
% to the number of cores in your machine. Note however that this makes 
% your computer unresponsive to other tasks. When you need responsiveness, 
% spare one core for other work. Default 1 (no parallelization).
%
% Each signal is a real (m x n)-matrix that contains n samples of an
% m-dimensional signal. The signals contained in X (or Y or W) 
% must all have equal dimensionality, but their number of samples may vary. 
% If the number of samples varies with trials, the function uses 
% the minimum sample count among the trials of X, Y, and W.
% The number of trials in X, Y, and W must be equal.

% Description: Transfer entropy estimation
% Documentation: tim_matlab.txt

function I = transfer_entropy(X, Y, W, ...
    xLag, yLag, wLag, k, threads)

if nargin < 3
    error('Not enough input arguments.');
end

if ~iscell(X) || ~iscell(Y) || ~iscell(W)
    error('X, Y, or W is not a cell-array.');
end

if nargin >= 4 && nargin < 6
    error('Lags must be specified all at once to avoid errors.');
end

if nargin < 4
    xLag = 0;
    yLag = 0;
    wLag = 0;
end

if nargin < 9
    k = 1;
end

if nargin < 10
    threads = 1;
end

I = entropy_combination(...
    [W(:), X(:), Y(:)]', ...
    [1, 2, 1; 2, 3, 1; 2, 2, -1], ...
    [wLag, xLag, yLag], k, threads);
