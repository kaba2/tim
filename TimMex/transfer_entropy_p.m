% TRANSFER_ENTROPY_P
% A multivariate transfer entropy estimate from samples.
%
% I = transfer_entropy_p(X, Y, Z, W, 
%       xLag, yLag, zLag, wLag, k, threads)
%
% where
%
% X, Y, Z, and W are cell arrays of arbitrary dimension whose linearization
% contains q trials of the signals X, Y, Z, and W, respectively. 
%
% XLAG, YLAG, ZLAG, and WLAG are the lags in samples applied to 
% signals X, Y, Z, and W, respectively.
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
% m-dimensional signal. The signals need not have the same dimension,
% but the dimension must be the same among the trials of a given signal.
% If the number of samples varies with each signal, the function uses 
% the minimum sample count among the signals. 

function I = transfer_entropy_p(X, Y, Z, W, ...
    xLag, yLag, zLag, wLag, k, threads)

if nargin < 4
    error('Not enough input arguments.');
end

if ~iscell(X) || ~iscell(Y) || ~iscell(Z) || ~iscell(W)
    error('X, Y, Z, or W is not a cell-array.');
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
    threads = 1;
end

I = entropy_combination(...
    [W(:), X(:), Z(:), Y(:)]', ...
    [1, 3, 1; 2, 4, 1; 2, 3, -1], ...
    [wLag, xLag, zLag, yLag], k, threads);
