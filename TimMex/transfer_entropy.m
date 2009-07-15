% TRANSFER_ENTROPY 
% A multivariate transfer entropy estimate from samples.
%
% I = transfer_entropy(X, W, Y, Z, sigma, k, threads)
%
% where
%
% X is a cell-array of arbitrary dimension which contains q trials of 
% the already embedded X-signal. Call the signal from which the X-signal
% was embedded from the x-signal.
%
% W is a cell-array of arbitrary dimension which contains q trials of
% the time-shifted x-signal, i.e., the future of the x-signal by one sample.
%
% Y is a cell-array of arbitrary dimension which contains q trials of 
% the already embedded Y-signal.
%
% Z is a 2-dimensional p x q cell-array which contains q trials 
% of each of the p already embedded Z-signals. 
% Default empty cell array.
%
% SIGMA determines the width of the time window in samples which is
% used for the estimation of multivariate transfer entropy at each time 
% instant. Larger windows give smaller errors, but less sensitivity
% to temporal changes in transfer entropy, and vice versa.
% Default 5.
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
% The cell arrays X and Y are treated as 1-dimensional arrays.
% Each signal is a real (m x n)-matrix that contains n samples of an
% m-dimensional signal. The signals need not have the same dimension,
% but the dimension must be the same among the trials of a given signal.
% If the number of samples varies with each signal, the function uses 
% the minimum sample count among the signals. The number of trials q 
% must coincide with all the signals in X, Y, and Z.

function I = transfer_entropy(X, W, Y, Z, sigma, k, threads)

% The limit for the dimension is arbitrary, but
% protects for the case when the user accidentally
% passes the transpose of the intended data.
maxDimension = 32;

if nargin < 3
    error('Not enough input arguments.');
end

if nargin > 7
    error('Too many input arguments.');
end

if nargout > 1
    error('Too many output arguments.');
end

if nargin < 4
    Z = cell(0);
end

if nargin < 5
    sigma = 5
end

if nargin < 6
    k = 1;
end

if nargin < 7
    threads = 1;
end

if ~iscell(X)
    error('X must be a cell array.');
end

if ~iscell(W)
    error('W must be a cell array.');
end

if ~iscell(Y)
    error('Y must be a cell array.');
end

if ~iscell(Z)
    error('Z must be a cell array.');
end

xSignals = prod(size(X));

if xSignals < 1
    error('X must contain at least one trial of at least one signal.');
end

wSignals = prod(size(W));

if wSignals < 1
    error('W must contain at least one trial of at least one signal.');
end

ySignals = prod(size(Y));

if ySignals < 1
    error('Y must contain at least one trial of at least one signal.');
end

for i = 1 : xSignals
    if ~isa(X{i}, 'double')
        error('Some signal of X is not of type double.');
    end

    if size(X{i}, 1) > maxDimension
        error(['Some signal of X has dimension greater than ', ...
            int2str(maxDimension), '.']);
    end
end

for i = 1 : wSignals
    if ~isa(W{i}, 'double')
        error('Some signal of W is not of type double.');
    end

    if size(W{i}, 1) > maxDimension
        error(['Some signal of W has dimension greater than ', ...
            int2str(maxDimension), '.']);
    end
end

for i = 1 : ySignals
    if ~isa(Y{i}, 'double')
        error('Some signal of Y is not of type double.');
    end

    if size(Y{i}, 1) > maxDimension
        error(['Some signal of Y has dimension greater than ', ...
            int2str(maxDimension), '.']);
    end
end

zSignals = prod(size(Z));

for i = 1 : zSignals
    if ~isa(Z{i}, 'double')
        error('Some signal of Z is not of type double.');
    end

    if size(Z{i}, 1) > maxDimension
        error(['Some signal of Z has dimension greater than ', ...
            int2str(maxDimension), '.']);
    end
end

if xSignals ~= wSignals || xSignals ~= ySignals || ...
        xSignals ~= wSignals || xSignals ~= size(zSignals, 2) 
    error('The number of trials does not match among X, W, Y, and Z.');
end

if size(k, 1) ~= 1 || ...
   size(k, 2) ~= 1
    error('K must be a scalar integer.');
end

if k < 1
    error('K must be at least 1.');
end

if size(sigma, 1) ~= 1 || ...
   size(sigma, 2) ~= 1
    error('SIGMA must be a scalar.');
end

if sigma < 0
    error('SIGMA must be a non-negative integer.');
end

if size(threads, 1) ~= 1 || ...
   size(threads, 2) ~= 1
    error('THREADS must be a scalar integer.');
end

if threads < 1
    error('THREADS must be at least 1.');
end

I = timTransferEntropy(X, W, Y, Z, sigma, k, threads);
