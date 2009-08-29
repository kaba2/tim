% PARTIAL_MUTUAL_INFORMATION 
% A partial mutual information estimate from samples.
%
% I = partial_mutual_information(X, Y, Z, yLag, zLag, k, threads)
%
% where
%
% X, Y, and Z are arbitrary-dimensional cell-arrays whose 
% linearizations contain q trials of signal x, y, and z, 
% respectively.
%
% YLAG and ZLAG are the lags in samples applied to signals
% y and z, respectively.
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
% m-dimensional signal. The dimensions of X, Y, and Z need not coincide.
% However, the number of trials has to coincide.
% If the number of samples varies with trials, the function uses 
% the minimum sample count among the trials of X and Y.

% Description: Partial mutual information estimation
% Documentation: tim_matlab.txt

function I = partial_mutual_information(...
    X, Y, Z, yLag, zLag, k, threads)

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
    yLag = 0;
end

if nargin < 5
    zLag = 0;
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

if ~iscell(Y)
    error('Y must be a cell array.');
end

if ~iscell(Z)
    error('Z must be a cell array.');
end

signals = prod(size(X));

if prod(size(Y)) ~= signals || prod(size(Z)) ~= signals
    error('X, Y, and Z must contain the same number of signals.');
end

for i = 1 : signals
    if ~isa(X{i}, 'double')
        error('Some signal of X is not of type double.');
    end

    if size(X{i}, 1) > maxDimension
        error(['Some signal of X has dimension greater than ', ...
            int2str(maxDimension), '.']);
    end
end

for i = 1 : signals
    if ~isa(Y{i}, 'double')
        error('Some signal of Y is not of type double.');
    end

    if size(Y{i}, 1) > maxDimension
        error(['Some signal of Y has dimension greater than ', ...
            int2str(maxDimension), '.']);
    end
end

for i = 1 : signals
    if ~isa(Z{i}, 'double')
        error('Some signal of Z is not of type double.');
    end

    if size(Z{i}, 1) > maxDimension
        error(['Some signal of Z has dimension greater than ', ...
            int2str(maxDimension), '.']);
    end
end

if size(yLag, 1) ~= 1 || ...
   size(yLag, 2) ~= 1
    error('YLAG must be a scalar integer.');
end

if yLag < 0
    error('YLAG must be non-negative.');
end

if size(zLag, 1) ~= 1 || ...
   size(zLag, 2) ~= 1
    error('ZLAG must be a scalar integer.');
end

if zLag < 0
    error('ZLAG must be non-negative.');
end

if size(k, 1) ~= 1 || ...
   size(k, 2) ~= 1
    error('K must be a scalar integer.');
end

if k < 1
    error('K must be at least 1.');
end

if size(threads, 1) ~= 1 || ...
   size(threads, 2) ~= 1
    error('THREADS must be a scalar integer.');
end

if threads < 1
    error('THREADS must be at least 1.');
end

%I = timPartialMutualInformation(X, Y, Z, yLag, zLag, k, threads);

