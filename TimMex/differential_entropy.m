% DIFFERENTIAL_ENTROPY
% A differential entropy estimate from samples.
%
% H = differential_entropy(S, epsilon, k, threads)
%
% where
%
% S is a real (m x n)-matrix that contains n samples of an
% m-dimensional signal.
%
% EPSILON is the maximum relative error in distance that
% nearest neighbor searching is allowed to result in.
% Higher tolerances result in enhanced performance, but
% increases errors in the estimate. Default 0.
%
% K determines which k:th nearest neighbor the algorithm
% uses for estimation. Default 1.
%
% THREADS determines the number of threads to use for parallelization.
% To fully take advantage of multiple cores in your machine, set this
% to the number of cores in your machine. Note however that this makes 
% your computer unresponsive to other tasks. When you need responsiveness, 
% spare one core for other work. Default 1 (no parallelization).

function H = differential_entropy(S, epsilon, k, threads)

% The limit for the dimension is arbitrary, but
% protects for the case when the user accidentally
% passes the transpose of the intended data.
maxDimension = 32;

if nargin < 1
    error('Not enough input arguments.');
end

if nargin > 4
    error('Too many input arguments.');
end

if nargout > 1
    error('Too many output arguments.');
end

if nargin < 2
    epsilon = 0;
end

if nargin < 3
    k = 1;
end

if nargin < 4
    threads = 1;
end

if ~isa(S, 'double')
    error('S must be of type double.');
end

if size(S, 1) > maxDimension
    error(['Some signal of S has dimension greater than ', ...
        int2str(maxDimension), '.']);
end

if size(epsilon, 1) ~= 1 || ...
   size(epsilon, 2) ~= 1
    error('EPSILON must be a scalar.');
end

if size(k, 1) ~= 1 || ...
   size(k, 2) ~= 1
    error('K must be a scalar integer.');
end

if epsilon < 0
    error('EPSILON must be positive'.');
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

H = timDifferentialEntropy(S, epsilon, k, threads);
