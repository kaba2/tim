% DIFFERENTIAL_ENTROPY A differential entropy estimate from samples.
% H = differential_entropy(S, epsilon, k)
% where
% S is a real (m x n)-matrix that contains n samples of an
% m-dimensional signal.
% EPSILON is the maximum relative error in distance that
% nearest neighbor searching is allowed to result in.
% Higher tolerances result in enhanced performance, but
% increases errors in the estimate. Default 0.
% K determines which k:th nearest neighbor the algorithm
% uses for estimation. Default 1.

function H = differential_entropy(S, epsilon, k)

if nargin < 1
    error('Not enough input arguments.');
end

if nargin > 3
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

if ~isfloat(S)
    error('S must be real and non-complex.');
end

% The limit for the dimension is arbitrary, but
% protects for the case when the user accidentally
% passes the transpose of the intended data.
if size(S, 1) > 32
    error('The height of S must be less than 32.');
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

H = timDifferentialEntropy(S, epsilon, k);
