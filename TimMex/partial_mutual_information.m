% PARTIAL_MUTUAL_INFORMATION
% A partial mutual information estimate I(X, Y | Z) from samples.
%
% I = partial_mutual_information(X, Y, Z, k)
%
% where
%
% X is a real (m x n)-matrix that contains n samples of an
% m-dimensional signal.
%
% Y is a real (p x n)-matrix that contains n samples of a
% p-dimensional signal.
%
% Z is a real (q x n)-matrix that contains n samples of a
% q-dimensional signal.
%
% K determines which k:th nearest neighbor the algorithm
% uses for estimation. Default 1.

function I = partial_mutual_information(X, Y, Z, k)

% The limit for the dimension is arbitrary, but
% protects for the case when the user accidentally
% passes the transpose of the intended data.
maxDimension = 32;

if nargin < 3
    error('Not enough input arguments.');
end

if nargin > 4
    error('Too many input arguments.');
end

if nargout > 1
    error('Too many output arguments.');
end

if nargin < 4
    k = 1;
end

if ~isfloat(X)
    error('X must be real and non-complex.');
end

if ~isfloat(Y)
    error('Y must be real and non-complex.');
end

if ~isfloat(Z)
    error('Z must be real and non-complex.');
end

% The limit for the dimension is arbitrary, but
% protects for the case when the user accidentally
% passes the transpose of the intended data.
if size(X, 2) > 32
    error('The width of X must be less than 32.');
end

if size(Y, 2) > 32
    error('The width of Y must be less than 32.');
end

if size(Z, 2) > 32
    error('The width of Z must be less than 32.');
end

if size(k, 1) ~= 1 || ...
   size(k, 2) ~= 1
    error('K must be a scalar integer.');
end

if k < 1
    error('K must be at least 1.');
end

I = timPartialMutualInformation(X, Y, Z, k);
