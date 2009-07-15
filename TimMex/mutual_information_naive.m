% MUTUAL_INFORMATION_NAIVE
% A mutual information estimate from samples.
%
% I = mutual_information_naive(S, bins)
%
% where
%
% S is a real (m x n)-matrix that contains n samples of an
% m-dimensional signal.
%
% BINS determines the number of bins to use for 1d
% distribution estimation. Default 100.

function I = mutual_information_naive(S, bins)

if nargin < 1
    error('Not enough input arguments.');
end

if nargin > 2
    error('Too many input arguments.');
end

if nargout > 1
    error('Too many output arguments.');
end

if nargin < 2
    bins = 100;
end

if ~isa(S, 'double')
    error('S must be of type double.');
end

if size(bins, 1) ~= 1 || ...
   size(bins, 2) ~= 1
    error('BINS must be a scalar integer.');
end

I = timMutualInformationNaive(S, bins);
