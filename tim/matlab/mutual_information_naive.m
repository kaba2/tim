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
%
% Type 'help tim' for more documentation.

% Description: Naive mutual information estimation
% Documentation: tim_matlab_matlab.txt

function I = mutual_information_naive(S, bins)

check(nargin, 'inputs', 1 : 2);
check(nargout, 'outputs', 0 : 1);

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
