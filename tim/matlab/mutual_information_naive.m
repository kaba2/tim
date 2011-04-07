% MUTUAL_INFORMATION_NAIVE
% A mutual information estimate from samples.
%
% I = mutual_information_naive(S)
% I = mutual_information_naive(S, 'key', value, ...)
%
% where
%
% S is a real (m x n)-matrix that contains n samples of an
% m-dimensional signal.
%
% Additional (optional) input arguments can be provided as 'key', value
% pairs. Accepted ('key', VALUE) pairs are:
%
% ('bins', BINS) is an integer that determines the number of bins to use
% for 1d distribution estimation. Default: 100.
%
%
% Type 'help tim' for more documentation.

% Description: Naive mutual information estimation
% Documentation: mutual_information.txt

function I = mutual_information_naive(S, varargin)

pkgname = regexpi(mfilename('fullpath'), ['+(TIM_.*)' filesep mfilename], ...
    'tokens', 'once');
import([pkgname{1} '.tim_matlab']);
import([pkgname{1} '.concept_check']);
import([pkgname{1} '.process_options']);

keySet = {'bins'};

% Default values for the optional input arguments
bins = 100;

eval(process_options(keySet, varargin));

concept_check(nargin, 'inputs', 1 : 2);
concept_check(nargout, 'outputs', 0 : 1);


if ~isa(S, 'double')
    error('S must be of type double.');
end

if size(bins, 1) ~= 1 || ...
   size(bins, 2) ~= 1
    error('BINS must be a scalar integer.');
end

I = tim_matlab('mutual_information_naive',  S, bins);
