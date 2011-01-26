% TSALLIS_ENTROPY_LPS
% A Tsallis entropy estimate from samples
% using Leonenko-Pronzato-Savani nearest neighbor estimator.
%
% H = tsallis_entropy_lps(S, q, kSuggestion)
%
% where
%
% S is a signal set.
%
% Q is the power in the definition Tsallis entropy.
% In case Q = 1, differential_entropy_kl() is used to
% compute the result instead. Default 2.
%
% KSUGGESTION is a suggestion for the k:th nearest neighbor 
% that should be used for estimation. The k can't
% be freely set because the estimation algorithm is only defined
% for k > q - 1. Value zero means an accurate (q-dependent) default 
% is used. For accurate results one should choose 
% kSuggestion >= 2 * ceil(q) - 1. Default 0.
%
% Type 'help tim' for more documentation.

% Description: Tsallis entropy estimation
% Detail: Leonenko-Pronzato-Savani nearest neighbor estimator
% Documentation: tim_matlab_matlab.txt

function H = tsallis_entropy_lps(S, q, kSuggestion)

check(nargin, 'inputs', 1 : 3);
check(nargout, 'outputs', 0 : 1);

if nargin < 2
	q = 2;
end

if nargin < 3
    kSuggestion = 0;
end

if isnumeric(S)
    S = {S};
end

check(S, 'signalSet');

if size(q, 1) ~= 1 || ...
   size(q, 2) ~= 1
    error('Q must be a scalar.');
end

if size(kSuggestion, 1) ~= 1 || ...
   size(kSuggestion, 2) ~= 1
    error('KSUGGESTION must be a scalar integer.');
end

if q <= 0
	error('Q must be positive');
end

if kSuggestion < 0
    error('KSUGGESTION must be non-negative.');
end

H = tim_matlab('tsallis_entropy_lps', ...
	S, q, kSuggestion);
