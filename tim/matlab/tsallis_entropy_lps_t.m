% TSALLIS_ENTROPY_LPS_T
% A temporal Tsallis entropy estimate from samples
% using Leonenko-Pronzato-Savani nearest neighbor estimator.
%
% H = tsallis_entropy_lps_t(
%     S, timeWindowRadius, q, kSuggestion, filter)
%
% where
%
% S is a signal set.
%
% Q is the power in the definition Tsallis entropy.
% In case Q = 1, differential_entropy_kl_t() is used to
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

% Description: Temporal Tsallis entropy estimation
% Detail: Leonenko-Pronzato-Savani nearest neighbor estimator
% Documentation: tsallis_entropy_lps.txt

function H = tsallis_entropy_lps_t(...
    S, timeWindowRadius, q, kSuggestion, filter)

check(nargin, 'inputs', 2 : 5);
check(nargout, 'outputs', 0 : 1);

if nargin < 3
    q = 2;
end

if nargin < 4
    kSuggestion = 0;
end

if nargin < 5
    filter = 1;
end

if isnumeric(S)
    S = {S};
end

check(S, 'signalSet');
check(timeWindowRadius, 'timeWindowRadius');

if size(q, 1) ~= 1 || ...
   size(q, 2) ~= 1
    error('Q must be a scalar.');
end

if q <= 0
	error('Q must be positive');
end

if size(kSuggestion, 1) ~= 1 || ...
   size(kSuggestion, 2) ~= 1
    error('KSUGGESTION must be a scalar integer.');
end

if kSuggestion < 0
    error('KSUGGESTION must be non-negative.');
end

check(filter, 'filter');

H = tim_matlab('tsallis_entropy_lps_t', ...
    S, timeWindowRadius, q, kSuggestion, filter);
