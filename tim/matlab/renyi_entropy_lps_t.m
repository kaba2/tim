% RENYI_ENTROPY_LPS_T
% A temporal Renyi entropy estimate from samples
% using Leonenko-Pronzato-Savani nearest neighbor estimator.
%
% H = renyi_entropy_lps_t(
%     S, timeWindowRadius, q, kSuggestion, filter)
%
% where
%
% S is a signal set.
%
% Q is the power in the definition Renyi entropy.
% If Q = 1, differential_entropy_kl_t() is used to
% compute the result instead. 
% If Q < 1, there are huge errors in the estimation.
% Default 2.
%
% KSUGGESTION is a suggestion for the k:th nearest neighbor 
% that should be used for estimation. The k can't
% be freely set because the estimation algorithm is only defined
% for k > q - 1. Value zero means an accurate (q-dependent) default 
% is used. For accurate results one should choose 
% kSuggestion >= 2 * ceil(q) - 1. Default 0.
%
% Type 'help tim' for more documentation.

% Description: Temporal Renyi entropy estimation
% Detail: Leonenko-Pronzato-Savani nearest neighbor estimator
% Documentation: renyi_entropy_lps.txt

function H = renyi_entropy_lps_t(...
    S, timeWindowRadius, q, kSuggestion, filter)

concept_check(nargin, 'inputs', 2 : 5);
concept_check(nargout, 'outputs', 0 : 1);

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

concept_check(S, 'signalSet');
concept_check(timeWindowRadius, 'timeWindowRadius');

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

concept_check(filter, 'filter');

H = tim_matlab('renyi_entropy_lps_t', ...
    S, timeWindowRadius, q, kSuggestion, filter);
