% TSALLIS_ENTROPY_LPS
% A Tsallis entropy estimate from samples
% using Leonenko-Pronzato-Savani nearest neighbor estimator.
%
% H = tsallis_entropy_lps(S)
% H = tsallis_entropy_lps(S, 'key', value, ...)
%
% where
%
% S is a signal set.
%
% H is the estimated Tsallis entropy.
%
% Optional input arguments in 'key'-value pairs:
%
% Q ('q') is the power in the definition Renyi entropy.
% If Q = 1, differential_entropy_kl() is used to
% compute the result instead. 
% If Q < 1, there are huge errors in the estimation.
% Default: 2.
%
% KSUGGESTION ('kSuggestion') is a suggestion for the k:th nearest
% neighbor that should be used for estimation. The k can't
% be freely set because the estimation algorithm is only defined
% for k > q - 1. Value zero means an accurate (q-dependent) default 
% is used. For accurate results one should choose 
% kSuggestion >= 2 * ceil(q) - 1. Default: 0.
%
% Type 'help tim' for more documentation.

% Description: Tsallis entropy estimation
% Detail: Leonenko-Pronzato-Savani nearest neighbor estimator
% Documentation: tsallis_entropy_lps.txt

function H = tsallis_entropy_lps(S, varargin)

import([tim_package, '.*']);

concept_check(nargin, 'inputs', 1);
concept_check(nargout, 'outputs', 0 : 1);

% Optional input arguments
q = 2;
kSuggestion = 0;
eval(process_options({'q', 'kSuggestion'}, varargin));

if isnumeric(S)
    S = {S};
end

pastelmatlab.concept_check(...
	S, tim_package('signal_set'), ...
	q, 'real', ...
	q, 'positive', ...
	kSuggestion, 'integer', ...
	kSuggestion, 'non_negative');

H = tim_matlab('tsallis_entropy_lps', ...
	S, q, kSuggestion);
