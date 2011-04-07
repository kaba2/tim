% DIFFERENTIAL_ENTROPY_KL_T
% A temporal differential entropy estimate from samples
% using Kozachenko-Leonenko nearest neighbor estimator.
%
% H = differential_entropy_kl_t(
%     S, timeWindowRadius, k, filter)
%
% where
%
% S is a signal set.
%
% Type 'help tim' for more documentation.

% Description: Temporal differential entropy estimation
% Detail: Kozachenko-Leonenko nearest neighbor estimator
% Documentation: differential_entropy_kl.txt

function H = differential_entropy_kl_t(...
    S, timeWindowRadius, k, filter)

concept_check(nargin, 'inputs', 2 : 4);
concept_check(nargout, 'outputs', 0 : 1);

if nargin < 3
    k = 1;
end

if nargin < 4
    filter = [1];
end

if isnumeric(S)
    S = {S};
end

concept_check(S, 'signalSet');
concept_check(timeWindowRadius, 'timeWindowRadius');
concept_check(k, 'k');
concept_check(filter, 'filter');

H = tim_matlab('differential_entropy_kl_t', ...
    S, timeWindowRadius, k, filter);
