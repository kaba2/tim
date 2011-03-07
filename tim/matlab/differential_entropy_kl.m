% DIFFERENTIAL_ENTROPY_KL
% A differential entropy estimate from samples
% using Kozachenko-Leonenko nearest neighbor estimator.
%
% H = differential_entropy_kl(S, k)
%
% where
%
% S is a signal set.
%
% Type 'help tim' for more documentation.

% Description: Differential entropy estimation
% Detail: Kozachenko-Leonenko nearest neighbor estimator
% Documentation: differential_entropy_kl.txt

function H = differential_entropy_kl(S, k)

check(nargin, 'inputs', 1 : 2);
check(nargout, 'outputs', 0 : 1);

if nargin < 2
    k = 1;
end

if isnumeric(S)
    S = {S};
end

check(S, 'signalSet');
check(k, 'k');

H = tim_matlab('differential_entropy_kl', ...
	S, k);
