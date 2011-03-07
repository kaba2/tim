% DIFFERENTIAL_ENTROPY_SP
% A differential entropy estimate from samples
% using Stowell-Plumbley recursive partition estimator.
%
% H = differential_entropy_sp(S)
%
% where
%
% S is a signal set.
%
% Type 'help tim' for documentation.

% Description: Differential entropy estimation
% Detail: Stowell-Plumbley recursive partition estimator
% Documentation: differential_entropy_sp.txt

function H = differential_entropy_sp(S)

check(nargin, 'inputs', 1);
check(nargout, 'outputs', 0 : 1);

if isnumeric(S)
    S = {S};
end

check(S, 'signalSet');

if (size(S{1}, 1) > 3)
	warning('This estimator has bad accuracy for dimensions > 3!');
end

H = tim_matlab('differential_entropy_sp', S);
