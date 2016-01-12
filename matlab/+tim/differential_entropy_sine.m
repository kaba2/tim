% DIFFERENTIAL_ENTROPY_SINE
% Differential entropy of Y = sin(2 * pi * X + beta), X ~ Uniform(0, 1).
%
% H = differential_entropy_sine()
%
% Here BETA is a real number. The differential entropy of Y is 
% not dependent of BETA.
%
% The cumulative distribution function of Y is 
% F : [-1, 1] --> [0, 1] such that
%
%     F(y) = [max(asin(y) + (pi / 2) - beta, 0) + 
%             max(asin(y) - (pi / 2) + beta, 0)] / pi
%
% Let 
%
%     a = sin(beta - pi / 2)
%     b = sin(pi / 2 - beta)
%
% The probability density function of Y is
% f : [-1, 1] --> R such that 
%
%     f(y) = 2 / (pi sqrt(1 - y^2)),    a <= y and b <= y,
%            1 / (pi sqrt(1 - y^2)),    a <= y xor b <= y,
%            0,                         a > y and b > y,
%
% The differential entropy of Y can then be solved as
%
%     H(Y) = log2(pi) - 1.

% Description: Differential entropy of Y = sin(2 * pi * X + beta), X ~ Uniform(0, 1).
% Documentation: differential_entropy.txt

function H = differential_entropy_sine()

import([tim_package, '.*']);

concept_check(nargin, 'inputs', 0);
concept_check(nargout, 'outputs', 0 : 1);

H = log2(pi) - 1;
