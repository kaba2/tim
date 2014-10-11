% DIFFERENTIAL_ENTROPY_SINE
% Differential entropy of Y = sin(2 * pi * X), where X ~ Uniform(0, 1).
%
% H = differential_entropy_sine()
%
% The cumulative distribution function of Y is 
% F : [-1, 1] --> [0, 1] such that
%
%     F(y) = [asin(y) + (pi / 2)] / pi
%
% The probability density function of Y is
% f : [-1, 1] --> R such that
%
%     f(y) = 1 / (pi sqrt(1 - y^2)).

function H = differential_entropy_sine()

H = log2(pi) - 1;
