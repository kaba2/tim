% EXAMPLE_DIFFERENTIAL_ENTROPY_KL
% Effect of k and dimension on differential_entropy_kl
%
% This example demonstrates numerically that the 
% differential_entropy_kl estimator is biased. Moreover,
%
% * as dimension increases, the bias increases,
% * as k increases, the bias increases, and the variance decreases.
%
% The test-distribution is the multi-variate standard normal
% distribution. The differential entropy of this distribution is given by 
% tim.differential_entropy_normal().
%
% The k and d are varied over [1, 2, 4]. This generates 9 combinations.
% For every (k, d)-pair a figure is drawn of the empirical distribution 
% of differential_entropy_kl(), together with sample mean and the
% analytical solution.

% Description: Effect of k and dimension on differential_entropy_kl

clear all;
close all;

for i = 0 : 2
    for j = 0 : 2
        p = 4 * i + j + 1;
        d = 2^i;
        k = 2^j;
        
        figure;
        tim.example.draw_differential_entropy_kl('k', k, 'd', d);
    end
end
