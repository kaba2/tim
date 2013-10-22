% Description: Visualization of bias in differential_entropy_kl

clear all;
close all;

% This example demonstrates numerically that the 
% differential_entropy_kl estimator is biased. 
% Further:
%
% 1) as d increases, the bias increases,
% 2) as k increases, the bias increases,
% 3) as k increases, the variance decreases.

for i = 0 : 2
    for j = 0 : 2
        p = 4 * i + j + 1;
        d = 2^i;
        k = 2^j;
        figure;
        example.draw_differential_entropy_kl('k', k, 'd', d);
    end
end
