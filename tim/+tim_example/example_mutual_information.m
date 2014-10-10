% Description: Visualization of bias in mutual_information

clear all;
close all;

% This example demonstrates numerically that the 
% mutual_information estimator is biased. 

for i = 2 : 2
    d = 2^i;
    R = pastelmath.random_orthogonal(d, 'orientation', 1);
    D = rand(d, 1);
    M = R * diag(D);
    cov = M * M';
    
    for j = 0 : 3
        k = 2^j;
        figure;
        tim_example.draw_mutual_information('k', k, 'cov', cov);
    end
end
