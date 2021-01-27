% EXAMPLE_MUTUAL_INFORMATION_SINES
%
% Let
%
%     f(x) = sin(x)
%     g(x) = sin(x + alpha)
%
% Suppose X is a random real number uniformly distributed on [0, 2 * pi].
% Are f(X) and g(X) correlated? To answer this question, we will compute 
% the mutual information between f(X) and g(X) while varying alpha in 
% the interval [-pi, pi].

% Description: Mutual information between sinusoids
% DocumentationOf: mutual_information.m

clear all;
close all;

m = 127;
n = 2500;

xSet = rand(1, n) * (2 * pi);
alphaSet = linspace(-pi, pi, m);

miSet = zeros(1, m);
for i = 1 : m
    alpha = alphaSet(i);

    aSet = sin(xSet);
    bSet = sin(xSet + alpha);
    
    miSet(i) = tim.mutual_information(aSet, bSet);
end

figure;
hold on;
plot(alphaSet, miSet);
title('Mutual information between sinusoids');
ylabel('Mutual information (bits)');
xlabel('Offset angle (radians)');
axis([-pi, pi, 0, ceil(max(miSet))]);

