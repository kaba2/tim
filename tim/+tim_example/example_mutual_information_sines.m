% Description: Mutual information between sinusoids
% DocumentationOf: mutual_information.m

clear all;
close all;

m = 100;
n = 256;

xMin = 0;
xMax = 10 * (2 * pi);

xSet = linspace(xMin, xMax, n);
alphaSet = linspace(0, 2 * pi, m);

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
axis([0, 2 * pi, 0, ceil(max(miSet))]);

