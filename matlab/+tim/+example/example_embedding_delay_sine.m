% Description: Uniform sampling with maximal information
% DocumentationOf: embedding_delay.m

clear all;
close all;

n = 256;

xMin = 0;
xMax = 10 * (2 * pi);

xSet = linspace(xMin, xMax, n);
ySet = sin(xSet);

maxLag = 256;

figure;
plot(xSet, ySet);

figure;
tim.auto_mi_plot([xSet; ySet], 'lagSet', 0 : maxLag);

delay = tim.embedding_delay([xSet; ySet], 'maxLag', maxLag);

disp(['A good embedding delay for the sine is ', ...
    int2str(delay), '.']);


