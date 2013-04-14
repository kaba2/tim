clear all;
close all;

pointSet = tim.lorenz('n', 500, 'tMax', 50);

timtest.visualize_embedding(pointSet, 'tDelta', 1);

distanceMatrix = sqrt(tim.distance_matrix(pointSet));

figure;
h = contourf(distanceMatrix);
shading flat;
colormap gray;
%set(h, 'ShowText', 'on')
colorbar();
