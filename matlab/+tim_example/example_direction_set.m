clear all;
close all;

pointSet = tim.lorenz_system();

figure;
tim.spatial_plot(pointSet);

directionSet = tim.compute_direction_set(pointSet);
figure;
tim.spatial_plot(directionSet);

d = size(pointSet, 1);

kdTree = pastelmatlab.PointKdTree(d);
idSet = kdTree.insert(directionSet);
kdTree.refine();
direction = pointSet(:, 3001) - pointSet(:, 3000);
direction = direction / norm(direction);
[indexSet, ignore] = kdTree.search_nearest(direction(:), 'maxDistanceSet', 2 - 2 * 0.98, 'kNearest', 128);

indexSet = sort(indexSet);

figure;
tim.spatial_plot(pointSet, 'subset', indexSet);

clear kdTree;


