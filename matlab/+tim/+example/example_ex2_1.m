% Non-linear time series analysis
% Exercise 2.1

clear all;
close all;

n = 1024;

dataSet = randn(1, n);
tim.example.visualize(dataSet);

% Take a moving-average. This should cut down
% on the high frequencies.
radius = 7;
filteredSet = zeros(1, n);
for i = 1 : n
    filteredSet(i) = mean(dataSet(1, ...
        max(i - radius, 1) : min(i + radius, n)));
end

tim.example.visualize(filteredSet);