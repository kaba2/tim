clear all;
close all;

n = 63;
filterOrder = 16;

dataSet = randn(1, n);
filterSet = fir1(filterOrder, 0.5);
filteredSet = filter(filterSet, 1, dataSet);
surrogateSet = tim.surrogate(filteredSet, ...
	'frequency_range', [0.4, 1]);

figure;
hold all;
plot(filteredSet);
plot(surrogateSet);
hold off;
