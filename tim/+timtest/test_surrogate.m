clear all;
close all;

n = 200;
filterOrder = 21;

dataSet = randn(1, n);
filterSet = fir1(filterOrder, 0.5);
filteredSet = filter(filterSet, 1, dataSet);
surrogateSet = tim.surrogate(filteredSet, ...
	'frequency_range', [0.4, 1], ...
	'algorithm', 'preserve_dynamics', ...
	'k', 1);

figure;
hold on;
area(filteredSet);
plot(surrogateSet, 'r');
hold off;

figure;
hold on;
area(fftshift(abs(fft(filteredSet))));
plot(fftshift(abs(fft(surrogateSet))), 'r');
hold off;
