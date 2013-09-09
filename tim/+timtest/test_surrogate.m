clear all;
close all;

n = 201;
filterOrder = 21;

dataSet = randn(1, n);
%windowSet = blackman(n)';
%dataSet = windowSet .* dataSet;

filterSet = fir1(filterOrder, 0.5);
filteredSet = filter(filterSet, 1, dataSet);

%surrogateSet = tim.surrogate(filteredSet, ...
%		'frequencyRange', [0, 1], ...
%		'algorithm', 'preserve_correlations');

surrogateSet = tim.surrogate(filteredSet, ...
		'algorithm', 'preserve_dynamics');

mean(abs(fft(filteredSet)) - abs(fft(surrogateSet)))

figure;
hold on;
area(filteredSet);
plot(surrogateSet, 'r');
hold off;

figure;
hold all;
pwelch(filteredSet);
pwelch(surrogateSet);
axis([0, 1, -60, 0]);
hold off;
