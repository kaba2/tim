% EXAMPLE_DIFFERENTIAL_ENTROPY_SINE
%
% Discontinuities in the probability density function break the 
% assumptions of the SP and KL estimators. As a result, they 
% perform badly. This example demonstrates this aspect by computing
% the differential entropy of Y = sin(2 * pi * X), where 
% X ~ Uniform(0, 1).

clear all;
close all;

% Number of samples to generate.
n = 200000;

% Let X ~ Uniform(0, 1).
xSet = rand(1, n);
% Let Y = sin(2 * pi * X).
ySet = sin(2 * pi * xSet);

% The probability density of Y is given by 
% f(x) = 1 / (pi * sqrt(1 - x^2)).
f = @(x) 1 ./ (pi * sqrt(1 - x.^2));

edgeSet = linspace(-1, 1, 1000);

figure;
hold on;
% Draw empirical probability density
[hPlot, pdfSet] = tim.density_plot(ySet, edgeSet);
% Draw analytic probability density
h = ezplot(f);
set(h, 'color', 'red');
axis([-1, 1, 0, 2]);
xlabel('Y = sin(2 * pi * X), X ~ U(0, 1)');
title('Probability density function of Y');
hold off;

disp(' ');
disp('Differential entropy of sin(2 * pi * X), when X ~ Uniform(0, 1)');
disp([int2str(n), ' samples']);

analyticDe = tim.differential_entropy_sine();
disp(' ');
disp(['differential_entropy_sine: ', num2str(analyticDe)]);

klDe = tim.differential_entropy_kl(ySet);
disp(['differential_entropy_kl (k = 1): ', num2str(klDe)]);

klDe = tim.differential_entropy_kl(ySet, 'k', 100);
disp(['differential_entropy_kl (k = 100): ', num2str(klDe)]);

spDe = tim.differential_entropy_sp(ySet);
disp(['differential_entropy_sp: ', num2str(spDe)]);

pmfSet = pdfSet / sum(pdfSet);
binningDe = -sum(pdfSet .* log2(pdfSet)) * (2 / (numel(edgeSet) - 1));
disp(['Binning estimator: ', num2str(binningDe)]);

