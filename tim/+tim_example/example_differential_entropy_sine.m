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
hold off;
xlabel('Y = sin(2 * pi * X), X ~ U(0, 1)');

disp(' ');
disp('Differential entropy of sin(2 * pi * X), when X ~ Uniform(0, 1)');
disp([int2str(n), ' samples']);

% The analytic differential entropy of Y is
% log2(pi) - 1.
analyticDe = log2(pi) - 1;
disp(' ');
disp(['Analytic: ', num2str(analyticDe)]);

klDe = tim.differential_entropy_kl(ySet);
disp(['Kozachenko-Leonenko estimator: ', num2str(klDe)]);

spDe = tim.differential_entropy_sp(ySet);
disp(['Stowell-Plumbley estimator: ', num2str(spDe)]);

pmfSet = pdfSet / sum(pdfSet);
binningDe = -sum(pdfSet .* log2(pdfSet)) * (2 / (numel(edgeSet) - 1));
disp(['Binning estimator: ', num2str(binningDe)]);

disp(' ');
disp(['Lesson: discontinuities in the probability density ', ...
    'break the assumptions of the SP and KL estimators, and ', ...
    'they perform badly.']);
