% EXAMPLE_LINEAR_DEPENDENCY
%
% This example demonstrates that Pearson's linear correlation,
% Spearman's rank correlation, and Kendall's tau are all unable
% to detect simple non-linear dependencies between random variables.
% This contrasts with mutual information.
%
% Let
%
% X ~ Uniform(-3, 3)
% E ~ Normal(0, 1)
% Y = X + 0.2 E
% Z = X^2 + 0.2 E
%
% The problem is to answer the following:
% Are X and Y correlated?
% Are X and Z correlated?
%
% All of the dependency measures detect a strong dependency between
% X and Y. However, the three first measures fail to detect any
% dependency between X and Z.
%
% Reference
% ---------
%
% This example has been adapted from:
%
% Statistical Validation of Mutual Information Calculations:
% Comparison of Alternative Numerical Algorithms,
% C. J. Cellucci, A. M. Albano, P. E. Rapp,
% Physical Review E, Volume 71,
% 2005.

clear all;
close all;

n = 300;

xSet = -3 + rand(1, n) * 6;

figure;
scatter(xSet, xSet, '.');
xlabel('X');
ylabel('X');

ySet = xSet + 0.2 * randn(size(xSet));

figure;
scatter(xSet, ySet, '.');
xlabel('X');
ylabel('Y = X + 0.2 E');

zSet = xSet.^2 + 0.2 * randn(size(xSet));

figure;
scatter(xSet, zSet, '.');
xlabel('X');
ylabel('Y = X^2 + 0.2 E');

disp(' ');
disp('X ~ Uniform(-3, 3)');
disp('E ~ Normal(0, 1)');
disp('Y = X + 0.2 E');
disp('Z = X^2 + 0.2 E');

disp(' ');
disp('rho = Pearson linear correlation');
disp('rank = Spearman rank correlation');
disp('tau = Kendall rank correlation');
disp('MI = tim.mutual_information');

rhoXy = corr(xSet', ySet', 'type', 'Pearson');
rhoXz = corr(xSet', zSet', 'type', 'Pearson');

rankXy = corr(xSet', ySet', 'type', 'Spearman');
rankXz = corr(xSet', zSet', 'type', 'Spearman');

tauXy = corr(xSet', ySet', 'type', 'Kendall');
tauXz = corr(xSet', zSet', 'type', 'Kendall');

miXy = tim.mutual_information(xSet, ySet);
miXz = tim.mutual_information(xSet, zSet);

disp(' ');
disp(['rho(X, Y) = ', num2str(rhoXy)]);
disp(['rank(X, Y) = ', num2str(rankXy)]);
disp(['tau(X, Y) = ', num2str(tauXy)]);
disp(['MI(X, Y) = ', num2str(miXy)]);

disp(' ');
disp(['rho(X, Z) = ', num2str(rhoXz)]);
disp(['rank(X, Z) = ', num2str(rankXz)]);
disp(['tau(X, Z) = ', num2str(tauXz)]);
disp(['MI(X, Z) = ', num2str(miXz)]);