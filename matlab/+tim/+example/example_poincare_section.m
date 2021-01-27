% Description: Testing for poincare sections
% DocumentationOf: poincare_section.m

% Non-linear time series analysis
% Exercise 3.3

clear all;
close all;

pointSet = tim.lorenz_system();

sectionSet = tim.poincare_section(pointSet, 'direction', -1);
tim.test.visualize_embedding(sectionSet);
