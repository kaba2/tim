% Non-linear time series analysis
% Exercise 3.3

clear all;
close all;

pointSet = tim.lorenz();

sectionSet = tim.poincare_section(pointSet, 'direction', -1);
timtest.visualize_embedding(sectionSet);
