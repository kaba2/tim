% Description: Testing for recurrence plots
% DocumentationOf: recurrence_plot.m

clear all;
close all;

pointSet = tim.lorenz_system();

tim.recurrence_plot(pointSet, 'deviation', 8);
