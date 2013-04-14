clear all;
close all;

pointSet = tim.lorenz();

tim.recurrence_plot(pointSet, 'deviation', 8);
