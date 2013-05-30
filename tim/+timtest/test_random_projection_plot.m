clear all;
close all;

pointSet = tim.lorenz_system();

tim.random_projection_plot(pointSet, ...
    'spatial_plot', {'draw_points', false});

