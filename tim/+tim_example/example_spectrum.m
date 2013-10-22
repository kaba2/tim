clear all;
close all;

n = 100;
angleWidth = 2 * pi;

angleSet = linspace(0, angleWidth, n);
wave = @sin;
%wave = @sawtooth;

dataSet = wave(10 * angleSet) + 0.5 * wave(5 * angleSet);
tim_example.visualize(dataSet);

% Zero phases.
modifiedDataSet = ifft(abs(fft(dataSet)));
tim_example.visualize(modifiedDataSet);

