function h = visualize(dataSet)

n = size(dataSet, 2);

fftSet = fft(dataSet);
spectrumSet = abs(fftSet);

h = figure;

% Draw the data.
subplot(2, 1, 1);
plot(dataSet);
title('Data');
xlabel('Time (s)');
ylabel('Amplitude');
axis([0, n, -1.5, 1.5]);

% Draw the power spectrum.
subplot(2, 1, 2);
plot(linspace(-0.5 * n, 0.5 * n, n), fftshift(spectrumSet));
title('Power spectrum');
xlabel('Frequency (Hz)');
ylabel('Power');
