clear all;
close all;

[aSet, bSet] = tim.coupled_lorenz_systems();

figure;
tim.auto_mi_plot(aSet);

figure;
tim.auto_mi_plot(tim.lorenz_system());

aFutureSet = tim.delay_embed_future(aSet, 1);

lagSet = -10 : 10;
abSet = tim.transfer_entropy(aSet, bSet, aFutureSet, ...
	'xLag', lagSet, 'wLag', lagSet);

figure;
plot(lagSet, abSet);
xlabel('Lag (samples)')
ylabel('TE(X --> Y)')

lagSet = -10 : 10;
abSet = tim.transfer_entropy(aSet, bSet, aFutureSet, ...
	'xLag', lagSet, 'wLag', lagSet);


