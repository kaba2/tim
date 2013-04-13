clear all;
close all;

pointSet = tim.lorenz();

figure;
plot3(pointSet(1, :), pointSet(2, :), pointSet(3, :));
title('Lorenz system as generated')

figure;
plot(pointSet(1, :));
title('Lorenz system X-coordinates')

maxLag = 17;
lags = 3;

for lag = int32(linspace(1, maxLag, lags))
    % 2D delay-embedding.
    embeddedSet = tim.delay_embed(pointSet(1, :), 2, lag);
    
    lagText = ['(lag ', num2str(lag), ')'];
    
    figure;
    plot(embeddedSet(1, :), embeddedSet(2, :));
    title(['2D delay-embedding of the Lorenz X-coordinates ', ...
        lagText]);

    % 3D delay-embedding.
    embeddedSet = tim.delay_embed(pointSet(1, :), 3, lag);

    figure;
    plot3(embeddedSet(1, :), embeddedSet(2, :), embeddedSet(3, :));
    title(['3D delay-embedding of the Lorenz X-coordinates ', ...
        lagText]);
end

% A good embedding lag is given by the first minimum of the
% auto mutual information. If the signal is sampled densely enough, 
% temporally close samples are related simply because of continuity. 
% Increasing the lag decreases the continuity-caused dependency
% between delayed samples. On the other hand, at some point the
% dependency starts to rise again, since a chaotic system usually
% has cyclic-like properties. Overall, the dependency will lower,
% since temporally distant samples have less and less dependency to
% each other, because of the chaotic pseudo-randomness.

lagSet = floor(linspace(1, 200, 50));
miSet = tim.mutual_information(...
    pointSet(1, :), pointSet(1, :), ...
    'yLag', lagSet);

figure;
plot(int32(lagSet), miSet);
title('Auto mutual information of the Lorenz X-coordinates');
xlabel('Lag');
ylabel('Mutual information');
