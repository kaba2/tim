function predictedSet = predict_all(signal)

import([tim_package, '.*']);

[d, n] = size(signal);

predictor = Predictor(d);
predictor.setKNearest(1);

predictedSet = zeros(d, n);
for i = 1 : (n - 1)
	predictor.insert(signal(:, i));
    prediction = predictor.predict(signal(:, i));
    predictedSet(:, i) = mean(prediction, 2);

    %[ignore, index] = min(prediction(end, :));
    %predictedSet(:, i) = prediction(:, index);
end

clear predictor;
