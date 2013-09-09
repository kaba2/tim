% SURROGATE
% Generates random data resembling the input data.
%
% surrogateSet = surrogate(inputSet)
%
% where
%
% INPUTSET is a real (d x n)-matrix, containing n samples of
% a d-dimensional signal.
%
% Return values
% -------------
%
% SURROGATESET is a real (d x n)-matrix, containing data that resembles
% the INPUTSET, where the sense of resemblance depends on the chosen 
% algorithm. 
%
% Optional arguments
% ------------------
%
% ALGORITHM ('algorithm') is a string which specifies the algorithm to
% use for generating the surrogate data set. Must be one of
%
% ### preserve_correlations
%
% Preserves the auto-correlations in a row and the 
% cross-correlations between the rows. If INPUTSET is
% a wide-sense stationary stochastic process, then so is
% SURROGATESET. If INPUTSET is a gaussian process, then so
% is SURROGATESET. This algorithm is described in
% "Generating Surrogate Data for Time series with Several 
% Simultaneously Measured Variables", Dean Prichard, James Theiler, 
% Physical Review Letters, Volume 73, Number 7, 1994. 
%
% ### preserve_distribution_and_correlations
%
% This algorithm is described for the 1-dimensional case in
% "Improved Surrogate Data for Nonlinearity Tests",
% Thomas Schreiber, Andreas Schmitz, Physical Review Letters,
% Volume 77, Number 4, 1996.
%
% ### preserve_dynamics
%
% Preserves the underlying dynamical system. The surrogate
% data is generated according to the examples set by the
% input set. It is assumed that the input data already resides
% in a state space which represents the dynamical system well
% (i.e. has been delay-embedded properly).
%
% Optional arguments for correlation-preserving algorithms
% --------------------------------------------------------
%
% FREQUENCY_RANGE ('frequencyRange') is a pair [a, b] of real numbers, 
% which denotes the interval [a, b] of normalized frequencies whose 
% phases to randomize. The normalized frequency 1 corresponds to half 
% the sampling rate. It must hold that 0 <= a <= b <= 1. This allows
% for truncated randomization, as described in xxx.
% Default: [0, 1]
%
% Optional arguments for dynamics-preserving algorithms
% -----------------------------------------------------
%
% K ('k') is a positive integer which specifies the
% number of nearest neighbors to use for prediction.

function surrogateSet = surrogate(inputSet, varargin)

import([tim_package, '.*']);

concept_check(nargin, 'inputs', 1);
concept_check(nargout, 'outputs', 0 : 1);

% Optional input arguments
algorithm = 'preserve_correlations';
frequencyRange = [0, 1];
k = 1;
eval(process_options({...
	'algorithm', ...
	'frequencyRange', ...
	'k'}, varargin));

[d, n] = size(inputSet);
minFrequency = frequencyRange(1);
maxFrequency = frequencyRange(2);
kNearest = k;

if strcmp(algorithm, 'preserve_dynamics')
	% Create a kd-tree from the points in the input set.
	kdTree = pastelgeometry.PointKdTree(d);
	idSet = kdTree.insert(inputSet);
	kdTree.refine();

	% Hide the last point in the time-series.
	% This is because it can not be used for prediction.
	kdTree.hide(idSet(end));

	% Pick a neighborhood of k points from the trajectory.
	startSet = kdTree.search_nearest(...
		kdTree.as_points(randi([1, n])), Inf, kNearest);

	% Choose the starting point as a random convex combination
	% of the points in the neighborhood.
	factorSet = rand(kNearest, 1);
	startPoint = kdTree.as_points(startSet) * (factorSet / sum(factorSet));

	% Generate surrogate data using non-linear prediction.
	currentPoint = startPoint;
	for i = 1 : n
		surrogateSet(:, i) = currentPoint;

		% Pick a neighborhood of k points from the trajectory
		% around the current point.
		neighborSet = kdTree.search_nearest(...
			currentPoint, Inf, kNearest);

		% Evolve the neighborhood in time. This is always
		% defined since we hided the last point in time.
		neighborSet = neighborSet + 1;

		% Choose the predicted point as a convex combination
		% of the predicted neighborhood.
		%factorSet = rand(kNearest, 1);
		factorSet = ones(kNearest, 1);
		currentPoint = kdTree.as_points(neighborSet) * ...
			(factorSet / sum(factorSet));
	end

	clear kdTree;
end

if strcmp(algorithm, 'preserve_correlations')
	% Let
	% 
	% 	X : [1, n] subset ZZ --> RR^d 
	%
	% be a stochastic process, such that
	%
	%	X(t) = [X_1(t), ... X_d(t)].
	%
	% Compute the discrete Fourier transform of X_i as 
	%
	%	F_i : [1, n] subset ZZ --> CC.
	%
	% Let R : [1, n] subset ZZ --> CC be a function such that
	% R(f) is a random unit complex number for all f in [1, n],
	% and R(mod(-f + 1, n) + 1) = -R(f). The surrogate data 
	% Y_i : [1, n] subset ZZ --> RR is then given by
	%
	%	Y_i = T^{-1}[abs(F_i) * R],
	%
	% where T^{-1} is the inverse discrete Fourier transform.

	nHalf = floor((n - 1) / 2);
	minIndex = 1 + round((nHalf - 1) * minFrequency);
	maxIndex = 1 + round((nHalf - 1) * maxFrequency);

	% Only the phases for the requested frequency-range
	% will be randomized.
	halfAngleSet = [...
		zeros(1, minIndex - 1), ...
		rand(1, maxIndex - minIndex + 1) * 2 * pi, ...
		zeros(1, nHalf - maxIndex)];

	% It can be proved that the result of the inverse discrete 
	% Fourier transform is real if and only if
	%
	% 	|F(N - k)| = |F(k)|, and
	% 	arg(F(N - k)) = -arg(F(k)), 
	%
	% for all k in [1, N - 1]. Since the discrete Fourier transform
	% is invertible, also the discrete Fourier transform of a 
	% real sequence has these properties.
	%
	% Since the surrogate data was requested to be real,
	% the phase-randomization must be done under this 
	% constraint. The value of arg(F(0)) does not matter;
	% it is the phase of the mean.

	% Enforce the real-condition on the angle-shifts.
	if mod(n, 2) == 1
		angleSet = [0, halfAngleSet, -fliplr(halfAngleSet)];
	else
		angleSet = [0, halfAngleSet, 0, -fliplr(halfAngleSet)];
	end

	surrogateSet = zeros(d, n);
	for i = 1 : d
		% The k:th component of the discrete Fourier transform of
		% inputSet(i, :) is given alpha_k addition in angle. This 
		% addition is independent of i, which guarantees that the 
		% cross-correlations are preserved.
		randomizedFftSet = fft(inputSet(i, :)) .* exp(j * angleSet);
		surrogateSet(i, :) = real(ifft(randomizedFftSet));
	end
end

