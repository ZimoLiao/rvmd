function tests = test_rvmd
% Regression tests for the real- and complex-valued RVMD implementations.
tests = functiontests(localfunctions);
end

function testFilterCoefficientMatchesObjective(testCase)
T = 16;
n = 0:(T - 1);
q = 2 + cos(2 * pi * 2 * n / T);
alpha = 37;

[~, ~, state] = rvmd(q, 1, alpha, ...
    'MaximumSteps', 1, 'Tolerance', 0, 'InitFreqType', 0, ...
    'FPPrecision', 'double');

qExtended = mirrorExtend(q);
qSpectrum = fft(qExtended, [], 2);
qSpectrum = qSpectrum(:, 1:(T + 1));
frequency = (0:T) / (2 * T);
expected = qSpectrum.' ./ (1 + 2 * alpha * frequency.^2).';

verifyLessThan(testCase, ...
    norm(state.c_spec_n - expected) / norm(expected), 1e-12);
end

function testOneSidedProjectionUsesEndpointMultiplicity(testCase)
T = 12;
n = 0:(T - 1);
q = [2 + cos(2 * pi * n / T); ...
     (-1).^n + 0.2 * sin(4 * pi * n / T); ...
     0.3 + cos(6 * pi * n / T)];
weight = [1, 3, 2];

[~, ~, state] = rvmd(q, 1, 0, ...
    'Weight', weight, 'MaximumSteps', 1, 'Tolerance', 0, ...
    'InitFreqType', 0, 'FPPrecision', 'double');

qSpectrum = fft(mirrorExtend(q), [], 2);
qSpectrum = qSpectrum(:, 1:(T + 1));
frequencyWeight = ones(T + 1, 1);
frequencyWeight([1, end]) = 0.5;
expected = real(qSpectrum * frequencyWeight);
normalizedWeight = weight(:) / mean(weight);
expected = expected / sqrt(sum(abs(expected).^2 .* normalizedWeight));
[~, pivot] = max(abs(expected));
if expected(pivot) < 0
    expected = -expected;
end

verifyEqual(testCase, state.phi_n, expected, 'AbsTol', 1e-12);
end

function testRowWeightUsesWeightedUnitNorm(testCase)
rng(11);
q = randn(3, 24);
weight = [1, 2, 5];

[mode, info] = rvmd(q, 2, 20, ...
    'Weight', weight, 'MaximumSteps', 2, 'Tolerance', 0, ...
    'FPPrecision', 'double');

expectedWeight = weight(:) / mean(weight);
weightedNorm = sqrt(sum(abs(mode.phi).^2 .* expectedWeight, 1));
verifyEqual(testCase, weightedNorm, ones(1, 2), 'AbsTol', 1e-12);
verifyEqual(testCase, info.weight, expectedWeight, 'AbsTol', 0);
end

function testZeroInputRemainsFinite(testCase)
[mode, info, state] = rvmd(zeros(3, 16), 2, 10, ...
    'MaximumSteps', 3, 'Tolerance', 0, 'FPPrecision', 'double');

verifyTrue(testCase, all(isfinite(mode.phi(:))));
verifyTrue(testCase, all(isfinite(mode.c(:))));
verifyTrue(testCase, all(isfinite(mode.omega(:))));
verifyTrue(testCase, all(isfinite(info.Iteration.difference(:))));
verifyTrue(testCase, all(isfinite(state.c_spec_n(:))));
verifyEqual(testCase, mode.c, zeros(size(mode.c)), 'AbsTol', 1e-14);
end

function testComplexModeHasCanonicalPhaseAndReconstructs(testCase)
T = 32;
n = 0:(T - 1);
spatial = [1 + 2i; -0.3 + 0.7i; 0.2 - 0.1i];
coefficient = exp(2i * pi * 3 * n / T) ...
    + 0.25 * exp(-2i * pi * 7 * n / T);
q = spatial * coefficient;

[mode, ~] = rvmd(q, 1, 0, ...
    'MaximumSteps', 2, 'Tolerance', 0, 'InitFreqType', 0, ...
    'FPPrecision', 'double');

[~, pivot] = max(abs(mode.phi(:, 1)));
verifyGreaterThan(testCase, real(mode.phi(pivot, 1)), 0);
verifyEqual(testCase, imag(mode.phi(pivot, 1)), 0, 'AbsTol', 1e-12);
verifyLessThan(testCase, ...
    norm(q - mode.phi * mode.c.', 'fro') / norm(q, 'fro'), 1e-12);
end

function testRestartMatchesSingleRun(testCase)
T = 48;
n = 0:(T - 1);
q = [cos(2 * pi * 3 * n / T) + 0.2 * cos(2 * pi * 9 * n / T); ...
     sin(2 * pi * 3 * n / T) - 0.1 * cos(2 * pi * 9 * n / T)];
common = {'Tolerance', 0, 'InitFreqType', 1, ...
    'InitFreqMaximum', 0.25, 'FPPrecision', 'double'};

[modeFull, infoFull] = rvmd(q, 2, 80, common{:}, 'MaximumSteps', 6);
[~, ~, state] = rvmd(q, 2, 80, common{:}, 'MaximumSteps', 3);
[modeRestart, infoRestart] = rvmd('Restart', state, ...
    'Tolerance', 0, 'MaximumSteps', 6);

verifyEqual(testCase, modeRestart.phi, modeFull.phi, 'AbsTol', 1e-12);
verifyEqual(testCase, modeRestart.c, modeFull.c, 'AbsTol', 1e-12);
verifyEqual(testCase, modeRestart.omega, modeFull.omega, 'AbsTol', 1e-12);
verifyEqual(testCase, infoRestart.Iteration.omega, ...
    infoFull.Iteration.omega, 'AbsTol', 1e-12);
verifyEqual(testCase, infoRestart.Iteration.difference, ...
    infoFull.Iteration.difference, 'AbsTol', 1e-12);
end

function testWeightedSingleRestartIsBitwiseEquivalent(testCase)
T = 40;
n = 0:(T - 1);
q = single([cos(2 * pi * 3 * n / T); ...
    sin(2 * pi * 7 * n / T); ...
    0.4 * cos(2 * pi * 11 * n / T)]);
common = {'Weight', [0.1, 0.3, 2.7], 'Tolerance', 0, ...
    'InitFreqMaximum', 0.3, 'FPPrecision', 'single'};

[modeFull, infoFull] = rvmd(q, 3, 40, common{:}, 'MaximumSteps', 9);
[~, ~, state] = rvmd(q, 3, 40, common{:}, 'MaximumSteps', 4);
[modeRestart, infoRestart] = rvmd('Restart', state, ...
    'Tolerance', 0, 'MaximumSteps', 9);

verifyEqual(testCase, modeRestart.phi, modeFull.phi, 'AbsTol', 0);
verifyEqual(testCase, modeRestart.c, modeFull.c, 'AbsTol', 0);
verifyEqual(testCase, modeRestart.omega, modeFull.omega, 'AbsTol', 0);
verifyEqual(testCase, infoRestart.Iteration.omega, ...
    infoFull.Iteration.omega, 'AbsTol', 0);
end

function testLegacyDevelopRestartIsAccepted(testCase)
T = 32;
n = 0:(T - 1);
q = [cos(2 * pi * 3 * n / T); sin(2 * pi * 7 * n / T)];
[modeFull, ~] = rvmd(q, 2, 30, ...
    'MaximumSteps', 5, 'Tolerance', 0, 'FPPrecision', 'double');
[~, ~, state] = rvmd(q, 2, 30, ...
    'MaximumSteps', 2, 'Tolerance', 0, 'FPPrecision', 'double');

legacy = state;
legacy.steps = state.Iteration.steps;
legacy.omega = state.Iteration.omega;
legacy.difference = state.Iteration.difference;
legacy = rmfield(legacy, {'version', 'Tolerance', 'MaximumSteps', ...
    'InitFreqType', 'InitFreqMaximum', 'Device', 'FPPrecision', ...
    'nDC', 'isRealInput', 'Display', 'Iteration'});

[modeLegacy, ~] = rvmd('Restart', legacy, ...
    'MaximumSteps', 5, 'Tolerance', 0);
verifyEqual(testCase, modeLegacy.phi, modeFull.phi, 'AbsTol', 1e-12);
verifyEqual(testCase, modeLegacy.c, modeFull.c, 'AbsTol', 1e-12);
verifyEqual(testCase, modeLegacy.omega, modeFull.omega, 'AbsTol', 1e-12);
end

function testRealTwoToneFrequencies(testCase)
T = 128;
n = 0:(T - 1);
q = [1; -0.4; 0.7] * cos(2 * pi * 8 * n / T) ...
    + [0.2; 1; -0.5] * sin(2 * pi * 24 * n / T);

[mode, ~] = rvmd(q, 2, 500, ...
    'MaximumSteps', 100, 'Tolerance', 1e-10, ...
    'InitFreqType', 1, 'InitFreqMaximum', 0.2, ...
    'FPPrecision', 'double');

verifyEqual(testCase, mode.omega, [8; 24] / T, 'AbsTol', 2 / T);
verifyTrue(testCase, isreal(mode.phi));
verifyTrue(testCase, isreal(mode.c));
end

function testComplexSignedFrequenciesUseMagnitude(testCase)
T = 128;
n = 0:(T - 1);
q = [1; 0.4i] * exp(2i * pi * 9 * n / T) ...
    + [0.2i; 1] * exp(-2i * pi * 27 * n / T);

[mode, ~] = rvmd(q, 2, 500, ...
    'MaximumSteps', 100, 'Tolerance', 1e-10, ...
    'InitFreqType', 1, 'InitFreqMaximum', 0.25, ...
    'FPPrecision', 'double');

verifyEqual(testCase, mode.omega, [9; 27] / T, 'AbsTol', 2 / T);
verifyFalse(testCase, isreal(mode.phi));
verifyFalse(testCase, isreal(mode.c));
end

function testSinglePrecisionAndFixedDCMode(testCase)
T = 32;
n = 0:(T - 1);
q = single([ones(1, T); cos(2 * pi * 4 * n / T)]);
[mode, info] = rvmd(q, 2, 40, ...
    'nDC', 1, 'MaximumSteps', 5, 'Tolerance', 0, ...
    'FPPrecision', 'single');

verifyClass(testCase, mode.phi, 'single');
verifyClass(testCase, mode.c, 'single');
verifyClass(testCase, mode.omega, 'single');
verifyEqual(testCase, info.Iteration.omega(1, :), ...
    zeros(1, info.Iteration.steps + 1, 'single'), 'AbsTol', 0);
verifyEqual(testCase, mode.omega(1), single(0), 'AbsTol', 0);
end

function testInvalidInputsAreRejected(testCase)
verifyError(testCase, @() rvmd(randn(2, 8), 2, 10, ...
    'Weight', [1, -1]), 'rvmd:InvalidWeight');
verifyError(testCase, @() rvmd(randn(2, 8), 2, 10, ...
    'nDC', 3), 'rvmd:InvalidNDC');
verifyError(testCase, @() rvmd(randn(2, 8), 2, 10, ...
    'FPPrecision', 'half'), 'rvmd:InvalidPrecision');
verifyError(testCase, @() rvmd(randn(2, 8), 2, 10, ...
    'Restart', struct()), 'rvmd:InvalidRestart');

[~, ~, state] = rvmd(randn(2, 8), 2, 10, ...
    'MaximumSteps', 2, 'Tolerance', 0);
verifyError(testCase, @() rvmd('Restart', state, ...
    'MaximumSteps', 1), 'rvmd:InvalidMaximumSteps');
end

function qExtended = mirrorExtend(q)
T = size(q, 2);
half = ceil(T / 2);
qExtended = zeros(size(q, 1), 2 * T, 'like', q);
qExtended(:, 1:half) = q(:, half:-1:1);
qExtended(:, (half + 1):(half + T)) = q;
qExtended(:, (half + T + 1):end) = q(:, T:-1:(half + 1));
end
