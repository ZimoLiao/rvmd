function run_octave_tests
%RUN_OCTAVE_TESTS Executable regression suite for GNU Octave.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);

testFilterCoefficient;
testEndpointMultiplicity;
testWeight;
testZeroInput;
testComplexPhase;
testRestart;
testLegacyRestart;
testRealFrequencies;
testComplexFrequencies;
testSinglePrecisionDC;
testInvalidInputs;

fprintf('All RVMD Octave regression tests passed.\n');
end

function testEndpointMultiplicity
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
assert(norm(state.phi_n - expected) < 1e-12);
end

function testFilterCoefficient
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
assert(norm(state.c_spec_n - expected) / norm(expected) < 1e-12);
end

function testWeight
rand('seed', 11);
randn('seed', 11);
q = randn(3, 24);
weight = [1, 2, 5];
[mode, info] = rvmd(q, 2, 20, ...
    'Weight', weight, 'MaximumSteps', 2, 'Tolerance', 0, ...
    'FPPrecision', 'double');
expectedWeight = weight(:) / mean(weight);
weightedNorm = sqrt(sum(abs(mode.phi).^2 .* expectedWeight, 1));
assert(max(abs(weightedNorm - 1)) < 1e-12);
assert(norm(info.weight - expectedWeight) == 0);
end

function testZeroInput
[mode, info, state] = rvmd(zeros(3, 16), 2, 10, ...
    'MaximumSteps', 3, 'Tolerance', 0, 'FPPrecision', 'double');
assert(all(isfinite(mode.phi(:))));
assert(all(isfinite(mode.c(:))));
assert(all(isfinite(mode.omega(:))));
assert(all(isfinite(info.Iteration.difference(:))));
assert(all(isfinite(state.c_spec_n(:))));
assert(norm(mode.c, 'fro') < 1e-14);
end

function testComplexPhase
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
assert(real(mode.phi(pivot, 1)) > 0);
assert(abs(imag(mode.phi(pivot, 1))) < 1e-12);
assert(norm(q - mode.phi * mode.c.', 'fro') / norm(q, 'fro') < 1e-12);
end

function testRestart
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
assert(norm(modeRestart.phi - modeFull.phi, 'fro') < 1e-12);
assert(norm(modeRestart.c - modeFull.c, 'fro') < 1e-12);
assert(norm(modeRestart.omega - modeFull.omega) < 1e-12);
assert(norm(infoRestart.Iteration.omega - infoFull.Iteration.omega, 'fro') < 1e-12);
assert(norm(infoRestart.Iteration.difference - ...
    infoFull.Iteration.difference) < 1e-12);
end

function testLegacyRestart
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
assert(norm(modeLegacy.phi - modeFull.phi, 'fro') < 1e-12);
assert(norm(modeLegacy.c - modeFull.c, 'fro') < 1e-12);
assert(norm(modeLegacy.omega - modeFull.omega) < 1e-12);
end

function testRealFrequencies
T = 128;
n = 0:(T - 1);
q = [1; -0.4; 0.7] * cos(2 * pi * 8 * n / T) ...
    + [0.2; 1; -0.5] * sin(2 * pi * 24 * n / T);
[mode, ~] = rvmd(q, 2, 500, ...
    'MaximumSteps', 100, 'Tolerance', 1e-10, ...
    'InitFreqType', 1, 'InitFreqMaximum', 0.2, ...
    'FPPrecision', 'double');
assert(max(abs(mode.omega - [8; 24] / T)) < 2 / T);
assert(isreal(mode.phi));
assert(isreal(mode.c));
end

function testComplexFrequencies
T = 128;
n = 0:(T - 1);
q = [1; 0.4i] * exp(2i * pi * 9 * n / T) ...
    + [0.2i; 1] * exp(-2i * pi * 27 * n / T);
[mode, ~] = rvmd(q, 2, 500, ...
    'MaximumSteps', 100, 'Tolerance', 1e-10, ...
    'InitFreqType', 1, 'InitFreqMaximum', 0.25, ...
    'FPPrecision', 'double');
assert(max(abs(mode.omega - [9; 27] / T)) < 2 / T);
assert(~isreal(mode.phi));
assert(~isreal(mode.c));
end

function testSinglePrecisionDC
T = 32;
n = 0:(T - 1);
q = single([ones(1, T); cos(2 * pi * 4 * n / T)]);
[mode, info] = rvmd(q, 2, 40, ...
    'nDC', 1, 'MaximumSteps', 5, 'Tolerance', 0, ...
    'FPPrecision', 'single');
assert(isa(mode.phi, 'single'));
assert(isa(mode.c, 'single'));
assert(isa(mode.omega, 'single'));
assert(all(info.Iteration.omega(1, :) == 0));
assert(mode.omega(1) == 0);
end

function testInvalidInputs
expectError(@() rvmd(randn(2, 8), 2, 10, 'Weight', [1, -1]), ...
    'rvmd:InvalidWeight');
expectError(@() rvmd(randn(2, 8), 2, 10, 'nDC', 3), ...
    'rvmd:InvalidNDC');
expectError(@() rvmd(randn(2, 8), 2, 10, 'FPPrecision', 'half'), ...
    'rvmd:InvalidPrecision');
end

function expectError(callable, identifier)
try
    callable();
catch exception
    assert(strcmp(exception.identifier, identifier));
    return
end
error('rvmd:TestFailure', 'Expected error %s was not thrown.', identifier);
end

function qExtended = mirrorExtend(q)
T = size(q, 2);
half = ceil(T / 2);
qExtended = zeros(size(q, 1), 2 * T, class(q));
qExtended(:, 1:half) = q(:, half:-1:1);
qExtended(:, (half + 1):(half + T)) = q;
qExtended(:, (half + T + 1):end) = q(:, T:-1:(half + 1));
end
