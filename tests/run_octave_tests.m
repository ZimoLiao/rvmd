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
testWeightedSingleRestart;
testRestartDefaults;
testLegacyRestart;
testLegacyComplexSingleSnapshot;
testRealFrequencies;
testComplexFrequencies;
testSinglePrecisionDC;
testOutputFunctionStop;
testTimeLimit;
testCheckpoint;
testInitialFrequencies;
testHilbertAnalysis;
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

function testWeightedSingleRestart
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
assert(isequal(modeRestart.phi, modeFull.phi));
assert(isequal(modeRestart.c, modeFull.c));
assert(isequal(modeRestart.omega, modeFull.omega));
assert(isequal(infoRestart.Iteration.omega, infoFull.Iteration.omega));
end

function testRestartDefaults
T = 32;
n = 0:(T - 1);
q = [cos(2 * pi * 3 * n / T); sin(2 * pi * 7 * n / T)];
common = {'Tolerance', 0, 'InitFreqType', 0, ...
    'FPPrecision', 'double'};
[modeFull, infoFull] = rvmd(q, 2, 30, common{:}, 'MaximumSteps', 5);
[~, ~, state] = rvmd(q, 2, 30, common{:}, 'MaximumSteps', 2);
[modeRestart, infoRestart] = rvmd('Restart', state, 'MaximumSteps', 5);
assert(isequal(modeRestart.phi, modeFull.phi));
assert(isequal(modeRestart.c, modeFull.c));
assert(infoRestart.Tolerance == infoFull.Tolerance);
assert(strcmp(infoRestart.FPPrecision, 'double'));
assert(infoRestart.InitFreqType == 0);
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
legacy.Q = cast(q, state.FPPrecision);
legacy.steps = state.Iteration.steps;
legacy.omega = state.Iteration.omega;
legacy.difference = state.Iteration.difference;
legacy = rmfield(legacy, {'version', 'dataSpectrumNorm', ...
    'InitialFrequencies', 'DisplayInterval', 'Tolerance', 'MaximumSteps', ...
    'InitFreqType', 'InitFreqMaximum', 'Device', 'FPPrecision', ...
    'nDC', 'isRealInput', 'Display', 'Iteration'});
[modeLegacy, ~] = rvmd('Restart', legacy, ...
    'MaximumSteps', 5, 'Tolerance', 0);
assert(norm(modeLegacy.phi - modeFull.phi, 'fro') < 1e-12);
assert(norm(modeLegacy.c - modeFull.c, 'fro') < 1e-12);
assert(norm(modeLegacy.omega - modeFull.omega) < 1e-12);
end

function testLegacyComplexSingleSnapshot
q = [1 + 2i; -0.3 + 0.7i];
[modeFull, ~] = rvmd(q, 1, 0, ...
    'MaximumSteps', 3, 'Tolerance', 0, 'FPPrecision', 'double');
[~, ~, state] = rvmd(q, 1, 0, ...
    'MaximumSteps', 1, 'Tolerance', 0, 'FPPrecision', 'double');
legacy = state;
legacy.Q = cast(q, state.FPPrecision);
legacy = rmfield(legacy, {'version', 'residual_n', 'dataSpectrumNorm', ...
    'InitialFrequencies', 'DisplayInterval', 'isRealInput'});
[modeLegacy, ~] = rvmd('Restart', legacy, ...
    'MaximumSteps', 3, 'Tolerance', 0);
assert(~isreal(modeLegacy.phi));
assert(norm(modeLegacy.phi - modeFull.phi, 'fro') < 1e-12);
assert(norm(modeLegacy.c - modeFull.c, 'fro') < 1e-12);
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

function testOutputFunctionStop
T = 48;
n = 0:(T - 1);
q = [cos(2 * pi * 3 * n / T) + 0.2 * cos(2 * pi * 9 * n / T); ...
     sin(2 * pi * 3 * n / T) - 0.1 * cos(2 * pi * 9 * n / T)];
common = {'Tolerance', 0, 'InitFreqType', 1, ...
    'InitFreqMaximum', 0.25, 'FPPrecision', 'double'};
stopAfterThree = @(progress, phase) ...
    strcmp(phase, 'iter') && progress.step >= 3;
[modeFull, ~] = rvmd(q, 2, 80, common{:}, 'MaximumSteps', 6);
[~, infoStopped, state] = rvmd(q, 2, 80, common{:}, ...
    'MaximumSteps', 10, 'OutputFcn', stopAfterThree);
[modeRestart, ~] = rvmd('Restart', state, ...
    'Tolerance', 0, 'MaximumSteps', 6);
assert(infoStopped.Iteration.steps == 3);
assert(infoStopped.ExitFlag == -1);
assert(strcmp(infoStopped.StopReason, 'outputFunction'));
assert(norm(modeRestart.phi - modeFull.phi, 'fro') < 1e-12);
assert(norm(modeRestart.c - modeFull.c, 'fro') < 1e-12);
end

function testTimeLimit
q = reshape(sin(2 * pi * (0:31) / 8), 2, 16);
common = {'Tolerance', 0, 'InitFreqType', 1, ...
    'FPPrecision', 'double'};
[modeFull, ~] = rvmd(q, 2, 20, common{:}, 'MaximumSteps', 4);
[~, infoStopped, state] = rvmd(q, 2, 20, common{:}, ...
    'MaximumSteps', 4, 'TimeLimit', 0);
[modeRestart, ~] = rvmd('Restart', state, 'MaximumSteps', 4);
assert(infoStopped.Iteration.steps == 0);
assert(infoStopped.ExitFlag == -2);
assert(norm(modeRestart.phi - modeFull.phi, 'fro') < 1e-12);
assert(norm(modeRestart.c - modeFull.c, 'fro') < 1e-12);
end

function testCheckpoint
checkpointFile = [tempname, '.mat'];
unwind_protect
    T = 40;
    n = 0:(T - 1);
    q = [cos(2 * pi * 3 * n / T); sin(2 * pi * 7 * n / T)];
    common = {'Tolerance', 0, 'FPPrecision', 'double'};
    [modeFull, ~] = rvmd(q, 2, 30, common{:}, 'MaximumSteps', 6);
    [~, infoPart] = rvmd(q, 2, 30, common{:}, ...
        'MaximumSteps', 3, 'CheckpointFile', checkpointFile, ...
        'CheckpointInterval', 2);
    saved = load(checkpointFile, 'restart');
    previous = load([checkpointFile, '.prev'], 'restart');
    [modeRestart, ~] = rvmd('Restart', saved.restart, ...
        'Tolerance', 0, 'MaximumSteps', 6);
    assert(infoPart.LastCheckpointStep == 3);
    assert(saved.restart.version == 4);
    assert(~isfield(saved.restart, 'Q'));
    assert(previous.restart.Iteration.steps == 2);
    assert(norm(modeRestart.phi - modeFull.phi, 'fro') < 1e-12);
    assert(norm(modeRestart.c - modeFull.c, 'fro') < 1e-12);
unwind_protect_cleanup
    deleteCheckpointFiles(checkpointFile);
end_unwind_protect
end

function testInitialFrequencies
q = randn(3, 24);
initial = [0; 0.07; 0.21];
[~, info] = rvmd(q, 3, 20, ...
    'nDC', 1, 'InitialFrequencies', initial, ...
    'MaximumSteps', 1, 'Tolerance', 0, 'FPPrecision', 'double');
assert(isequal(info.InitialFrequencies, initial));
assert(isequal(info.Iteration.omega(:, 1), initial));
end

function testHilbertAnalysis
sampleRate = 256;
sampleCount = 256;
frequency = 32;
time = (0:(sampleCount - 1)).' / sampleRate;
mode.c = cos(2 * pi * frequency * time);
mode.omega = frequency / sampleRate;
analysis = rvmdhilbert(mode, sampleRate, ...
    'MirrorExtension', false, 'FrequencyBins', 128);
interior = 3:(sampleCount - 2);
assert(max(abs(analysis.Amplitude(interior) - 1)) < 1e-10);
assert(max(abs(analysis.InstantaneousFrequency(interior) - frequency)) < 1e-10);
assert(issparse(analysis.HilbertSpectrum));
assert(isequal(size(analysis.HilbertSpectrum), [128, sampleCount]));

complexFrequency = -25;
complexMode.c = exp(2i * pi * complexFrequency * time);
complexMode.omega = abs(complexFrequency) / sampleRate;
complexAnalysis = rvmdhilbert(complexMode, sampleRate, ...
    'FrequencyBins', 128);
assert(strcmp(complexAnalysis.Method, 'complex-coefficient'));
assert(max(abs(complexAnalysis.InstantaneousFrequency(interior) - ...
    complexFrequency)) < 1e-10);
end

function testInvalidInputs
expectError(@() rvmd(randn(2, 8), 2, 10, 'Weight', [1, -1]), ...
    'rvmd:InvalidWeight');
expectError(@() rvmd(randn(2, 8), 2, 10, 'nDC', 3), ...
    'rvmd:InvalidNDC');
expectError(@() rvmd(randn(2, 8), 2, 10, 'FPPrecision', 'half'), ...
    'rvmd:InvalidPrecision');
expectError(@() rvmd(randn(2, 8), 2, 10, 'Restart', struct()), ...
    'rvmd:RestartSyntax');
[~, ~, state] = rvmd(randn(2, 8), 2, 10, ...
    'MaximumSteps', 2, 'Tolerance', 0);
expectError(@() rvmd('Restart', state, 'MaximumSteps', 1), ...
    'rvmd:InvalidMaximumSteps');
expectError(@() rvmd(randn(2, 8), 2, 10, ...
    'Restart', state), 'rvmd:RestartSyntax');
state.version = 999;
expectError(@() rvmd('Restart', state), ...
    'rvmd:UnsupportedRestartVersion');
end

function deleteCheckpointFiles(checkpointFile)
files = {checkpointFile, [checkpointFile, '.prev'], ...
    [checkpointFile, '.tmp']};
for index = 1:numel(files)
    if exist(files{index}, 'file') == 2
        delete(files{index});
    end
end
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
