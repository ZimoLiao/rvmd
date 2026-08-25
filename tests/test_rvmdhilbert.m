function tests = test_rvmdhilbert
% Regression tests for RVMD Hilbert spectral analysis.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.OriginalPath = path;
repositoryRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repositoryRoot);
end

function teardownOnce(testCase)
path(testCase.TestData.OriginalPath);
end

function testRealToneInstantaneousFrequency(testCase)
sampleRate = 256;
sampleCount = 256;
frequency = 32;
time = (0:(sampleCount - 1)).' / sampleRate;
mode.c = cos(2 * pi * frequency * time);
mode.omega = frequency / sampleRate;

analysis = rvmdhilbert(mode, sampleRate, ...
    'MirrorExtension', false, 'FrequencyBins', 128);

interior = 3:(sampleCount - 2);
verifyEqual(testCase, analysis.Amplitude(interior), ...
    ones(numel(interior), 1), 'AbsTol', 1e-10);
verifyEqual(testCase, analysis.InstantaneousFrequency(interior), ...
    frequency * ones(numel(interior), 1), 'AbsTol', 1e-10);
verifyEqual(testCase, size(analysis.HilbertSpectrum), [128, sampleCount]);
verifyTrue(testCase, issparse(analysis.HilbertSpectrum));
[~, peak] = max(analysis.MarginalSpectrum);
verifyLessThanOrEqual(testCase, ...
    abs(analysis.FrequencyBins(peak) - frequency), sampleRate / 128);
end

function testComplexCoefficientKeepsSignedFrequency(testCase)
sampleRate = 200;
sampleCount = 200;
frequency = -25;
time = (0:(sampleCount - 1)).' / sampleRate;
mode.c = exp(2i * pi * frequency * time);
mode.omega = abs(frequency) / sampleRate;

analysis = rvmdhilbert(mode, sampleRate, 'FrequencyBins', 100);

verifyEqual(testCase, analysis.Method, 'complex-coefficient');
verifyEqual(testCase, analysis.FrequencyLimits, [-100, 100], 'AbsTol', 0);
verifyEqual(testCase, analysis.InstantaneousFrequency(3:end-2), ...
    frequency * ones(sampleCount - 4, 1), 'AbsTol', 1e-10);
end

function testZeroAmplitudeIsMasked(testCase)
mode.c = zeros(32, 2);
mode.omega = [0; 0.1];

analysis = rvmdhilbert(mode, 10);

verifyTrue(testCase, all(isnan(analysis.InstantaneousFrequency(:))));
verifyEqual(testCase, nnz(analysis.HilbertSpectrum), 0);
verifyEqual(testCase, analysis.ModeEnergy, [0, 0], 'AbsTol', 0);
end

function testModeSelection(testCase)
sampleRate = 64;
time = (0:63).' / sampleRate;
mode.c = [cos(2 * pi * 4 * time), cos(2 * pi * 12 * time)];
mode.omega = [4; 12] / sampleRate;

analysis = rvmdhilbert(mode, sampleRate, 'ModeIndices', 2, ...
    'MirrorExtension', false);

verifyEqual(testCase, analysis.ModeIndices, 2);
verifyEqual(testCase, size(analysis.AnalyticSignal), [64, 1]);
verifyEqual(testCase, analysis.CenterFrequencies, 12, 'AbsTol', 1e-12);
end

function testInvalidInputsAreRejected(testCase)
mode.c = randn(16, 2);
mode.omega = [0.1; 0.2];
verifyError(testCase, @() rvmdhilbert(mode, 0), ...
    'rvmdhilbert:InvalidSampleRate');
verifyError(testCase, @() rvmdhilbert(mode, 1, 'ModeIndices', [1, 1]), ...
    'rvmdhilbert:InvalidModeIndices');
verifyError(testCase, @() rvmdhilbert(mode, 1, 'FrequencyLimits', [1, 0]), ...
    'rvmdhilbert:InvalidFrequencyLimits');
end
