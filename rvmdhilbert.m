function analysis = rvmdhilbert(mode, sampleRate, varargin)
%RVMDHILBERT Hilbert spectral analysis of RVMD time coefficients.
%   ANALYSIS = RVMDHILBERT(MODE, FS) analyzes MODE.C returned by RVMD. MODE.C
%   is T-by-K, MODE.OMEGA has K elements, and FS is a finite positive sample
%   rate in samples per time unit. Frequencies in ANALYSIS use the resulting
%   inverse-time unit.
%
%   ANALYSIS = RVMDHILBERT(MODE, FS, NAME, VALUE) accepts:
%       'ModeIndices'        Unique mode indices. Default: all modes.
%       'MirrorExtension'    Logical or finite numeric scalar; zero is false
%                            and nonzero is true. For real coefficients with
%                            more than two samples, use a half-length
%                            reflection on each side. Default: true.
%       'AmplitudeThreshold' Finite nonnegative fraction of each mode's
%                            maximum amplitude. Lower-amplitude instantaneous
%                            frequencies are set to NaN and excluded from the
%                            Hilbert spectrum. Default: 1e-6.
%       'FrequencyBins'      Positive integer bin count. Default: 256.
%       'FrequencyLimits'    Increasing finite two-vector. The default is
%                            [0,FS/2] for real coefficients and
%                            [-FS/2,FS/2] for complex coefficients.
%
%   ANALYSIS contains the following fields, where M is the number of selected
%   modes and B is the FrequencyBins option:
%       SampleRate             Scalar FS.
%       Time                   T-by-1 sample times beginning at zero.
%       ModeIndices            Selected mode indices.
%       CenterFrequencies      Selected MODE.OMEGA values multiplied by FS.
%       Method                 'hilbert', 'hilbert-mirror', or
%                              'complex-coefficient'.
%       AnalyticSignal         T-by-M selected analytic coefficients.
%       Amplitude              ABS(AnalyticSignal).
%       Phase                  Time-unwrapped phase in radians.
%       InstantaneousFrequency T-by-M phase derivative divided by 2*pi.
%       InstantaneousEnergy    Amplitude.^2.
%       FrequencyEdges         (B+1)-by-1 bin edges.
%       FrequencyBins          B-by-1 bin centers.
%       FrequencyLimits        Two-element frequency interval.
%       HilbertSpectrum        Sparse B-by-T matrix. At each time,
%                              selected-mode energies in the same bin are summed.
%       MarginalSpectrum       SUM(HilbertSpectrum,2)/FS.
%       ModeEnergy             SUM(ABS(selected MODE.C).^2,1)/FS.
%       HilbertEnergy          SUM(InstantaneousEnergy,1)/FS.
%       AmplitudeThreshold     Applied relative threshold.
%       MirrorExtension        Whether reflection was actually applied.
%
%   No Signal Processing Toolbox is required. For real coefficients the
%   analytic signal is constructed directly in the Fourier domain. Complex
%   coefficients are treated as already analytic-valued signals, so their
%   instantaneous frequencies may be signed.
%
%   RVMDHILBERT(MODE, FS, ...) with no output calls RVMDHILBERTPLOT.

p = inputParser;
p.FunctionName = 'rvmdhilbert';
addRequired(p, 'mode');
addRequired(p, 'sampleRate');
addParameter(p, 'ModeIndices', []);
addParameter(p, 'MirrorExtension', true);
addParameter(p, 'AmplitudeThreshold', 1e-6);
addParameter(p, 'FrequencyBins', 256);
addParameter(p, 'FrequencyLimits', []);
parse(p, mode, sampleRate, varargin{:});

validateMode(mode);
sampleRate = validatePositiveScalar(p.Results.sampleRate, 'SampleRate');

coefficient = mode.c;
[sampleCount, modeCount] = size(coefficient);
modeIndices = validateModeIndices(p.Results.ModeIndices, modeCount);
coefficient = coefficient(:, modeIndices);

mirrorExtension = p.Results.MirrorExtension;
if ~(islogical(mirrorExtension) || isnumeric(mirrorExtension)) || ...
        ~isscalar(mirrorExtension) || ~isfinite(double(mirrorExtension))
    error('rvmdhilbert:InvalidMirrorExtension', ...
        'MirrorExtension must be a logical scalar.');
end
mirrorExtension = logical(mirrorExtension);

amplitudeThreshold = p.Results.AmplitudeThreshold;
if ~isnumeric(amplitudeThreshold) || ~isreal(amplitudeThreshold) || ...
        ~isscalar(amplitudeThreshold) || ~isfinite(amplitudeThreshold) || ...
        amplitudeThreshold < 0
    error('rvmdhilbert:InvalidAmplitudeThreshold', ...
        'AmplitudeThreshold must be a finite nonnegative scalar.');
end
amplitudeThreshold = double(amplitudeThreshold);

frequencyBinCount = p.Results.FrequencyBins;
if ~isnumeric(frequencyBinCount) || ~isreal(frequencyBinCount) || ...
        ~isscalar(frequencyBinCount) || ~isfinite(frequencyBinCount) || ...
        frequencyBinCount < 1 || frequencyBinCount ~= fix(frequencyBinCount)
    error('rvmdhilbert:InvalidFrequencyBins', ...
        'FrequencyBins must be a positive integer.');
end
frequencyBinCount = double(frequencyBinCount);

appliedMirrorExtension = false;
if isreal(coefficient)
    if mirrorExtension && sampleCount > 2
        padLength = min(floor(sampleCount / 2), sampleCount - 1);
        extended = [flipud(coefficient(2:(padLength + 1), :)); ...
                    coefficient; ...
                    flipud(coefficient((end - padLength):(end - 1), :))];
        analyticExtended = analyticSignal(extended);
        analytic = analyticExtended(padLength + (1:sampleCount), :);
        method = 'hilbert-mirror';
        appliedMirrorExtension = true;
    else
        analytic = analyticSignal(coefficient);
        method = 'hilbert';
    end
else
    analytic = coefficient;
    method = 'complex-coefficient';
end

amplitude = abs(analytic);
phase = unwrap(angle(analytic), [], 1);
instantaneousFrequency = phaseDerivative(phase, sampleRate);

perModeMaximum = max(amplitude, [], 1);
threshold = amplitudeThreshold .* perModeMaximum;
validAmplitude = amplitude >= threshold;
validAmplitude(:, perModeMaximum == 0) = false;
instantaneousFrequency(~validAmplitude) = nan;
instantaneousEnergy = amplitude .^ 2;

frequencyLimits = p.Results.FrequencyLimits;
if isempty(frequencyLimits)
    if isreal(coefficient)
        frequencyLimits = [0, sampleRate / 2];
    else
        frequencyLimits = [-sampleRate / 2, sampleRate / 2];
    end
else
    if ~isnumeric(frequencyLimits) || ~isreal(frequencyLimits) || ...
            numel(frequencyLimits) ~= 2 || ...
            ~all(isfinite(frequencyLimits(:))) || ...
            frequencyLimits(1) >= frequencyLimits(2)
        error('rvmdhilbert:InvalidFrequencyLimits', ...
            'FrequencyLimits must be an increasing finite two-vector.');
    end
    frequencyLimits = double(frequencyLimits(:).');
end

frequencyEdges = linspace( ...
    frequencyLimits(1), frequencyLimits(2), frequencyBinCount + 1);
frequencyCenters = (frequencyEdges(1:end-1) + frequencyEdges(2:end)) / 2;
hilbertSpectrum = buildHilbertSpectrum(instantaneousFrequency, ...
    instantaneousEnergy, frequencyLimits, frequencyBinCount);

analysis.SampleRate = sampleRate;
analysis.Time = (0:(sampleCount - 1)).' / sampleRate;
analysis.ModeIndices = modeIndices;
analysis.CenterFrequencies = double(mode.omega(modeIndices)) * sampleRate;
analysis.Method = method;
analysis.AnalyticSignal = analytic;
analysis.Amplitude = amplitude;
analysis.Phase = phase;
analysis.InstantaneousFrequency = instantaneousFrequency;
analysis.InstantaneousEnergy = instantaneousEnergy;
analysis.FrequencyEdges = frequencyEdges(:);
analysis.FrequencyBins = frequencyCenters(:);
analysis.FrequencyLimits = frequencyLimits;
analysis.HilbertSpectrum = hilbertSpectrum;
analysis.MarginalSpectrum = full(sum(hilbertSpectrum, 2)) / sampleRate;
analysis.ModeEnergy = sum(abs(double(coefficient)) .^ 2, 1) / sampleRate;
analysis.HilbertEnergy = sum(double(instantaneousEnergy), 1) / sampleRate;
analysis.AmplitudeThreshold = amplitudeThreshold;
analysis.MirrorExtension = appliedMirrorExtension;

if nargout == 0
    rvmdhilbertplot(analysis);
    clear analysis
end
end

function analytic = analyticSignal(signal)
sampleCount = size(signal, 1);
multiplier = zeros(sampleCount, 1, 'like', signal);
multiplier(1) = 1;
if mod(sampleCount, 2) == 0
    multiplier(2:(sampleCount / 2)) = 2;
    multiplier(sampleCount / 2 + 1) = 1;
else
    multiplier(2:((sampleCount + 1) / 2)) = 2;
end
analytic = ifft(fft(signal, [], 1) .* multiplier, [], 1);
end

function frequency = phaseDerivative(phase, sampleRate)
sampleCount = size(phase, 1);
frequency = zeros(size(phase), 'like', phase);
if sampleCount == 1
    frequency(:) = nan;
    return
end
scale = sampleRate / (2 * pi);
frequency(1, :) = (phase(2, :) - phase(1, :)) * scale;
frequency(end, :) = (phase(end, :) - phase(end - 1, :)) * scale;
if sampleCount > 2
    frequency(2:(end - 1), :) = ...
        (phase(3:end, :) - phase(1:(end - 2), :)) * (scale / 2);
end
end

function spectrum = buildHilbertSpectrum(frequency, energy, limits, binCount)
[sampleCount, modeCount] = size(frequency);
span = limits(2) - limits(1);
bin = floor((double(frequency) - limits(1)) / span * binCount) + 1;
bin(double(frequency) == limits(2)) = binCount;
valid = isfinite(frequency) & bin >= 1 & bin <= binCount & isfinite(energy);

column = repmat((1:sampleCount).', 1, modeCount);
spectrum = sparse(bin(valid), column(valid), double(energy(valid)), ...
    binCount, sampleCount);
end

function validateMode(mode)
if ~isstruct(mode) || ~all(isfield(mode, {'c', 'omega'})) || ...
        ~isnumeric(mode.c) || isempty(mode.c) || ~ismatrix(mode.c) || ...
        ~all(isfinite(mode.c(:))) || ...
        ~isnumeric(mode.omega) || numel(mode.omega) ~= size(mode.c, 2)
    error('rvmdhilbert:InvalidMode', ...
        'Mode must be a finite RVMD mode struct with c and omega fields.');
end
end

function indices = validateModeIndices(indices, modeCount)
if isempty(indices)
    indices = 1:modeCount;
    return
end
if ~isnumeric(indices) || ~isreal(indices) || ~isvector(indices) || ...
        ~all(isfinite(indices(:))) || any(indices(:) < 1) || ...
        any(indices(:) > modeCount) || any(indices(:) ~= fix(indices(:))) || ...
        numel(unique(indices(:))) ~= numel(indices)
    error('rvmdhilbert:InvalidModeIndices', ...
        'ModeIndices must contain unique valid mode indices.');
end
indices = double(indices(:).');
end

function value = validatePositiveScalar(value, name)
if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
        ~isfinite(value) || value <= 0
    error(['rvmdhilbert:Invalid', name], ...
        '%s must be a finite positive real scalar.', name);
end
value = double(value);
end
