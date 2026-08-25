function figureHandle = rvmdhilbertplot(analysis)
%RVMDHILBERTPLOT Plot an RVMD Hilbert spectral analysis result.
%   RVMDHILBERTPLOT(ANALYSIS) creates one figure from the struct returned by
%   RVMDHILBERT. It displays the log10 Hilbert energy spectrum, the marginal
%   spectrum, and the selected modes' instantaneous-frequency trajectories.
%   The time and frequency coordinates are taken from ANALYSIS without
%   additional scaling.
%
%   FIGUREHANDLE = RVMDHILBERTPLOT(ANALYSIS) returns the created figure
%   handle. Plotting expands ANALYSIS.HilbertSpectrum to a full matrix for
%   display.

required = {'Time', 'FrequencyBins', 'HilbertSpectrum', ...
    'MarginalSpectrum', 'InstantaneousFrequency', 'ModeIndices'};
if ~isstruct(analysis) || ~all(isfield(analysis, required))
    error('rvmdhilbertplot:InvalidAnalysis', ...
        'Input must be a struct returned by RVMDHILBERT.');
end

figureHandle = figure('Color', 'white');

subplot(2, 2, [1, 3]);
energy = full(analysis.HilbertSpectrum);
positiveEnergy = energy(energy > 0);
if isempty(positiveEnergy)
    logEnergy = zeros(size(energy));
else
    floorEnergy = max(min(positiveEnergy), max(positiveEnergy) * 1e-12);
    logEnergy = log10(max(energy, floorEnergy));
end
imagesc(analysis.Time, analysis.FrequencyBins, logEnergy);
axis xy;
xlabel('Time');
ylabel('Frequency');
title('Hilbert energy spectrum');
colorbar;

subplot(2, 2, 2);
plot(analysis.MarginalSpectrum, analysis.FrequencyBins, 'LineWidth', 1.2);
grid on;
xlabel('Integrated energy');
ylabel('Frequency');
title('Marginal spectrum');

subplot(2, 2, 4);
plot(analysis.Time, analysis.InstantaneousFrequency, 'LineWidth', 1.0);
grid on;
xlabel('Time');
ylabel('Frequency');
title('Instantaneous frequencies');
end
