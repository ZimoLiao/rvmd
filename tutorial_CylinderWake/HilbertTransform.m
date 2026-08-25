clear
close all

addpath '../'
load('vortex_shedding_mode.mat');

sampleRate = 4;
modeIndex = 3;
analysis = rvmdhilbert(mode, sampleRate, ...
    'ModeIndices', modeIndex, 'MirrorExtension', false);
rvmdhilbertplot(analysis);

figure; hold on;
plot(analysis.Time, mode.c(:, modeIndex), '-k');
plot(analysis.Time, real(analysis.AnalyticSignal), '--r');
plot(analysis.Time, imag(analysis.AnalyticSignal), '--b');
xlabel('Time');
legend('coefficient', 'analytic real part', 'analytic imaginary part');
