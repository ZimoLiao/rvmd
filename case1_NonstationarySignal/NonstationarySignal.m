% 1D nonstationary signal decomposition (RVMD reduces to VMD in this case)
clear
close all

addpath '../'

%% signal generation
T = 1000;
fs = 1 / T;
t = (1:T) / T;
% center frequencies of components
f_1 = 2;
f_2 = 24;
f_3 = 288;
% components
v_1 = (cos(2 * pi * f_1 * t));
v_2 = 1/4 * (cos(2 * pi * f_2 * t)) .* (10 * t.^2) .* [ones(1, ceil(T / 2)), zeros(1, floor(T / 2))];
v_3 = 1/16 * (cos(2 * pi * f_3 * t));
noiseAmp = 1/24;
% composite signal, including noise
f = v_1 + v_2 + v_3 + noiseAmp * randn(size(v_1));
f_spec = fft(f);

% plot signal
figure;
subplot(5, 1, 1:2);
plot(t, f, '-k');
xticklabels([])
title('Original data', 'Interpreter', 'latex')
subplot(5, 1, 3);
plot(t, v_1, '-k');
title(['component 1\quad$\omega_1=', num2str(f_1 * fs), '$'], 'Interpreter', 'latex')
subplot(5, 1, 4);
plot(t, v_2, '-k');
title(['component 2\quad$\omega_2=', num2str(f_2 * fs), '$'], 'Interpreter', 'latex')
subplot(5, 1, 5);
plot(t, v_3, '-k');
title(['component 3\quad$\omega_3=', num2str(f_3 * fs), '$'], 'Interpreter', 'latex')
xlabel('time','Interpreter', 'latex')
set(gca, 'Box', 'on');

%% parameters for RVMD
K = 3;              % number of modes
alpha = 1000;       % filtering parameter
tol = 1e-5;         % tolerance
N = 1000;           % maximum steps
init = 1;           % frequency initialization (1: uniformly distributed)
initFreqMax = 0.3;  % frequency initialization

%% computation
tic
[mode, info] = ...
    rvmd(f, K, alpha, 'Tolerance', tol, 'MaximumSteps', N, ...
    'InitFreqType', init, 'InitFreqMaximum', initFreqMax, ...
    'Device', 'cpu', 'FPPrecision', 'double');
toc

%% post-processing
% RVMD results (time-evolution coefficients)
f_rec = mode.phi * mode.c.';
figure;
subplot(K + 2, 1, 1:2);
plot(f_rec', '-b');
xticklabels([])
title('RVMD reconstructed data', ...
    'Interpreter', 'latex');

for k = 1:K
    subplot(K + 2, 1, k + 2); hold on;
    plot(t,mode.c(:, k), '-b');
    title(['mode ', num2str(k), '\quad$\omega_', num2str(k), '$=', ...
        num2str(mode.omega(k))], 'Interpreter', 'latex')
    if (k ~= K)
        xticklabels([])
    else
        xlabel('time','Interpreter', 'latex')
    end
    set(gca, 'Box', 'on');
end

% convergence curve
figure;
plot(info.Iteration.omega.')
xlabel('iteration step $n$', 'Interpreter', 'latex');
ylabel('$\omega_k^n$', 'Interpreter', 'latex')

% spectrum
figure; hold on;
f_spec = sum(abs(fft(f,[],2)).^2,1);
f_spec_rec = sum(abs(fft(f_rec,[],2)).^2,1);
nt = length(f_spec);
nk = ceil(nt/2);
freq = (0:(nk-1))/nt;
c_spec = abs(fft(mode.c,[],1)).^2;
plot(freq,f_spec(1:nk),'LineStyle','-','Color','b')
plot(freq,f_spec_rec(1:nk),'LineStyle','--','Color','m')
plot(freq,c_spec(1:nk,:))
legend('original','reconstructed','mode 1','mode 2','mode 3')
set(gca,'XScale','log','YScale','log')