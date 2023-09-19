% Lorenz Attractor
clear
close all

addpath '../'

%% load data
load('data_lorenz.mat');
% data from Brunton, S. L., Brunton, B. W., Proctor, J. L., Kaiser, E.,
%   & Kutz, J. N. (2017). Chaos as an intermittently forced linear system.
%   Nature communications, 8(1), 19.

x = x(1:10:end,:).';

%% parameters for RVMD
K = 5;              % number of modes
alpha = 1000;       % filtering parameter
tol = 1e-4;         % tolerance
N = 8000;           % maximum steps
init = 1;           % frequency initialization (1: uniformly distributed)
initFreqMax = 0.2;  % frequency initialization

%% computation
% [phi, c, omega, energy, omega_iter, diff_iter] = ...
%     rvmd(x, K, alpha, tol, N, init, initFreqMax);

tic
[mode, info] = ...
    rvmd(x, K, alpha, 'Tolerance', tol, 'MaximumSteps', N, 'InitFreqType', init, 'InitFreqMaximum', initFreqMax, ...
    'FPPrecision', 'single');
toc

%% post-processing
% energy rank
[~,index] = maxk(mode.energy,K);
for k = 1:K
    rank(index(k)) = k;
end

% time-evolution coefficients
x_rec = mode.phi * mode.c.';
figure;
subplot(K + 2, 1, 1:2);
plot(x_rec', '-k');
xticklabels([])
title('RVMD reconstructed data', ...
    'Interpreter', 'latex');

for k = 1:K
    subplot(K + 2, 1, k + 2); hold on;
    plot(mode.c(:, k), '-k');
    title(['mode ', num2str(k), '\quad$\omega_', num2str(k), '$=', ...
        num2str(mode.omega(k)),'\quad energy rank ',num2str(rank(k)), ...
        ' $E_', num2str(k), '$=',num2str(mode.energy(k))], 'Interpreter', 'latex')
    if (k ~= K)
        xticklabels([])
    else
        xlabel('time','Interpreter', 'latex')
    end
    ylim([-35,35])
    set(gca, 'Box', 'on');
end
set(gcf,'Color','white')

% modes
figure;
for k = 1:K
    subplot(1, K, k); hold on;
    plot(mode.phi(:, k), '-k');
    title(['mode ', num2str(k), '\quad$\omega_', num2str(k), '$=', ...
        num2str(mode.omega(k))], 'Interpreter', 'latex')
    xlabel('x','Interpreter', 'latex')
    ylim([-1,1])
    set(gca, 'Box', 'on');
end
set(gcf,'Color','white')

% spectrum
x_spec = sum(abs(fft(x,[],2)).^2,1);
x_spec_rec = sum(abs(fft(x_rec,[],2)).^2,1);
nt = length(x_spec);
nk = ceil(nt/2);
freq = (0:(nk-1))/nt;
c_spec = abs(fft(mode.c,[],1)).^2;

figure;
for k = 1:K
    subplot(K, 1, k); hold on;
    plot(freq,x_spec(1:nk),'LineStyle','-','Color','b')
    plot(freq,x_spec_rec(1:nk),'LineStyle','--','Color','m')
    plot(freq,c_spec(1:nk,k),'Color','k')
    set(gca,'XScale','log','YScale','log')
    title(['mode ', num2str(k), '\quad$\omega_', num2str(k), '$=', ...
        num2str(mode.omega(k))], 'Interpreter', 'latex')
    xlabel('x','Interpreter', 'latex')
end
set(gcf,'Color','white')

%% convergence curve
figure;
plot(info.Iteration.omega.')
xlabel('iteration step $n$', 'Interpreter', 'latex');
ylabel('$\omega_k^n$', 'Interpreter', 'latex')
set(gcf,'Color','white')

%% mode recontruction
% inividual mode
for k = 1:K
    figure;
    for i = 1:3
        subplot(3, 1, i); hold on;
        yline(8.485,'Color',[0.5,0.5,0.5])
        yline(-8.485,'Color',[0.5,0.5,0.5])
        plot(x(i,:), '-b');
        plot(x_rec(i,:), '--m');
        plot(mode.phi(i,k)*mode.c(:,k), '-k');
        if i==1
            title(['mode ', num2str(k), '\quad$\omega_', num2str(k), '$=', ...
                num2str(mode.omega(k)),'\quad energy rank ',num2str(rank(k)), ...
                ' $E_', num2str(k), '$=',num2str(mode.energy(k))], 'Interpreter', 'latex')
        end
    end
    set(gcf,'Color','white')
end

% all modes
figure; hold on;
plot3(x(1,:),x(2,:),x(3,:), '-b');
plot3(x_rec(1,:),x_rec(2,:),x_rec(3,:), '--m');
set(gcf,'Color','white')