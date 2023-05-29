% Lorenz system
clear
close all

%% load data
load('data/data_lorenz.mat');
% data from Brunton, S. L., Brunton, B. W., Proctor, J. L., Kaiser, E.,
%   & Kutz, J. N. (2017). Chaos as an intermittently forced linear system.
%   Nature communications, 8(1), 19.

x = x(1:10:end,:).';

% plot spectrum
figure;
pspectrum(x.');
set(gca,'XScale','log')

%% parameters for RVMD
K = 5;              % number of modes
alpha = 3000;       % filtering parameter
tol = 1e-4;         % tolerance
N = 2000;           % maximum steps
init = 1;           % frequency initialization (1: uniformly distributed)
initFreqMax = 0.2;  % frequency initialization

%% computation
[phi, c, omega, energy, omega_iter, diff_iter] = ...
    rvmd(x, K, alpha, tol, N, init, initFreqMax);

%% post-processing
% energy rank
[~,rank] = maxk(energy,K);

% RVMD results (time-evolution coefficients)
x_rec = phi * c.';
figure;
subplot(K + 2, 1, 1:2);
plot(x_rec', '-b');
xticklabels([])
title('RVMD reconstructed data', ...
    'Interpreter', 'latex');

for k = 1:K
    subplot(K + 2, 1, k + 2); hold on;
    plot(c(:, k), '-b');
    title(['mode ', num2str(k), '\quad$\omega_', num2str(k), '$=', ...
        num2str(omega(k)),'\quad energy rank ',num2str(rank(k)), ...
        ' $E_', num2str(k), '$=',num2str(energy(k))], 'Interpreter', 'latex')
    if (k ~= K)
        xticklabels([])
    else
        xlabel('time','Interpreter', 'latex')
    end
    ylim([-35,35])
    set(gca, 'Box', 'on');
end

% mode
figure;
for k = 1:K
    subplot(1, K, k); hold on;
    plot(phi(:, k), '-b');
    title(['mode ', num2str(k), '\quad$\omega_', num2str(k), '$=', ...
        num2str(omega(k))], 'Interpreter', 'latex')
    if (k ~= K)
        xticklabels([])
    else
        xlabel('x','Interpreter', 'latex')
    end
    ylim([-1,1])
    set(gca, 'Box', 'on');
end

% convergence curve
figure;
plot(omega_iter.')
xlabel('iteration step $n$', 'Interpreter', 'latex');
ylabel('$\omega_k^n$', 'Interpreter', 'latex')

% spectrum
figure; hold on;
pspectrum(c);
set(gca,'XScale','log')

%% recontruction
for k = 1:K
    figure;
    for i = 1:3
        subplot(3, 1, i); hold on;
        plot(x(i,:), '-b');
        plot(x_rec(i,:), '--m');
        %         if k ~=1 && k~=2
        %             plot(phi(i,k)*c(:,k)+phi(i,1)*c(:,1)+phi(i,2)*c(:,2), '-k');
        %         else
        plot(phi(i,k)*c(:,k), '-k');
        %         end
    end
    title(['mode ', num2str(k), '\quad$\omega_', num2str(k), '$=', ...
        num2str(omega(k)),'\quad energy rank ',num2str(rank(k)), ...
        ' $E_', num2str(k), '$=',num2str(energy(k))], 'Interpreter', 'latex')
end