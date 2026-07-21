% The transient cylinder wakes
clear
close all

addpath '../'

%% load data
load('data_cylinder.mat');

q = reshape(velocity,[S,T]);
q_base = q(:,1);
q_fluc = bsxfun(@minus,q,q_base); 
velocity_mean = reshape(q_base, [V,I,J]);

%% parameters for RVMD
K = 1;                  % number of modes
alpha = 1000;           % filtering parameter
tol = 5e-3;             % tolerance
N = 500;                % maximum steps
init = 0;               % frequency initialization (1: uniformly distributed)
initFreqMax = 0.2;      % frequency initialization
Device = 'gpu';         % device on which computation is performed
FPPrecision = 'single'; % floating-point precision
nDC = 1;                % number of DC components

%% computation
tic
[mode, info] = ...
    rvmd(q_fluc, K, alpha, 'Tolerance', tol, 'MaximumSteps', ceil(N/2), ...
    'InitFreqType', init, 'InitFreqMaximum', initFreqMax, 'Device', Device, ...
    'FPPrecision', FPPrecision, 'nDC', nDC);
toc

save('shift_mode.mat','mode','info');

% bandwidth
c1 = mode.c; Ec1 = trapz(c1.^2);
dc1 = (c1(3:end)-c1(1:end-2))/2; dc1 = [c1(2)-c1(1);dc1;c1(end)-c1(end-1)];
Bsquare = trapz(dc1.^2)/Ec1;
alpha = (2*sqrt(2)-2)/(Bsquare*4);

%% oscillating mode
q_fluc = q_fluc - mode.phi*mode.c.';
% q_fluc = hilbert(q_fluc.').'; % TODO: analytic representation of data

% parameters for RVMD
K = 4;                  % number of modes
alpha = alpha; %1000;         	% filtering parameter
tol = 5e-3;             % tolerance
N = 500;                % maximum steps
init = 1;               % frequency initialization (1: uniformly distributed)
initFreqMax = 0.05;      % frequency initialization
Device = 'gpu';         % device on which computation is performed
FPPrecision = 'single'; % floating-point precision
nDC = 0;                % number of DC components

tic
[mode, info] = ...
    rvmd(q_fluc, K, alpha, 'Tolerance', tol, 'MaximumSteps', ceil(N/2), ...
    'InitFreqType', init, 'InitFreqMaximum', initFreqMax, 'Device', Device, ...
    'FPPrecision', FPPrecision, 'nDC', nDC);
toc

save('vortex_shedding_mode.mat','mode','info');



%% post-processing
phi_velocity = reshape(mode.phi,[V,I,J,K]);
[xx,yy] = meshgrid(x,y);
vorticity_mean = curl(xx, yy, ...
    squeeze(velocity_mean(1,:,:)).', squeeze(velocity_mean(2,:,:)).');
for k = 1:K
    phi_vorticity(:,:,k) = curl(xx, yy, ...
        squeeze(phi_velocity(1,:,:,k)).', squeeze(phi_velocity(2,:,:,k)).');
end

% velocity field
for k = 1:K
    figure;
    subplot(4,1,1);
    pcolor(gridx, gridy, real(squeeze(phi_velocity(1,:,:,k)))); colorbar
    shading interp; axis equal tight;
    title('$u$','Interpreter','latex')
    subplot(4,1,2);
    pcolor(gridx, gridy, real(squeeze(phi_velocity(2,:,:,k)))); colorbar
    shading interp; axis equal tight;
    title('$v$','Interpreter','latex')
    subplot(4,1,3);
    pcolor(xx, yy, real(squeeze(phi_vorticity(:,:,k)))); colorbar
    shading interp; axis equal tight;
    title('$\omega_z$','Interpreter','latex')
    subplot(4,1,4); hold on;
    plot(real(mode.c(:,k)))
    ylabel('$c(t)$','Interpreter','latex')
    title(['RVMD mode ',num2str(k)],'Interpreter','latex')
    set(gcf,'Position',[100,200,500,800])
end

% velocity field
for k = 1:K
    figure;
    subplot(4,1,1);
    pcolor(gridx, gridy, imag(squeeze(phi_velocity(1,:,:,k)))); colorbar
    shading interp; axis equal tight;
    title('$u$','Interpreter','latex')
    subplot(4,1,2);
    pcolor(gridx, gridy, imag(squeeze(phi_velocity(2,:,:,k)))); colorbar
    shading interp; axis equal tight;
    title('$v$','Interpreter','latex')
    subplot(4,1,3);
    pcolor(xx, yy, imag(squeeze(phi_vorticity(:,:,k)))); colorbar
    shading interp; axis equal tight;
    title('$\omega_z$','Interpreter','latex')
    subplot(4,1,4); hold on;
    plot(imag(mode.c(:,k)))
    ylabel('$c(t)$','Interpreter','latex')
    title(['RVMD mode ',num2str(k)],'Interpreter','latex')
    set(gcf,'Position',[600,200,500,800])
end
