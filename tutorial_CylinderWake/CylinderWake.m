clear
close all

addpath '../'

%% load data
load('../case3_CylinderWake/data_cylinder.mat');

q = reshape(velocity,[S,T]);
q_base = q(:,1); % approximation
q_fluc = bsxfun(@minus,q,q_base);
velocity_mean = reshape(q_base, [V,I,J]);

%% spectrum
Ts = 1/4;
fs = 1/Ts;
[freq,q_E] = frequencyspectrum(q_fluc,fs); % TODO: 注意，一般情况下，需要加权！！！！
freqs = fftshift(freq);
q_Es = fftshift(q_E);

figure; hold on;
plot(freqs,q_Es,'-k','LineWidth',1.0);
xlim([0,floor(T/2-1)/T*fs]);
set(gca,'LineWidth',0.8,'Box','on','XScale','log','YScale','log')
set(gcf,'Color','white','Position',[100,450,480,240])

figure; hold on;
plot(freqs,q_Es,'-k','LineWidth',1.0);
xlim([0,floor(T/2-1)/T*fs]);
set(gca,'LineWidth',0.8,'Box','on','XScale','linear','YScale','log')
set(gcf,'Color','white','Position',[100,100,480,240])

%% shift mode
% filter bandwidth (prior)
b1 = 0.02*4;
alpha1 = (2*sqrt(2)-2)/(b1/fs)^2;

% parameters for RVMD
K = 1;                  % number of modes
alpha = alpha1;%1000;   % filtering parameter
tol = 5e-3;             % tolerance
N = 500;                % maximum steps
init = 0;               % frequency initialization (1: uniformly distributed)
initFreqMax = 0.2;      % frequency initialization
Device = 'gpu';         % device on which computation is performed
FPPrecision = 'single'; % float-pointing precision
nDC = 1;                % number of DC components

% computation
tic
[mode, info] = ...
    rvmd(q_fluc, K, alpha, 'Tolerance', tol, 'MaximumSteps', ceil(N/2), ...
    'InitFreqType', init, 'InitFreqMaximum', initFreqMax, 'Device', Device, ...
    'FPPrecision', FPPrecision, 'nDC', nDC);
toc

save('shift_mode.mat','mode','info');

% rms bandwidth (posterior)
c1 = mode.c; Ec1 = trapz(c1.^2);
dc1 = (c1(3:end)-c1(1:end-2))/2; dc1 = [c1(2)-c1(1);dc1;c1(end)-c1(end-1)]/Ts;
B1 = sqrt(trapz(dc1.^2)/Ec1*Ts);

% visualization
q_rec_1 = mode.phi*mode.c.';
[~,q_rec_1_E] = frequencyspectrum(q_rec_1,fs);
q_rec_1_Es = fftshift(q_rec_1_E);

figure; hold on;
plot(freqs,q_Es,'-k','LineWidth',1.0);
plot(freqs,q_rec_1_Es,'--r','LineWidth',1.0);
xline(b1/2,'-.r')
xline(B1/2,'-.b')
xlim([0,floor(T/2-1)/T*fs]);
set(gca,'LineWidth',0.8,'Box','on','XScale','log','YScale','log')
set(gcf,'Color','white','Position',[600,450,480,240])

figure; hold on;
plot(freqs,q_Es,'-k','LineWidth',1.0);
plot(freqs,q_rec_1_Es,'--r','LineWidth',1.0);
xline(b1/2,'-.k')
xline(B1/2,'-.r')
xlim([0,floor(T/2-1)/T*fs]);
set(gca,'LineWidth',0.8,'Box','on','XScale','linear','YScale','log')
set(gcf,'Color','white','Position',[600,100,480,240])

%% vortex shedding mode
% filter bandwidth
bk = B1*4;
alpha2 = (2*sqrt(2)-2)/((bk/fs).^2);

%% oscillating mode
q_fluc = q_fluc - mode.phi*mode.c.';
% q_fluc = hilbert(q_fluc.').'; % TODO: analytic representation of data

% parameters for RVMD
K = 4;                  % number of modes
alpha = alpha2; %1000;  % filtering parameter
tol = 5e-3;             % tolerance
N = 500;                % maximum steps
init = 1;               % frequency initialization (1: uniformly distributed)
initFreqMax = 0.05;     % frequency initialization
Device = 'gpu';         % device on which computation is performed
FPPrecision = 'single'; % float-pointing precision
nDC = 0;                % number of DC components

tic
[mode, info] = ...
    rvmd(q_fluc, K, alpha, 'Tolerance', tol, 'MaximumSteps', ceil(N/2), ...
    'InitFreqType', init, 'InitFreqMaximum', initFreqMax, 'Device', Device, ...
    'FPPrecision', FPPrecision, 'nDC', nDC);
toc

save('vortex_shedding_mode.mat','mode','info');

% visualization
for k = 1:K
    q_rec_k = mode.phi(:,k)*mode.c(:,k).';
    [~,q_rec_k_E] = frequencyspectrum(q_rec_k,fs);
    q_rec_k_Es = fftshift(q_rec_k_E);
    
    % rms bandwidth (posterior)
    ck = mode.c(:,k); Eck = trapz(ck.^2);
    ck_E = fftshift(abs(fft(ck)).^2); % TODO: 可以简化流程
    omega = ((0:T-1)-ceil(T/2))/T*fs;
    Bk = sqrt(trapz((abs(omega)-mode.omega(k)*fs).^2.*ck_E.')/Eck*fs/T);
    
    figure; hold on;
    plot(freqs,q_Es,'-k','LineWidth',1.0);
    plot(freqs,q_rec_k_Es,'--r','LineWidth',1.0);
    xline(mode.omega(k)*fs + bk/2,'-.r')
    xline(mode.omega(k)*fs - bk/2,'-.r')
    xline(mode.omega(k)*fs + Bk/2,'-.b')
    xline(mode.omega(k)*fs - Bk/2,'-.b')
    xlim([0,floor(T/2-1)/T*fs]);
    ylim([1e2,1e8]);
    set(gca,'LineWidth',0.8,'Box','on','XScale','log','YScale','log')
    set(gcf,'Color','white','Position',[600,450,480,240])
    
end

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
%
% % velocity field
% for k = 1:K
%     figure;
%     subplot(4,1,1);
%     pcolor(gridx, gridy, imag(squeeze(phi_velocity(1,:,:,k)))); colorbar
%     shading interp; axis equal tight;
%     title('$u$','Interpreter','latex')
%     subplot(4,1,2);
%     pcolor(gridx, gridy, imag(squeeze(phi_velocity(2,:,:,k)))); colorbar
%     shading interp; axis equal tight;
%     title('$v$','Interpreter','latex')
%     subplot(4,1,3);
%     pcolor(xx, yy, imag(squeeze(phi_vorticity(:,:,k)))); colorbar
%     shading interp; axis equal tight;
%     title('$\omega_z$','Interpreter','latex')
%     subplot(4,1,4); hold on;
%     plot(imag(mode.c(:,k)))
%     ylabel('$c(t)$','Interpreter','latex')
%     title(['RVMD mode ',num2str(k)],'Interpreter','latex')
%     set(gcf,'Position',[600,200,500,800])
% end