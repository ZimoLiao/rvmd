% The transient cylinder wakes
clear
close all

addpath '../'

%% load data
load('data_cylinder.mat');

q = reshape(velocity,[S,T]);
q_mean = mean(q,2);
q_fluc = bsxfun(@minus,q,q_mean);
velocity_mean = reshape(q_mean, [V,I,J]);

%% parameters for RVMD
K = 5;          	% number of modes
alpha = 1000;       % filtering parameter
tol = 5e-3;         % tolerance
N = 500;            % maximum steps
init = 1;           % frequency initialization (1: uniformly distributed)
initFreqMax = 0.05; % frequency initialization

%% computation
tic
[mode, info] = ...
    rvmd(q_fluc, K, alpha, 'Tolerance', tol, 'MaximumSteps', N, ...
    'InitFreqType', init, 'InitFreqMaximum', initFreqMax, 'Device', 'gpu', ...
    'FPPrecision', 'single');
toc

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
    pcolor(gridx, gridy, squeeze(phi_velocity(1,:,:,k)));
    shading interp; axis equal tight;
    title('$u$','Interpreter','latex')
    subplot(4,1,2);
    pcolor(gridx, gridy, squeeze(phi_velocity(2,:,:,k)));
    shading interp; axis equal tight;
    title('$v$','Interpreter','latex')
    subplot(4,1,3);
    pcolor(xx, yy, squeeze(phi_vorticity(:,:,k)));
    shading interp; axis equal tight;
    title('$\omega_z$','Interpreter','latex')
    subplot(4,1,4);
    plot(mode.c(:,k))
    ylabel('$c(t)$','Interpreter','latex')
    title(['RVMD mode ',num2str(k)],'Interpreter','latex')
    set(gcf,'Position',[200,200,500,800])
end

%% reconstruction
% mean flow + shift-mode
for t = 1:T
    shift_mode_vel_rec(:,:,:,t) = squeeze(phi_velocity(:,:,:,1))*mode.c(t,1);
    shift_mode_vor_rec(:,:,t) = squeeze(phi_vorticity(:,:,1))*mode.c(t,1);
end
shift_mode_vel_rec = bsxfun(@plus, shift_mode_vel_rec, velocity_mean);
shift_mode_vor_rec = bsxfun(@plus,shift_mode_vor_rec, vorticity_mean);

% mean flow + shift-mode + intermediate vortex shedding modes
for t = 1:T
    inter_mode_vel_rec(:,:,:,t) = squeeze(phi_velocity(:,:,:,2))*mode.c(t,2)+ ...
        squeeze(phi_velocity(:,:,:,3))*mode.c(t,3)+ ...
        squeeze(phi_velocity(:,:,:,1))*mode.c(t,1);
    inter_mode_vor_rec(:,:,t) = squeeze(phi_vorticity(:,:,2))*mode.c(t,2)+ ...
        squeeze(phi_vorticity(:,:,3))*mode.c(t,3)+ ...
        squeeze(phi_vorticity(:,:,1))*mode.c(t,1);
end
inter_mode_vel_rec = bsxfun(@plus, inter_mode_vel_rec, velocity_mean);
inter_mode_vor_rec = bsxfun(@plus,inter_mode_vor_rec, vorticity_mean);

% mean flow + shift-mode + post-transient vortex shedding modes
for t = 1:T
    post_mode_vel_rec(:,:,:,t) = squeeze(phi_velocity(:,:,:,4))*mode.c(t,4)+ ...
        squeeze(phi_velocity(:,:,:,5))*mode.c(t,5)+ ...
        squeeze(phi_velocity(:,:,:,1))*mode.c(t,1);
    post_mode_vor_rec(:,:,t) = squeeze(phi_vorticity(:,:,4))*mode.c(t,4)+ ...
        squeeze(phi_vorticity(:,:,5))*mode.c(t,5)+ ...
        squeeze(phi_vorticity(:,:,1))*mode.c(t,1);
end
post_mode_vel_rec = bsxfun(@plus, post_mode_vel_rec, velocity_mean);
post_mode_vor_rec = bsxfun(@plus,post_mode_vor_rec, vorticity_mean);

% animation
% figure;
% for t = 1:5:T
%     subplot(3, 1, 1);
%     pcolor(xx, yy, squeeze(shift_mode_vor_rec(:,:,t))); hold on;
%     %     quiver(gridx, gridy, squeeze(shift_mode_vel_rec(1,:,:,t)), ...
%     %         squeeze(inter_mode_vel_rec(2,:,:,t)),'Color','k');
%     hold off;
%     shading interp; axis equal tight;
%     caxis([-3,3])
%     title(['t = ',num2str(t),' (',num2str(T),') shift-mode'])
%     
%     subplot(3, 1, 2);
%     pcolor(xx, yy, squeeze(inter_mode_vor_rec(:,:,t))); hold on;
%     %     quiver(gridx, gridy, squeeze(inter_mode_vel_rec(1,:,:,t)), ...
%     %         squeeze(inter_mode_vel_rec(2,:,:,t)),'Color','k');
%     hold off;
%     shading interp; axis equal tight;
%     caxis([-3,3])
%     title('intermediate vortex shedding modes')
%     
%     subplot(3, 1, 3);
%     pcolor(xx, yy, squeeze(post_mode_vor_rec(:,:,t))); hold on;
%     %     quiver(gridx, gridy, squeeze(post_mode_vel_rec(1,:,:,t)), ...
%     %         squeeze(inter_mode_vel_rec(2,:,:,t)),'Color','k');
%     hold off;
%     shading interp; axis equal tight;
%     caxis([-3,3])
%     set(gcf,'Position',[200,200,500,800])
%     title('post-transient vortex shedding modes')
%     
%     pause(0.01)
% end