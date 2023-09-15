% Mode decomposition for motion capture data
clear
close all

addpath '../'

%% load data
[skeleton,time] = loadbvh('49_03.bvh');
Njoints = numel(skeleton);

frame_first = 2;
frame_last = 1501;
S = Njoints*3;
T = frame_last-frame_first+1;

% convert bvh data into snapshots
[q, parent] = bvh2snapshot(skeleton,frame_first,frame_last);
q_mean = mean(q,2);
q_fluc = bsxfun(@minus, q, q_mean);

%% parameters for RVMD
K = 15;          	% number of modes
alpha = 2000;       % filtering parameter
tol = 1e-2;         % tolerance
N = 1500;            % maximum steps
init = 1;           % frequency initialization (1: uniformly distributed)
initFreqMax = 0.2;  % frequency initialization

%% computation
tic
[mode, info] = ...
    rvmd(q_fluc, K, alpha, 'Tolerance', tol, 'MaximumSteps', N, ...
    'InitFreqType', init, 'InitFreqMaximum', initFreqMax, ...
    'Device', 'gpu', 'FloatPrecision', 'single');
toc

%% post-processing
% reconstruction
rec_fluc = mode.phi*mode.c.';
rec = rec_fluc;
% rec = bsxfun(@plus,rec_fluc,q_mean);
for k = 1:K
    rec_mode(:,:,k) = mode.phi(:,k)*mode.c(:,k).';
end

err_fluc = norm(rec_fluc-q_fluc,'fro')/norm(q_fluc,'fro')
err = norm(bsxfun(@plus,rec_fluc,q_mean)-q,'fro')/norm(q,'fro')

%% plot spectrum
T_spec = floor(T/2+1);
freq_list = (0:(T_spec-1))/T;
q_flux_spec = sum(abs(fft(q_fluc,[],2)).^2,1);
q_flux_spec = q_flux_spec(1:T_spec);
rec_spec = sum(abs(fft(rec,[],2)).^2,1);
rec_spec = rec_spec(1:T_spec);
for k = 1:K
    temp = sum(abs(fft(rec_mode(:,:,k),[],2)).^2,1);
    rec_mode_spec(k,:) = temp(1:T_spec);
end

figure; hold on;
plot(freq_list, q_flux_spec, '--k', 'LineWidth',1.);
plot(freq_list, rec_spec, '-r', 'LineWidth',1.);
xlim([-inf,0.5])
ylim([1e2,1e11])
set(gca,'XScale','log','YScale','log');
legend({'data','reconstruction'})
xlabel('frequency')
ylabel('Spatially averaged spectrum')

%% RVMD mode: energy_k - omega_k
figure;
stem(mode.omega,mode.energy.^2);
set(gca,'XScale','log','YScale','log');
xlim([-inf,0.5])
ylim([1e2,1e11])

% figure;
% for t = 1:10:info.T
%     for nn = 1:Njoints
%         plot3(rec(3*nn-2,t),rec(3*nn-1,t),rec(3*nn,t),'.','MarkerSize',20); hold on;
%         plot3(q(3*nn-2,t),q(3*nn-1,t),q(3*nn,t),'.','MarkerSize',20,'Color','k'); hold on;
%     end
%     axis equal
%     pause(0.01);
%     hold off;
% end

%% mode
% k_low = [1:6];
% k_high = [11:K];
% 
% vidObj = VideoWriter('MotionCapture'); open(vidObj);
% 
% figure;
% for t = 1:2:info.T
%     subplot(4,3,[1,4,7])
%     for nn = 1:Njoints
%         plot3(q(3*nn,t),q(3*nn-2,t),q(3*nn-1,t),'.','MarkerSize',20,'Color','k'); hold on;
%         if parent(nn) ~= 0
%             plot3([q(3*parent(nn),t), q(3*nn,t)],...
%                 [q(3*parent(nn)-2,t), q(3*nn-2,t)],...
%                 [q(3*parent(nn)-1,t), q(3*nn-1,t)],'Color','k');
%         end
%         
%         zz = rec(3*nn,t)+q_mean(3*nn);
%         xx = rec(3*nn-2,t)+q_mean(3*nn-2);
%         yy = rec(3*nn-1,t)+q_mean(3*nn-1);
%         plot3(zz,xx,yy,'.','MarkerSize',20,'Color','m'); hold on;
%         if parent(nn) ~= 0
%             zzp = rec(3*parent(nn),t)+q_mean(3*parent(nn));
%             xxp = rec(3*parent(nn)-2,t)+q_mean(3*parent(nn)-2);
%             yyp = rec(3*parent(nn)-1,t)+q_mean(3*parent(nn)-1);
%             if parent(nn) ~= 0
%                 plot3([zzp, zz],...
%                     [xxp, xx],...
%                     [yyp, yy],'Color','m');
%             end
%         end
%     end
%     xlabel('Z')
%     ylabel('X')
%     zlabel('Y')
%     axis equal
%     xlim([0,14])
%     ylim([-12,12])
%     zlim([0,35])
%     v = [5 2 3];
%     [caz,cel] = view(v);
%     hold off;
%     title('Motion')
%     
%     subplot(4,3,[2,5,8])
%     for nn = 1:Njoints
%         zz = sum(rec_mode(3*nn,t,k_low),3)+q_mean(3*nn);
%         xx = sum(rec_mode(3*nn-2,t,k_low),3)+q_mean(3*nn-2);
%         yy = sum(rec_mode(3*nn-1,t,k_low),3)+q_mean(3*nn-1);
%         plot3(zz,xx,yy,'.','MarkerSize',20,'Color','r'); hold on;
%         if parent(nn) ~= 0
%             zzp = sum(rec_mode(3*parent(nn),t,k_low),3)+q_mean(3*parent(nn));
%             xxp = sum(rec_mode(3*parent(nn)-2,t,k_low),3)+q_mean(3*parent(nn)-2);
%             yyp = sum(rec_mode(3*parent(nn)-1,t,k_low),3)+q_mean(3*parent(nn)-1);
%             if parent(nn) ~= 0
%                 plot3([zzp, zz],...
%                     [xxp, xx],...
%                     [yyp, yy],'Color','r');
%             end
%         end
%     end
%     xlabel('Z')
%     ylabel('X')
%     zlabel('Y')
%     axis equal
%     xlim([0,14])
%     ylim([-12,12])
%     zlim([0,35])
%     v = [5 2 3];
%     [caz,cel] = view(v);
%     hold off;
%     title('Low frequency (Leg lifting modes)')
%     subplot(4,3,11)
%     plot(mode.c(:,k_low),'Color',[0.5,0.5,0.5]); hold on;
%     plot(mode.c(1:t,k_low)); hold off;
%     
%     subplot(4,3,[3,6,9])
%     for nn = 1:Njoints
%         zz = sum(rec_mode(3*nn,t,k_high),3)+q_mean(3*nn);
%         xx = sum(rec_mode(3*nn-2,t,k_high),3)+q_mean(3*nn-2);
%         yy = sum(rec_mode(3*nn-1,t,k_high),3)+q_mean(3*nn-1);
%         plot3(zz,xx,yy,'.','MarkerSize',20,'Color','b'); hold on;
%         if parent(nn) ~= 0
%             zzp = sum(rec_mode(3*parent(nn),t,k_high),3)+q_mean(3*parent(nn));
%             xxp = sum(rec_mode(3*parent(nn)-2,t,k_high),3)+q_mean(3*parent(nn)-2);
%             yyp = sum(rec_mode(3*parent(nn)-1,t,k_high),3)+q_mean(3*parent(nn)-1);
%             if parent(nn) ~= 0
%                 plot3([zzp, zz],...
%                     [xxp, xx],...
%                     [yyp, yy],'Color','b');
%             end
%         end
%     end
%     xlabel('Z')
%     ylabel('X')
%     zlabel('Y')
%     axis equal
%     xlim([0,14])
%     ylim([-12,12])
%     zlim([0,35])
%     v = [5 2 3];
%     [caz,cel] = view(v);
%     hold off;
%     title('High frequency (Jumping modes)')
%     subplot(4,3,12)
%     plot(mode.c(:,k_high),'Color',[0.5,0.5,0.5]); hold on;
%     plot(mode.c(1:t,k_high)); hold off;
%     
%     set(gcf,'Position',[100,100,1200,600],'Color','white')
%     %     pause(0.001);
%     
%     drawnow
%     writeVideo(vidObj,getframe(gcf));
% end
% 
% close(vidObj)