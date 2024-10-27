clear
close all

t = 1:1000;
Q = cos(t*pi/100)+1i*sin(t*pi/50);
K = 2;
Alpha = 2000;

[mode, info] = vmd(Q, K, Alpha);

%% plot results
figure; hold on;
plot(real(Q),'Color','k')
for k = 1:K
    plot(real(mode.c(:,k)))
end
set(gcf,'Position',[100,100,480,400])

figure; hold on;
plot(imag(Q),'Color','k')
for k = 1:K
    plot(imag(mode.c(:,k)))
end
set(gcf,'Position',[600,100,480,400])