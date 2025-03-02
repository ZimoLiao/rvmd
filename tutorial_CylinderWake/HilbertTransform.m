clear
close all

addpath '../'

load('vortex_shedding_mode.mat');

Ts = 1/4;
fs = 1/Ts;

k = 3;

%%
T = length(mode.c(:,k));
times = (1:T)*Ts;
c = mode.c(:,k);
ca = hilbert(c);

%% 
[freq,c_E] = frequencyspectrum(c.',fs);
[~,ca_E] = frequencyspectrum(ca.',fs);

fmax = 0.5;
freqs = fftshift(freq);
c1s = fftshift(c_E);
c1as = fftshift(ca_E);

figure; hold on;
plot(freqs,c1s,'-k','LineWidth',1.0);
plot(freqs,[c1as(1:ceil(T/2))/4,c1as(ceil(T/2)+1),c1as(ceil(T/2)+2:end)/4],'--r','LineWidth',1.0);
xlim([-1,1]*fmax);
ylim([1e0,1e8]);
set(gca,'LineWidth',0.8,'Box','on','XScale','linear','YScale','log')
set(gcf,'Color','white','Position',[100,100,480,240])

%%
figure; hold on;
plot(times,c,'-k')
plot(times,real(ca),'--r')
plot(times,imag(ca),'--b')