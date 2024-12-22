clear 
close all

load('shift_mode.mat','mode');
c(:,1) = mode.c;

load('votex_shedding_mode.mat','mode');
c(:,2) = real(mode.c(:,1));
c(:,3) = imag(mode.c(:,1));
c(:,4) = real(mode.c(:,2));
c(:,5) = imag(mode.c(:,2));

c = c(51:500,:);

%% 
figure; hold on;
plot(c);

%%
save('c.mat','c');