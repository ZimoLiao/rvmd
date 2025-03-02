function [freq,E,spec] = frequencyspectrum(q,fs)

T = size(q,2);
spec = fft(q,[],2);
E = sum(abs(spec).^2,1);
freq = ifftshift((0:(T-1))-ceil(T/2))/T*fs;

end