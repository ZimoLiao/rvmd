function [qim] = imirrorextension(q)

t = size(q,2);
ts = floor(t/4);
qim = circshift(q,-ts,2);
qim = qim(1:t/2);

end