function [qm] = mirrorextension(q)

t = size(q,2);
ts = floor(t/2);
qm = circshift([q,fliplr(q)],ts,2);

end