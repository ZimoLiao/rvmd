function [q, parent] = bvh2snapshot(skeleton,frame_first,frame_last)

Njoints = numel(skeleton);
S = Njoints*3;
T = frame_last-frame_first+1;

% convert bvh data into snapshots
q = zeros(S, T);
for nn = 1:Njoints
    q((3*nn-2):(3*nn), :) = skeleton(nn).Dxyz(:,frame_first:frame_last);
end

parent = zeros(Njoints);
for nn = 1:Njoints
    parent(nn) = skeleton(nn).parent;
end

end

