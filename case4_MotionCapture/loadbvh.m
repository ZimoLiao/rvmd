function [skeleton,time] = loadbvh(fname,varargin)
%% LOADBVH  Load a .bvh (Biovision) file.
%
% Loads BVH file specified by FNAME (with or without .bvh extension)
% and parses the file, calculating joint kinematics and storing the
% output in SKELETON.
%
% Optional argument 'delim' allows for setting the delimiter between
% fields. E.g., for a tab-separated BVH file:
%
%    skeleton = loadbvh('louise.bvh','delim','\t')
%
% By default 'delim' is set to the space character.
%
% SKELETON is a structure with N elements, where N is the number of joints.
% Each element consists of the following fields:
%
%     * `name` -- human-readable description of the joint
%     * `nestdepth` -- how many joints away from the origin
%     * `parent` -- index to parent joint
%     * `offset` -- translation from previous joint to current joint
%     * `Nchannels` -- number of channels describing pose
%                      (usually 3 for rotations only or 6 for 6DOF pose)
%     * `Nframes` -- number of time samples
%     * `order` -- the Euler angle order; e.g. [3 1 2] = ZXY
%     * `Dxyz` -- XYZ displacements (directly from `trans` matrix)
%     * `rxyz` -- XYZ rotations (to calculate `trans` matrix)
%     * `trans` -- transformation matrix 
%
% Some details on the BVH file structure are given in "Motion Capture File
% Formats Explained": http://www.dcs.shef.ac.uk/intranet/research/resmes/CS0111.pdf
% But most of it is fairly self-evident.

%% Options

p = inputParser;
p.addParameter('delim',' ');
parse(p,varargin{:});
opt = p.Results;

%% Load and parse header data
%
% The file is opened for reading, primarily to extract the header data (see
% next section). However, I don't know how to ask Matlab to read only up
% until the line "MOTION", so we're being a bit inefficient here and
% loading the entire file into memory. Oh well.

% add a file extension if necessary:
if ~strncmpi(fliplr(fname),'hvb.',4)
  fname = [fname,'.bvh'];
end

fid = fopen(fname);
if fid == -1
  error(['File "',fname,'" not found.'])
end
C = textscan(fid,'%s');
fclose(fid);
C = C{1};


%% Parse data
%
% This is a cheap tokeniser, not particularly clever.
% Iterate word-by-word, counting braces and extracting data.

% Initialise:
skeleton = [];
ii = 1;
nn = 0;
brace_count = 1;

while ~strcmp( C{ii} , 'MOTION' )
  
  ii = ii+1;
  token = C{ii};
  
  if strcmp( token , '{' )
    
    brace_count = brace_count + 1;
    
  elseif strcmp( token , '}' )
    
    brace_count = brace_count - 1;
    
  elseif strcmp( token , 'OFFSET' )
    
    skeleton(nn).offset = [str2double(C(ii+1)) ; str2double(C(ii+2)) ; str2double(C(ii+3))];
    ii = ii+3;
    
  elseif strcmp( token , 'CHANNELS' )
    
    skeleton(nn).Nchannels = str2double(C(ii+1));
    
    % The 'order' field is an index corresponding to the order of 'X' 'Y' 'Z'.
    % Subtract 87 because char numbers "X" == 88, "Y" == 89, "Z" == 90.
    if skeleton(nn).Nchannels == 3
      skeleton(nn).order = [C{ii+2}(1),C{ii+3}(1),C{ii+4}(1)]-87;
    elseif skeleton(nn).Nchannels == 6
      skeleton(nn).order = [C{ii+5}(1),C{ii+6}(1),C{ii+7}(1)]-87;
    else
      error('Not sure how to handle not (3 or 6) number of channels.')
    end
    
    if ~all(sort(skeleton(nn).order)==[1 2 3])
      error('Cannot read channels order correctly. Should be some permutation of [''X'' ''Y'' ''Z''].')
    end

    ii = ii + skeleton(nn).Nchannels + 1;

  elseif strcmp( token , 'JOINT' ) || strcmp( token , 'ROOT' )
    % Regular joint
    
    nn = nn+1;
    
    skeleton(nn).name = C{ii+1};
    skeleton(nn).nestdepth = brace_count;

    if brace_count == 1
      % root node
      skeleton(nn).parent = 0;
    elseif skeleton(nn-1).nestdepth + 1 == brace_count;
      % if I am a child, the previous node is my parent:
      skeleton(nn).parent = nn-1;
    else
      % if not, what is the node corresponding to this brace count?
      prev_parent = skeleton(nn-1).parent;
      while skeleton(prev_parent).nestdepth+1 ~= brace_count
        prev_parent = skeleton(prev_parent).parent;
      end
      skeleton(nn).parent = prev_parent;
    end
    
    ii = ii+1;
            
  elseif strcmp( [C{ii},' ',C{ii+1}] , 'End Site' )
    % End effector; unnamed terminating joint
    %
    % N.B. The "two word" token here is why we don't use a switch statement
    % for this code.
    
    nn = nn+1;
    
    skeleton(nn).name = ' ';
    skeleton(nn).offset = [str2double(C(ii+4)) ; str2double(C(ii+5)) ; str2double(C(ii+6))];
    skeleton(nn).parent = nn-1; % always the direct child
    skeleton(nn).nestdepth = brace_count;
    skeleton(nn).Nchannels = 0;
        
  end
  
end

%% Initial processing and error checking

Nnodes = numel(skeleton);
Nchannels = sum([skeleton.Nchannels]);
Nchainends = sum([skeleton.Nchannels]==0);

% Calculate number of header lines:
%  - 5 lines per joint
%  - 4 lines per chain end
%  - 2 additional lines ('HIERARCHY' and 'MOTION')
Nheaderlines = (Nnodes-Nchainends)*5 + Nchainends*4 + 2;

fid = fopen(fname);
for ii = 1:Nheaderlines
  tline = fgetl(fid);
end
if ~strcmp(tline,'MOTION')
   error('Could not parse BVH file %s. Number of header lines (before "MOTION" appears incorrect.',fname)
end

% get next non-blank line
tline = char([]);
while isempty(tline)
  tline = fgetl(fid);
  Nheaderlines = Nheaderlines + 1;
end
Nframes = sscanf(tline,'Frames: %f');

% get next non-blank line
tline = char([]);
while isempty(tline)
  tline = fgetl(fid);
  Nheaderlines = Nheaderlines + 1;
end
frame_time = sscanf(tline,'Frame Time: %f');

first_line = fgetl(fid);
fclose(fid);

% get the full data array from here:
rawdata = importdata(fname,opt.delim,Nheaderlines);

if ~isstruct(rawdata)
   error('Could not parse BVH file %s. Check the delimiter. First line of data appears to be:\n\n> %s',fname,first_line)
end

time = frame_time*(0:Nframes-1);

if size(rawdata.data,2) ~= Nchannels
  error('Error reading BVH file: channels count does not match.')
end

if size(rawdata.data,1) ~= Nframes
  warning('LOADBVH:frames_wrong','Error reading BVH file: frames count does not match; continuing anyway.')
  Nframes = size(rawdata.data,1);
end

%% Load motion data into skeleton structure
%
% We have three possibilities for each node we come across:
% (a) a root node that has displacements already defined,
%     for which the transformation matrix can be directly calculated;
% (b) a joint node, for which the transformation matrix must be calculated
%     from the previous points in the chain; and
% (c) an end effector, which only has displacement to calculate from the
%     previous node's transformation matrix and the offset of the end
%     joint.
%
% These are indicated in the skeleton structure, respectively, by having
% six, three, and zero "channels" of data.
% In this section of the code, the channels are read in where appropriate
% and the relevant arrays are pre-initialised for the subsequent calcs.

channel_count = 0;

for nn = 1:Nnodes
    
  if skeleton(nn).Nchannels == 6 % root node
    
    % assume translational data is always ordered XYZ
    skeleton(nn).Dxyz = repmat(skeleton(nn).offset,[1 Nframes]) + rawdata.data(:,channel_count+[1 2 3])';
    skeleton(nn).rxyz(skeleton(nn).order,:) = rawdata.data(:,channel_count+[4 5 6])';
        
    % Kinematics of the root element:
    skeleton(nn).trans = nan(4,4,Nframes);
    for ff = 1:Nframes
      skeleton(nn).trans(:,:,ff) = transformation_matrix(skeleton(nn).Dxyz(:,ff) , skeleton(nn).rxyz(:,ff) , skeleton(nn).order);
    end
    
  elseif skeleton(nn).Nchannels == 3 % joint node
        
    skeleton(nn).rxyz(skeleton(nn).order,:) = rawdata.data(:,channel_count+[1 2 3])';
    skeleton(nn).Dxyz  = nan(3,Nframes);
    skeleton(nn).trans = nan(4,4,Nframes);
    
  elseif skeleton(nn).Nchannels == 0 % end node
    skeleton(nn).Dxyz  = nan(3,Nframes);
  end
  
  channel_count = channel_count + skeleton(nn).Nchannels;
  skeleton(nn).Nframes = Nframes;
  
end


%% Calculate kinematics
%
% No calculations are required for the root nodes.

% For each joint, calculate the transformation matrix and for convenience
% extract each position in a separate vector.
for nn = find([skeleton.parent] ~= 0 & [skeleton.Nchannels] ~= 0)
  
  parent = skeleton(nn).parent;
  
  for ff = 1:Nframes
    transM = transformation_matrix( skeleton(nn).offset , skeleton(nn).rxyz(:,ff) , skeleton(nn).order );
    skeleton(nn).trans(:,:,ff) = skeleton(parent).trans(:,:,ff) * transM;
    skeleton(nn).Dxyz(:,ff) = skeleton(nn).trans([1 2 3],4,ff);
  end

end

% For an end effector we don't have rotation data;
% just need to calculate the final position.
for nn = find([skeleton.Nchannels] == 0)
  
  parent = skeleton(nn).parent;
  
  for ff = 1:Nframes
    transM = skeleton(parent).trans(:,:,ff) * [eye(3), skeleton(nn).offset; 0 0 0 1];
    skeleton(nn).Dxyz(:,ff) = transM([1 2 3],4);
  end

end

end



function transM = transformation_matrix(displ,rxyz,order)
% Constructs the transformation for given displacement, DISPL, and
% rotations RXYZ. The vector RYXZ is of length three corresponding to
% rotations around the X, Y, Z axes.
%
% The third input, ORDER, is a vector indicating which order to apply
% the planar rotations. E.g., [3 1 2] refers applying rotations RXYZ
% around Z first, then X, then Y.
%
% Years ago we benchmarked that multiplying the separate rotation matrices
% was more efficient than pre-calculating the final rotation matrix
% symbolically, so we don't "optimise" by having a hard-coded rotation
% matrix for, say, 'ZXY' which seems more common in BVH files.
% Should revisit this assumption one day.
%
% Precalculating the cosines and sines saves around 38% in execution time.

c = cosd(rxyz);
s = sind(rxyz);

RxRyRz(:,:,1) = [1 0 0; 0 c(1) -s(1); 0 s(1) c(1)];
RxRyRz(:,:,2) = [c(2) 0 s(2); 0 1 0; -s(2) 0 c(2)];
RxRyRz(:,:,3) = [c(3) -s(3) 0; s(3) c(3) 0; 0 0 1];

rotM = RxRyRz(:,:,order(1))*RxRyRz(:,:,order(2))*RxRyRz(:,:,order(3));

transM = [rotM, displ; 0 0 0 1];

end
