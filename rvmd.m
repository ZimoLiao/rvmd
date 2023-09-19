function [mode, info] = rvmd(Q, K, Alpha, varargin)
% Reduced-order variational mode decomposition (RVMD)
%
% Description
%   [mode, info] = rvmd(Q, K, Alpha) return the RVMD results of the data
%   matrix Q whose first dimension is space, i.e., size(Q)=[S, T], where
%   S and T denote the numbers of sampling points in space and time,
%   respectively. K is the number of modes, and Alpha is the filtering
%   parameter. The columns of mode.phi are the spatial modes, and the
%   columns of mode.c are the time-evolution coefficients. mode.omega and
%   mode.energy are the central frequencies and mode energy, respectively.
%   The struct info contains information about the present decomposition,
%   including parameter setups and options.
%
%   [...] = rvmd(Q, K, Alpha, Name, Value) specifies the decomposition
%   options using one or more Name, Value pair arguments. Here is the list
%   of properties:
%   Properties:
%       'Tolerance' - tolerance value for the stopping criterion of RVMD
%       iteration
%       'MaximumSteps' - maximum number of iteration steps for the stopping
%       criterion of RVMD iteration
%       'InitFreqType' - specifies how to initialize the central frequencies
%       in domain [0, min(InitFreqMaximum, 0.5)]
%           Value	Strategy
%           -1    	randomly distributed
%       	0      	all zero
%           1     	uniformly distributed
%       'InitFreqMaximum' - specifies the domain of initial central
%       frequencies [0, min(InitFreqMaximum, 0.5)]
%       'Device' - device on which computation is performed: `cpu' or `gpu'
%       'FPPrecision' - float-pointing precision for RVMD computation:
%       `single' or `double'
%
% Reference
%   [1] Liao, Z.-M., Zhao, Z., Chen, L.-B., Wan, Z.-H., Liu, N.-S.
%       & Lu, X.-Y. 2023 Reduced-order variational mode decomposition
%       to reveal transient and non-stationary dynamics in fluid flows.
%       J. Fluid Mech., 966, A7. doi:10.1017/jfm.2023.435

% constants
FREQUNIFORM = 1;
FREQALLZERO = 0;
FREQRANDOM = -1;
FREQMAXIMUM = 0.5;

% default values/options
default_Tolerance = 5e-3;
default_MaximumSteps = 500;
default_InitFreqType = FREQUNIFORM;
default_InitFreqMaximum = FREQMAXIMUM;
default_Device = 'cpu';
default_FPPrecision = 'single';

% input parser
p = inputParser;
addRequired(p, 'Q');
validScalar = @(x) isnumeric(x) && (x>0);
addRequired(p, 'K', validScalar);
addRequired(p, 'alpha', validScalar);
addParameter(p, 'Tolerance', default_Tolerance);
addParameter(p, 'MaximumSteps', default_MaximumSteps);
addParameter(p, 'InitFreqType', default_InitFreqType);
addParameter(p, 'InitFreqMaximum', default_InitFreqMaximum);
addParameter(p, 'Device', default_Device);
addParameter(p, 'FPPrecision', default_FPPrecision);
parse(p, Q, K, Alpha, varargin{:});

% parameters
info.S = size(Q, 1);
info.T = size(Q, 2);
info.K = p.Results.K;
info.alpha = p.Results.alpha;
info.Tolerance = p.Results.Tolerance;
info.MaximumSteps = p.Results.MaximumSteps;
info.InitFreqType = p.Results.InitFreqType;
info.InitFreqMaximum = min(p.Results.InitFreqMaximum, FREQMAXIMUM);
info.Device = p.Results.Device;
info.FPPrecision = p.Results.FPPrecision;

Type = p.Results.FPPrecision;
N = p.Results.MaximumSteps;

%% pre-processing
S = info.S;
T = info.T*2;

% single precision ?
EPS = eps;
if (Type == 'single')
    Q = single(Q);
    Alpha = single(Alpha);
    EPS = 1e-6;
end

% mirror extension
Q_origin = Q; Q = zeros(S, T, Type);
T_half = ceil(info.T/2); % first half
Q(:,1:T_half) = Q_origin(:,T_half:-1:1);
Q(:,(T_half+1):(T_half+info.T)) = Q_origin;
Q(:,(T_half+info.T+1):end) = Q_origin(:,(info.T):-1:(T_half+1));

% variables initialization
Q_spec = fft(Q, [], 2);
T_spec = floor(T/2 + 1);
Q_spec = Q_spec(:, 1:T_spec); % positive half in the frequency domain

phi_n = zeros(S, K, Type) + EPS; % spatial modes
c_spec_n = zeros(T_spec, K, Type) + EPS; % time-evolution coefficients
omega_k = zeros(K, N, Type); % central frequencies
mode_k = zeros(S, T_spec, Type);

omega = (0:(T_spec-1))/T; % frequency list

% convert to GPU array
if (info.Device == 'gpu')
    Q_spec = gpuArray(Q_spec);
    
    phi_n = gpuArray(phi_n);
    c_spec_n = gpuArray(c_spec_n);
    omega_k = gpuArray(omega_k);
    mode_k = gpuArray(mode_k);
    
    omega = gpuArray(omega);
end

% central frequency initialization
switch info.InitFreqType
    case FREQRANDOM
        omega_k(:,1) = rand(K,1) * info.InitFreqMaximum;
    case FREQALLZERO
        omega_k(:,1) = 0;
    case FREQUNIFORM
        omega_k(:,1) = (0:1/(K-1):1) * info.InitFreqMaximum;
end

%% main loop
n = 1;
diff = info.Tolerance + EPS; % iteration difference initialization
residual_k = Q_spec - phi_n * c_spec_n.'; % residual function initialization

while (n <= N && diff > info.Tolerance)
    diff = 0;
    
    for k = 1:K
        % calculate residual matrix
        mode_k = phi_n(:, k) * c_spec_n(:, k).';
        residual_k = residual_k + mode_k;
        
        % update phi_k
        phi_n(:, k) = real(conj(residual_k) * c_spec_n(:, k));
        phi_n(:, k) = phi_n(:, k) / norm(phi_n(:, k), 'fro');
        
        % update c_k
        c_spec_n(:, k) = (residual_k.' * phi_n(:, k)) ./ ...
            (1 + 2 * Alpha * (omega - omega_k(k, n)).^2).';
        
        % update omega_k
        omega_k(k, n + 1) = omega * (conj(c_spec_n(:, k)) .* c_spec_n(:, k)) ...
            / norm(c_spec_n(:, k), 'fro')^2;
        
        mode_k_new = phi_n(:, k) * c_spec_n(:, k).';
        residual_k = residual_k - mode_k_new;
        
        % calculate iteration difference
        diff = diff + norm(mode_k_new - mode_k, 'fro') / norm(mode_k, 'fro');
    end
    
    diff_iter(n) = diff;
    disp(['iteration step: ', num2str(n), '    differences: ', num2str(diff)]);
    
    n = n + 1;
end

%% post-processing
info.Iteration.difference = diff_iter;

% convert to CPU array
if (info.Device == 'gpu')
    Q_spec = gather(Q_spec);
    
    phi_n = gather(phi_n);
    c_spec_n = gather(c_spec_n);
    omega_k = gather(omega_k);
    mode_k = gather(mode_k);
end

% central frequencies
info.Iteration.omega = omega_k(:,1:n);
omega = omega_k(:,n);

% reconstruct time-evolution coefficients
c_spec = c_spec_n;
c_spec((T_spec+1):T,:) = conj(c_spec_n((T_spec-1):-1:2,:));
c = ifft(c_spec);
c = c((T_half+1):(T_half+info.T),:);

% sort the rvmd modes according to central frequencies (from low to high)
[~, index] = sort(omega);

for k = 1:K
    mode.phi(:, k) = phi_n(:, index(k));
    c_sort(:, k) = c(:, index(k));
    omega_sort(k, :) = omega(index(k), :);
    mode.energy(k) = norm(c(:, index(k)))^2;
end

mode.c = c_sort;
mode.omega = omega_sort;

end
