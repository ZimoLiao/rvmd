function [mode, info, restart] = rvmd_new(Q, K, Alpha, varargin)
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
%   including parameter setups and options. The struct restart contains
%   data for restart computation.
%
%   [...] = rvmd(Q, K, Alpha, Name, Value) specifies the decomposition
%   options using one or more Name, Value pair arguments. Here is the list
%   of properties:
%   Properties:
%       'Restart' - input restart data to continue computation
%       'Weight' - weight vector for computing RVMD, default: identity
%       matrix
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
% addRequired(p, 'Q');
% validScalar = @(x) isnumeric(x) && (x>0);
% addRequired(p, 'K', validScalar);
% addRequired(p, 'alpha', validScalar);
addOptional(p, 'Q', 1);
addOptional(p, 'K', 1);
addOptional(p, 'alpha', 1);
addParameter(p, 'Restart', 0);
addParameter(p, 'Weight', 1);
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
info.weight = p.Results.Weight;
info.Tolerance = p.Results.Tolerance;
info.MaximumSteps = p.Results.MaximumSteps;
info.InitFreqType = p.Results.InitFreqType;
info.InitFreqMaximum = min(p.Results.InitFreqMaximum, FREQMAXIMUM);
info.Device = p.Results.Device;
info.FPPrecision = p.Results.FPPrecision;

Type = p.Results.FPPrecision;
N = p.Results.MaximumSteps;

% restart data
opt_restart = 0;
if (isfield(p.Results.Restart,'Q') ~= 0)
    restart = p.Results.Restart;
    opt_restart = 1;

    info.S = restart.S;
    info.T = restart.T;
    info.K = restart.K;
    info.alpha = restart.alpha;
    K = restart.K;
    Alpha = restart.alpha;
    info.weight = restart.weight;
    info.Iteration = restart.Iteration;

    Q = restart.Q;
end

%% pre-processing
S = info.S;
T = info.T*2;

% normalize weight vector
info.weight = info.weight/mean(info.weight);

% single precision ?
EPS = eps;
if (Type == 'single')
    Q = single(Q);
    Alpha = single(Alpha);
    EPS = 1e-6;
    info.weight = single(info.weight);
end

% weight vector specified ?
opt_weight = 0;
if (info.weight ~= 1 && length(info.weight) == S)
    opt_weight = 1;
    if (size(info.weight,1) == 1)
        info.weight = info.weight.'; % convert to column vector
    end
end

% mirror extension
restart.Q = Q; Q = zeros(S, T, Type);
T_half = ceil(info.T/2); % first half
Q(:,1:T_half) = restart.Q(:,T_half:-1:1);
Q(:,(T_half+1):(T_half+info.T)) = restart.Q;
Q(:,(T_half+info.T+1):end) = restart.Q(:,(info.T):-1:(T_half+1));

% variables initialization
Q_spec = fft(Q, [], 2);
T_spec = floor(T/2 + 1);
Q_spec = Q_spec(:, 1:T_spec); % positive half in the frequency domain

if (opt_restart)
    phi_n = restart.phi_n;
    c_spec_n = restart.c_spec_n;
    omega_k = zeros(K, N, Type);
    omega_k(:,1:info.Iteration.steps+1) = info.Iteration.omega;
else
    phi_n = zeros(S, K, Type) + EPS; % spatial modes
    c_spec_n = zeros(T_spec, K, Type) + EPS; % time-evolution coefficients
    omega_k = zeros(K, N, Type); % central frequencies
end
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
if (~opt_restart)
    switch info.InitFreqType
        case FREQRANDOM
            omega_k(:,1) = rand(K,1) * info.InitFreqMaximum;
        case FREQALLZERO
            omega_k(:,1) = 0;
        case FREQUNIFORM
            omega_k(:,1) = (0:1/(K-1):1) * info.InitFreqMaximum;
    end
end

%% main loop
if (opt_restart)
    n = info.Iteration.steps+1;
else
    n = 1;
end
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
        if (opt_weight)
            c_spec_n(:, k) = (residual_k.' * (phi_n(:, k).*info.weight)) ./ ...
                (1 + 2 * Alpha * (omega - omega_k(k, n)).^2).';
        else
            c_spec_n(:, k) = (residual_k.' * phi_n(:, k)) ./ ...
                (1 + 2 * Alpha * (omega - omega_k(k, n)).^2).';
        end

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
info.Iteration.steps = n-1;
info.Iteration.omega = omega_k(:,1:n);
omega = omega_k(:,n);

% reconstruct time-evolution coefficients
c_spec = c_spec_n;
c_spec((T_spec+1):T,:) = conj(c_spec_n((T_spec-1):-1:2,:));
c = ifft(c_spec);
c = c((T_half+1):(T_half+info.T),:);

% sort the rvmd modes according to central frequencies (from low to high)
[~, index] = sort(omega);

c_sort = c;
omega_sort = omega;
for k = 1:K
    mode.phi(:, k) = phi_n(:, index(k));
    c_sort(:, k) = c(:, index(k));
    omega_sort(k, :) = omega(index(k), :);
    mode.energy(k) = norm(c(:, index(k)))^2;
end

mode.c = c_sort;
mode.omega = omega_sort;

% restart data
restart.S = info.S;
restart.T = info.T;
restart.K = info.K;
restart.alpha = info.alpha;
restart.weight = info.weight;
restart.Iteration = info.Iteration;
restart.c_spec_n = c_spec_n;
restart.phi_n = phi_n;

end
