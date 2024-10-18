function [mode, info] = vmd(Q, K, Alpha, varargin)

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
addOptional(p, 'Q', 1);
addOptional(p, 'K', 1);
addOptional(p, 'alpha', 1);
addParameter(p, 'Tolerance', default_Tolerance);
addParameter(p, 'MaximumSteps', default_MaximumSteps);
addParameter(p, 'InitFreqType', default_InitFreqType);
addParameter(p, 'InitFreqMaximum', default_InitFreqMaximum);
addParameter(p, 'Device', default_Device);
addParameter(p, 'FPPrecision', default_FPPrecision);
parse(p, Q, K, Alpha, varargin{:});

% parameters
info.S = 1;
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
S = 1;
T = info.T*2;

% single precision ?
EPS = eps;
if (Type == 'single')
    Q = single(Q);
    Alpha = single(Alpha);
    EPS = 1e-6;
end

% mirror extension
restart.Q = Q; Q = zeros(S, T, Type);
T_half = ceil(info.T/2); % first half
Q(:,1:T_half) = restart.Q(:,T_half:-1:1);
Q(:,(T_half+1):(T_half+info.T)) = restart.Q;
Q(:,(T_half+info.T+1):end) = restart.Q(:,(info.T):-1:(T_half+1));

% variables initialization
Q_spec = fft(Q, [], 2);
T_spec = T;
Q_spec = fftshift(Q_spec);

c_spec_n = zeros(T_spec, K, Type) + EPS; % time-evolution coefficients
omega_k = zeros(K, N, Type); % central frequencies
mode_k = zeros(S, T_spec, Type);
omega = (-T_spec/2:(T_spec/2-1))/T; % frequency list

% convert to GPU array
if (info.Device == 'gpu')
    Q_spec = gpuArray(Q_spec);

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
residual_k = Q_spec - sum(c_spec_n,2).'; % residual function initialization

while (n <= N && diff > info.Tolerance)
    diff = 0;

    for k = 1:K
        % calculate residual matrix
        mode_k = c_spec_n(:, k).';
        residual_k = residual_k + mode_k;

        % update c_k
        c_spec_n(:, k) = (residual_k.') ./ ...
            (1 + 4 * Alpha * (abs(omega) - omega_k(k, n)).^2).';

        % update omega_k
        omega_k(k, n + 1) = abs(omega) * (conj(c_spec_n(:, k)) .* c_spec_n(:, k)) ...
            / norm(c_spec_n(:, k), 'fro')^2;

        mode_k_new = c_spec_n(:, k).';
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
c = ifft(ifftshift(c_spec,1));
c = c((T_half+1):(T_half+info.T),:);

% sort the rvmd modes according to central frequencies (from low to high)
[~, index] = sort(omega);

c_sort = c;
omega_sort = omega;
for k = 1:K
    c_sort(:, k) = c(:, index(k));
    omega_sort(k, :) = omega(index(k), :);
    mode.energy(k) = norm(c(:, index(k)))^2;
end

mode.c = c_sort;
mode.omega = omega_sort;

end